
#' Structure of the S4 class "growthSimulation"
#'
#' Structure of the S4 class \code{growthSimulation} as framework for the grwoth environment and container for agent-based flux balance analysis.
#' @import Matrix Rcpp
#' @exportClass growthSimulation
#'
#' @slot n_rounds integer for the number of simulation rounds that were already performed for this \code{growthSimulation} object.
#' @slot deltaTime double for length of each time step for the simulation in hours.
#' @slot diffusionNIter integer. Number of diffusion steps per simulation round.
#' @slot rMotion double. Maximum x- and y distance a cell can travel by means of random movement per simulation round.
#' @slot models list for \code{Organism} objects to represent the different strains in the simulation.
#' @slot history list with recordings of simulation status information at each simulation round.
#' @slot universePolygon Matrix specifiying the corners of the polygon that defined the growth environment boundaries. 2-dimensional: x and y.
#' @slot environ Object of S4-class \code{grwothEnvironment}, that specifies the environment mesh layout, compounds, and their concentrations.
setClass("growthSimulation",

         slots = c(

           n_rounds = "numeric",

           # simulation slots
           deltaTime      = "numeric",
           diffusionNIter = "numeric", # diffusion steps per simulation round
           rMotion        = "numeric",
           models         = "list",
           history        = "list",

           # ind. cell info
           cellDT = "data.table",

           # Environment slots
           universePolygon = "matrix",
           environ         = "growthEnvironment"

         )
)



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#  Initialize Simulation Object #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Initialise a growth simulation.
#'
#' @description Method to initialise a \code{growthSimulation} object
#'
#' @export
#' @rdname growthSimulation-Initiation
#' @param universePolygon A two column matrix specifiying the x and y coordinates of the polygon, that describes the growth environment boundaries.
#' @param gridFieldSize double. Distance between neighboring environments 3D mesh field elements (rhombic dodecahedrons) in µm.
#' @param gridFieldLayers integer. z-dimension (height) as the number of layers of field elements.
#' @param deltaTime double specifying the length of each time step for the simulation in hours.
#' @param diffusion.alpha double. This number should be between 0 and 1 and specifies the rate of naive diffusion, but specifying the fraction of compound in a field that is equally distributed to the neighboring fields. Default: 12/13.
#' @param diffusion.niter integer. Number of diffusion steps per simulation round.
#' @param pFBAcoeff integer. pFBA coefficient. Default: 1e-6
#' @param rMotion double. Maximum x- and y distance a cell can travel by means of random movement per minute. Default: 1/6 µm
#'
#' @return Object of class \code{growthSimulation}.
#'
#' @examples
#' # Contruction a square environment of dimensions 300µm x 300µm x 5µm
#' sim <- init.simulation(cbind(c(-150, -150, 150, 150), c(-150, 150, 150, -150)),
#'                        gridFieldSize = 1.75, gridFieldLayers = 5)
setGeneric(name="init.simulation",
           def=function(universePolygon,
                        gridFieldSize = 1,
                        gridFieldLayers = 3,
                        deltaTime = 1/6,
                        diffusion.alpha = 6/7,
                        diffusion.niter = 60,
                        pFBAcoeff = 1e-6,
                        rMotion = 1/6, ...)
           {
             standardGeneric("init.simulation")
           }
)

setMethod(f          = "init.simulation",
          signature  = signature(universePolygon = "matrix"),
          definition = function(universePolygon,
                                gridFieldSize   = 1,
                                gridFieldLayers = 3,
                                deltaTime       = 1/6,
                                diffusion.alpha = 12/13, # each grid field has 12 neighbors. Plus itself: 13. 12/13 means that the concentration in one field is equally distributed among all its neigbors and itself uin each diffusion iteration step
                                diffusion.niter = floor(deltaTime * 360 * 1/gridFieldSize)*2,
                                rMotion         = deltaTime, ...) {


            # init growth environment
            environ <- new("growthEnvironment",
                           polygon.coords  = universePolygon,
                           field.size      = gridFieldSize,
                           diffusion.alpha = diffusion.alpha,
                           field.layers    = gridFieldLayers)

            # construct empty cell info table
            cellDT <- data.table(cell = numeric(),
                                 type = character(),
                                 x = numeric(), y = numeric(),
                                 x.vel = numeric(), y.vel = numeric(),
                                 mass = numeric(),
                                 size = numeric(),
                                 parent = numeric())

            # init actual growth simulation object
            simobj <- new("growthSimulation",
                          n_rounds = 0,
                          deltaTime = deltaTime,
                          diffusionNIter = diffusion.niter,
                          rMotion = rMotion,
                          models = list(),
                          history = list(),
                          cellDT = cellDT,
                          universePolygon = universePolygon,
                          environ = environ
            )

            # setting the solver to cplex
            sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <<- 1

            return(simobj)
          }
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#  Add organism cells to simulation object  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Add a model/organism to simulation
#'
#' @description Adds an organism an its genome-scale metabolic network model to the growth simulation object.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param model The organisms metabolic model of S4-type \code{sybil::modelorg}
#' @param name Character for the name of the model, that will also be used for plotting.
#' @param ncells integer. Number of initial cells to be added to the growth simulation.
#' @param coords (optional) A two column numerical matrix specifying the coordinates (1st column x, 2nd column y) of the initial cells. If provided, the number of rows should be equal to \code{ncells}. Default: NULL
#' @param distribution.method If \code{coords} is \code{NULL}, this parameter specifies the distribution method for initial cells. Default: "random_centroid"
#' @param distribution.center Numeric vector of length 2, which specifies the coordinates of the centre for the \code{distribution.method}.
#' @param distribution.radius double. Spcifies the radius (in µm) in which inital cells are distributed.
#' @param cellDiameter double. Diameter in µm of initial cells.
#' @param cellMassInit double. Mass in pg of initial cells. Default is 0.28 pg
#' @param cellMassAtDivision double. Cell mass at which a cell devides in two daughter cells. Default: 0.55 pg
#' @param cellShape character. Shape of cells. Currently only "coccus" is supported.
#' @param vmax double. Maximum velocity of a cell in µm per minute.
#' @param scavengeDist double. Distance in µm a cell can scavenge nutrients from its surrounding.
#' @param rm.deadends If TRUE, dead-end metabolites and reactions are removed from the \code{model}, which reduces the computation time for FBA, but has otherwise no effect on the flux distribution solutions.
#'
#' @return Object of class \code{growthSimulation}.
#'
#' @references
#'  \url{https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100008} \cr
#'  \url{http://book.bionumbers.org/how-big-is-an-e-coli-cell-and-what-is-its-mass/}
#'
#' @examples
#' # add two bacterial models (Eubacterium rectale, Bifidobacterium longum)
#' # to the environment; each with 15 initial cells
#' models <- list()
#' models[['eure']] <- readRDS(system.file("extdata", "eure.RDS", package="EcoAgents"))
#' models[['bilo']] <- readRDS(system.file("extdata", "bilo.RDS", package="EcoAgents"))
#'
#' sim <- init.simulation(cbind(c(-150, -150, 150, 150), c(-150, 150, 150, -150)),
#'                        gridFieldSize = 1.75, gridFieldLayers = 5)
#'
#' sim <- add.organism(sim, model = models[["eure"]], name = "E. rectale",
#'                     ncells = 15, distribution.radius = 30)
#' sim <- add.organism(sim, model = models[["bilo"]], name = "B. longum",
#'                     ncells = 15, distribution.radius = 30)
#'
#' plot.cells(sim, xlim = c(-50,50), ylim= c(-50,50))
#'
#' @export
setGeneric(name="add.organism",
           def=function(object,
                        model,
                        name,
                        ncells,
                        coords = NULL,
                        distribution.method = "random_centroid",
                        distribution.center = NULL,
                        distribution.radius = NULL,
                        cellDiameter = (3 * 1 / (4 * pi))^(1/3) * 2, # diameter of sphere with 1 µm^3
                        cellMassInit = 0.28,
                        cellMassAtDivision = 0.55,
                        cellShape = "coccus",
                        vmax = 1,
                        scavengeDist = 5,
                        rm.deadends = T, ...
           )
           {
             standardGeneric("add.organism")
           }
)

setMethod(f          = "add.organism",
          signature  = signature(object = "growthSimulation",
                                 model  = "modelorg",
                                 name   = "character",
                                 ncells = "numeric"),
          definition = function(object,
                                model,
                                name,
                                ncells,
                                coords = NULL,
                                distribution.method = "random_centroid",
                                distribution.center = NULL,
                                distribution.radius = NULL,
                                cellDiameter = (3 * 1 / (4 * pi))^(1/3) * 2, # diameter of sphere with 1 µm^3
                                cellMassInit = 0.28,
                                cellMassAtDivision = 0.56,
                                cellShape = "coccus",
                                vmax = 1,
                                scavengeDist = cellDiameter*2.5,
                                rm.deadends = T,
                                ...) {

            # init new organism object
            object@models[[name]] <- new("Organism",
                                         cellDiameter = cellDiameter,
                                         cellMassInit = cellMassInit,
                                         cellMassAtDivision = cellMassAtDivision,
                                         vmax = vmax,
                                         scavengeDist = scavengeDist,
                                         mod = model,
                                         rm.deadends = rm.deadends)

            #if not provided -> assign cell positions:
            if(is.null(coords)) {
              if(!(distribution.method %in% c("random","random_centroid")))
                stop("Distribution method not supported. Choose one of \"random\",\"random_centroid\".")
              universePG <- Polygon(object@universePolygon)

              if(distribution.method == "random") {
                coords  <- spsample(universePG, ncells, type = "random")@coords
              }
              if(distribution.method == "random_centroid") {

                PGcent <- distribution.center
                if(is.null(PGcent))
                  PGcent <- SpatialPoints(matrix(universePG@labpt,ncol = 2))

                maxRadius <- distribution.radius
                if(is.null(maxRadius)) {
                  PGline <- Line(universePG@coords)
                  PGlines <- SpatialLines(list(Lines(list(PGline), ID = "one")))
                  maxRadius <- gDistance(PGlines, PGcent) / 3
                }

                coords <- matrix(0, ncol = 2, nrow = ncells)

                r_angle  <- runif(ncells, max = 2*pi)
                r_radius <- runif(ncells, max = maxRadius)

                coords[,1] <- cos(r_angle) * r_radius + PGcent@coords[1,1]
                coords[,2] <- sin(r_angle) * r_radius + PGcent@coords[1,2]
              }

            }

            # add cells to cell info table
            newcells <- data.table(cell = (1:ncells)+nrow(object@cellDT),
                                   type = name,
                                   x = coords[,1], y = coords[,2],
                                   x.vel = 0, y.vel = 0,
                                   mass = cellMassInit,
                                   size = cellDiameter,
                                   parent = NA_integer_)
            object@cellDT <- rbind(object@cellDT, newcells)

            # add compounds, whose exchange reactions hab a LB < 0, to environment
            mod <- object@models[[name]]@mod

            dt_exr <- data.table(id = mod@react_id,
                                 lb = mod@lowbnd,
                                 name = mod@react_name)
            dt_exr <- dt_exr[grepl("^EX_", id)]
            dt_exr[, id := gsub("^EX_","",id)]
            dt_exr <- dt_exr[id != "cpd11416_c0"]
            dt_exr[, name := gsub("-e0-e0 Exchange","",name)]
            dt_exr[, name := gsub("-e0 Exchange","",name)]
            dt_exr[, name := gsub(" Exchange","",name)]

            object <- add.compounds(object,
                                    compounds = dt_exr$id,
                                    concentrations = rep(0, nrow(dt_exr)),
                                    compound.names = dt_exr$name,
                                    is.constant = rep(FALSE, nrow(dt_exr)))

            return(object)
          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# show method for small summary #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Print a short summary of a growth simulation
#'
#' @description Displays a few numbers to describe the current status of a growth simulation.
#'
#' @param object S4-object of type \code{growthSimulation}.
#'
#' @export
setMethod(f          = "show",
          signature  = signature(object = "growthSimulation"),
          definition = function(object) {
            cat("Cell community growth simulation\n")
            cat("Passed simulation time:\t",object@n_rounds*object@deltaTime," hours (",
                object@n_rounds," rounds)\n", sep= '')
            cat("\n")

            cat("Cell counts:\n")
            print(object@cellDT[, .(cells = .N, `mass in pg` = sum(mass)), by = "type"])
            cat("\n")

            # Environment
            cat("Cell growth environment\n")
            cat("    Universe dimensions (µm):\t\t",
                abs(round(min(object@environ@field.pts@coords[,1])-max(object@environ@field.pts@coords[,1]), digits = 2))," x ",
                abs(round(min(object@environ@field.pts@coords[,2])-max(object@environ@field.pts@coords[,2]), digits = 2))," x ",
                abs(round(min(object@environ@field.pts@coords[,3])-max(object@environ@field.pts@coords[,3]), digits = 2)),"\n")
            cat("    Universe volume (µm^3):\t\t", round(object@environ@fieldSize * object@environ@nfields, digits = 2) , "\n")
            cat("    Number of rhombic dodecahedrons:\t",object@environ@nfields,"\n")
            cat("    Number of variable compounds:\t",length(object@environ@compounds),"\n")

          }
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot cell distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Plot spatial distribution of cells
#'
#' @description Plots the distribution of cells in a growth simulation.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param xlim Numeric vector of length 2, specifying the x-range to be displayed.
#' @param ylim Numeric vector of length 2, specifying the y-range to be displayed.
#' @param iter Positive integer number of the simulation step/iteration to plot the cell distribution. If past simulation are displayed, the cell postitions needed to be recorded when running the simulation before (see \code{link{run.simulation}}). If NULL, current distribution is displayed.
#'
#' @return A ggplot object.
#'
#' @export
setGeneric(name="plot.cells",
           def=function(object, xlim = NULL, ylim = NULL, ...)
           {
             standardGeneric("plot.cells")
           }
)

setMethod(f = "plot.cells",
          signature = signature(object = "growthSimulation"),
          definition = function(object, xlim = NULL, ylim = NULL, iter = NULL) {

            # If no limits are defined use polygon universe limits
            if(is.null(xlim)) {
              xlim <- c(min(object@universePolygon[,1]),
                        max(object@universePolygon[,1]))
            }
            if(is.null(ylim)) {
              ylim <- c(min(object@universePolygon[,2]),
                        max(object@universePolygon[,2]))
            }

            cellposDT <- data.table()
            # current
            if(is.null(iter)) {
              cellposDT <- object@cellDT
            }

            # from history
            if(!is.null(iter)) {
              if(iter > object@n_rounds) {
                warning(paste0("Simulation did not run ",iter," iterations yet. Displaying results for last iteration (",object@n_rounds,")"))
                iter <- object@n_rounds
              }

              cellposDT <- object@history[[iter]]$cells
            }

            # get expansion factor so scale bar is not to close to panel margins
            x_exp_fac <- (xlim[2]-xlim[1]) * 0.05
            y_exp_fac <- (ylim[2]-ylim[1]) * 0.05

            # get scale bar length to span approx 10% of x-axis
            x_span <- xlim[2]-xlim[1]
            x_magn <- floor(log(x_span, base = 10))
            bar_wd <- 10^x_magn / 10
            bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd * 2.5, bar_wd) # e.g. 10 -> 25
            bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd / 2.5 * 5, bar_wd) # e.g. 10 -> 50

            p <- ggplot(cellposDT, aes(x0 = x,y0 = y, r = size/2, fill = type, col = type)) +
              geom_circle(alpha = 0.3, show.legend = T, key_glyph = draw_key_point) +
              coord_equal(xlim = xlim, ylim = ylim) +
              geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac, xend = xlim[2]-x_exp_fac,
                               y = ylim[1]+y_exp_fac, yend = ylim[1]+y_exp_fac), color = "white") +
              annotate("text", x = xlim[2]-bar_wd/2-x_exp_fac, y = ylim[1]+y_exp_fac, label = paste0(bar_wd," µm"),
                       color = "white", hjust = 0.5, vjust = -1, size = 2.5) +
              scale_color_brewer(palette = "Set1") +
              scale_y_continuous(sec.axis = sec_axis(~ .)) + scale_x_continuous(sec.axis = sec_axis(~ .)) +
              theme_bw() +
              theme(plot.background = element_blank(),
                    axis.line.x.top = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.x.bottom = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.y.left = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.y.right = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.line = element_blank(),
                    axis.text = element_blank(),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_rect(fill = "black"),
                    legend.background = element_blank(),
                    legend.text = element_text(color = "black", face = "italic"),
                    legend.box.background = element_blank(),
                    legend.key = element_blank(),
                    legend.position = "bottom", legend.direction="vertical"
              ) +
              guides(color = guide_legend(title = "Organism",
                                          override.aes = list(alpha = 1, size = 5)),
                     fill  = guide_legend(title = "Organism"))

            return(p)

          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot growth curves      #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Plot growth curves
#'
#' @description Plot population growth over time
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param tlim Numeric vector of length 2, specifying the x-range (Time) to be displayed.
#'
#' @return A ggplot object.
#'
#' @export
setGeneric(name="plot.growth",
           def=function(object, tlim = NULL, ...)
           {
             standardGeneric("plot.growth")
           }
)


setMethod(f = "plot.growth",
          signature = signature(object = "growthSimulation"),
          definition = function(object, tlim = NULL) {

            # Sanity checks
            if(object@n_rounds < 2)
              stop("Simulation did not yet run for at least 2 iterations. No growth can be plotted.")
            if(is.null(object@history[[1]]$cells))
              stop("No cell/population data has been recorded.")

            dt_cellsg <- rbindlist(lapply(object@history, function(x) x$cells), idcol = "iteration")
            dt_cellsg <- dt_cellsg[, .(mass = sum(mass)), by = .(iteration, type)]
            dt_cellsg[, time := iteration * object@deltaTime]

            p_growth <- ggplot(dt_cellsg, aes(time, mass, col = type, group = type)) +
              geom_line() +
              #geom_point(cex = 0.5) +
              scale_color_brewer(palette = "Set1") +
              labs(x = "Time [hr]", y = "Mass [pg]",col = "Organism") +
              scale_x_continuous(expand = c(0,0)) +
              theme_bw() +
              theme(panel.grid = element_blank(),
                    panel.border = element_blank(),
                    axis.title = element_text(color = "black"),
                    axis.text = element_text(color = "black"),
                    axis.line = element_line(color = "black"),
                    legend.text = element_text(color = "black", face = "italic"))

            return(p_growth)
          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot compound time curves #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Plot compound time curves
#'
#' @description Plot concentrations of metabolites over time
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param compounds Vector of compound IDs whose concentration is plotted. If no IDs are provided,
#' the top variable (SD) compounds are plotted (max. 8 compounds).
#' @param tlim Numeric vector of length 2, specifying the x-range (Time) to be displayed.
#'
#' @return A ggplot object.
#'
#' @export
setGeneric(name="plot.compounds",
           def=function(object, compounds = NULL, tlim = NULL, ...)
           {
             standardGeneric("plot.compounds")
           }
)


setMethod(f = "plot.compounds",
          signature = signature(object = "growthSimulation"),
          definition = function(object, compounds = NULL, tlim = NULL) {

            # Sanity checks
            if(object@n_rounds < 2)
              stop("Simulation did not yet run for at least 2 iterations. No growth can be plotted.")
            if(is.null(object@history[[1]]$global_compounds))
              stop("No compound concentration data has been recorded.")

            dt_cpdg <- rbindlist(lapply(object@history, function(x) x$global_compounds), idcol = "iteration")
            dt_cpdg[, time := iteration * object@deltaTime]


            if(is.null(compounds)) {
              # identify most variable compounds
              rel_cpds <- dt_cpdg[, sd(global_concentration), by = cpd.id][V1 > 1e-4][order(-V1), cpd.id]
              rel_cpds <- rel_cpds[1:min(c(8, length(rel_cpds)))]
            } else {
              # check if slected compounds are part of models/simulation
              compounds <- compounds[compounds %in% object@environ@compounds]
              if(length(compounds) == 0)
                stop("Selected compounds not part of the simulation.")
              rel_cpds <- compounds
            }

            # filter to relevant cpds
            dt_cpdg <- dt_cpdg[cpd.id %in% rel_cpds]

            p_cpdth <- ggplot(dt_cpdg, aes(time, global_concentration, col = cpd.name, group = cpd.name)) +
              geom_line() +
              #geom_point(cex = 0.5) +
              scale_color_brewer(palette = "Set1") +
              labs(x = "Time [hr]", y = "Concentration [mM]",col = "Compound") +
              scale_x_continuous(expand = c(0,0)) +
              theme_bw() +
              theme(panel.grid = element_blank(),
                    panel.border = element_blank(),
                    axis.title = element_text(color = "black"),
                    axis.text = element_text(color = "black"),
                    axis.line = element_line(color = "black"),
                    legend.text = element_text(color = "black"))

            return(p_cpdth)

          }
)

#' @title Add compounds to the growth environment
#'
#' @description The function can be used to add substances to the growth environment by providing concentrations in mM.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param compounds Character vector with the compound IDs of substances to add to the environment. Compound IDs should correspond to the models' exchange reactions IDs ("EX_[cpdid]"), without the "EX_" prefix.
#' @param concentrations Numeric vector with the concentrations of the compounds from \code{compounds} in the same order. Values in mM.
#' @param compound.names Character vector with the compound names.
#' @param is.constant Logical vector that indicates if the compound should remain constant over time despite of potential uptake or production of cells.
#'
#' @details Compound concentration are equally distributed across the whole growth environment.
#' If compound is already present, old and new concentrations are added. More options are planned.
#'
#' @return Return a S4-object of type \code{growthSimulation}.
#'
#' @examples
#' sim <- init.simulation(cbind(c(-100, -100, 100, 100), c(-100, 100, 100, -100)),
#'                        gridFieldSize = 2, gridFieldLayers = 5)
#'
#' sim <- add.compounds(sim,
#'                      compounds = c("cpd00027_e0","cpd00029_e0","cpd00047_e0",
#'                                    "cpd00159_e0","cpd00211_e0"),
#'                      concentrations = c(50,0,0,0,0),
#'                      compound.names = c("D-Glucose","Acetate","Formate",
#'                                         "L-Lactate","Butyrate"),
#'                      is.constant = c(F, F, F, F, F))
#'
#' @export
setGeneric(name="add.compounds",
           def=function(object, compounds, concentrations,
                        compound.names = NULL, is.constant = NULL)
           {
             standardGeneric("add.compounds")
           }
)

setMethod(f          = "add.compounds",
          signature  = signature(object         = "growthSimulation",
                                 compounds      = "character",
                                 concentrations = "numeric",
                                 compound.names = "character",
                                 is.constant    = "logical"),
          definition = function(object, compounds, concentrations,
                                compound.names = NULL,
                                is.constant = NULL) {

            if(length(compounds) != length(concentrations))
              stop("Lengths of 'compounds' and 'concentrations' differ.")

            # make compounds non-constant if nothing else is specified
            if(is.null(is.constant)) {
              is.constant <- rep(F, length(compounds))
            }

            # is.constant should be NULL (see above) or of lenfth 1 or n (nr. of compounds)
            if(length(is.constant) != 1 & length(is.constant) != length(compounds)) {
              stop("Length of 'is.constant' is not 1 or the same length als 'compounds'")
            }

            # extent is.constant vector if only one value
            if(length(is.constant) == 1) {
              is.constant <- rep(is.constant, length(compounds))
            }

            # if no names are supplied -> NA
            if(is.null(compound.names)) {
              compound.names <- rep(NA_character_, length(compounds))
            }

            # if length of names is unequal the length of ids -> something's odd
            if(length(compound.names) != length(compounds)) {
              stop("Length of 'compound.names' is not the same length als 'compounds'")
            }

            for(i in 1:length(compounds)) {
              # Metabolite is already present in environment
              if(compounds[i] %in% object@environ@compounds) {
                c_ind <- which(object@environ@compounds == compounds[i])
                object@environ@concentrations[, c_ind] <- object@environ@concentrations[, c_ind] + concentrations[i]
                object@environ@conc.isConstant[c_ind]  <- is.constant[i]
                if(!is.na(compound.names[i]))
                  object@environ@compound.names[c_ind] <- compound.names[i]

              # Metabolite is new
              } else {
                object@environ@compounds       <- c(object@environ@compounds, compounds[i])
                object@environ@concentrations  <- cbind(object@environ@concentrations,
                                                       matrix(concentrations[i], ncol = 1, nrow = object@environ@nfields))
                object@environ@conc.isConstant <- c(object@environ@conc.isConstant, is.constant[i])
                if(!is.na(compound.names[i])) {
                  if(compound.names[i] %in% object@environ@compound.names) {
                    warning(paste0("Compound with name '",compound.names[i],"' already exists. Replacing with '",compounds[i],"'."))
                  }
                  object@environ@compound.names <- c(object@environ@compound.names, compound.names[i])
                } else {
                  object@environ@compound.names <- c(object@environ@compound.names, compounds[i])
                }
              }
            }

            return(object)
          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot compound distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Plot spatial distribution of compounds
#'
#' @description Plots the spatial distribution of compound concentrations in a heatmap.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param compounds Character string of compound ids to plot
#' @param compound.names Character string of compound names that should be displayed as facet header instead of compound names from the simulation object.
#' @param xlim Numeric vector of length 2 specifying the limits (left and right) for x axis; i.e. the horizontal dimension.
#' @param ylim Numeric vector of length 2 specifying the limits (top and bottom) for y axis; i.e. the vertical dimension.
#' @param iter Positive integer number of the simulation step/iteration to plot the distribution. Works only if the respective metabolite
#' concentrations were prior recorded (see \code{link{run.simulation}}). If \code{NULL}, current distribution of compound concentrations.
#'
#' @return A ggplot object.
#'
#' @export
setGeneric(name="plot.environment",
           def=function(object, compounds, compound.names = NULL, xlim = NULL, ylim = NULL, ...)
           {
             standardGeneric("plot.environment")
           }
)

setMethod(f = "plot.environment",
          signature = signature(object = "growthSimulation",
                                compounds = "character"),
          definition = function(object, compounds, compound.names = NULL, xlim = NULL, ylim = NULL, iter = NULL) {

            # Argument sanity check
            if(!all(compounds %in% object@environ@compounds)) {
              compounds <- compounds[compounds %in% object@environ@compounds]

              if(length(compounds) == 0)
                stop("Selected compound(s) not in list of variable environment compounds.")

              warning("Not all selected compounds are in the list of variable environment compounds. Continuing with the rest...")
            }

            # If no limits are defined use polygon universe limits
            if(is.null(xlim)) {
              xlim <- c(min(object@universePolygon[,1]),
                        max(object@universePolygon[,1]))
            }
            if(is.null(ylim)) {
              ylim <- c(min(object@universePolygon[,2]),
                        max(object@universePolygon[,2]))
            }

            # If no compound names are defined, use compound origunal compound names instead
            cpd_nameDT <- data.table(Compound = compounds,
                                     Compound.name = object@environ@compound.names[match(compounds, object@environ@compounds)])
            if(!is.null(compound.names)) {
              if(length(compound.names) != length(compounds)){
                warning("Number of compound names not equal to number of compound IDs. Using original compound names instead.")
              } else if(any(duplicated(compound.names))) {
                warning("Duplicate compound names. Using original compound names instead.")
              } else {
                cpd_nameDT$Compound.name <- compound.names
              }
            }

            if(is.null(iter)) {
              envDT <- data.table(object@environ@concentrations)
            } else {
              # TODO: When a former simulation step is selected
            }

            names(envDT) <- object@environ@compounds
            envDT <- envDT[, ..compounds]

            envDT <- cbind(object@environ@field.pts@coords, envDT)
            envDT <- envDT[z == 0] # baseplane
            envDT <- melt(envDT, id.vars = c("x","y"), variable.name = "Compound", value.name = "mM")
            envDT <- merge(envDT, cpd_nameDT)

            # get expansion factor so scale bar is not to close to panel margins
            x_exp_fac <- (xlim[2]-xlim[1]) * 0.05
            y_exp_fac <- (ylim[2]-ylim[1]) * 0.05

            # get scale bar length to span approx 10% of x-axis
            x_span <- xlim[2]-xlim[1]
            x_magn <- floor(log(x_span, base = 10))
            bar_wd <- 10^x_magn / 10
            bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd * 2.5, bar_wd) # e.g. 10 -> 25
            bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd / 2.5 * 5, bar_wd) # e.g. 10 -> 50

            p <- ggplot(envDT, aes(x,y, fill = mM)) +
              #geom_hex(aes(colour = mM), stat = "identity") +
              geom_bin2d(aes(colour = mM), stat = "identity") +
              coord_equal(xlim = xlim, ylim = ylim) +
              geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac, xend = xlim[2]-x_exp_fac,
                               y = ylim[1]+y_exp_fac, yend = ylim[1]+y_exp_fac), color = "white") +
              annotate("text", x = xlim[2]-bar_wd/2-x_exp_fac, y = ylim[1]+y_exp_fac, label = paste0(bar_wd," µm"),
                       color = "white", hjust = 0.5, vjust = -1, size = 2.5) +
              theme_bw() +
              scale_fill_viridis_c(guide = guide_colourbar(direction = "vertical")) + scale_color_viridis_c() + facet_wrap("Compound.name") +
              scale_y_continuous(sec.axis = sec_axis(~ .)) + scale_x_continuous(sec.axis = sec_axis(~ .)) +
              theme(legend.position = "right",
                    axis.line.x.top = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.x.bottom = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.y.left = element_line(color = "white", size = 1.5, lineend = "round"),
                    axis.line.y.right = element_line(color = "white", size = 1.5, lineend = "round"),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                    panel.background = element_blank(), axis.line = element_blank(),
                    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                    strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))


            return(p)
            #return(plot.environment(object@environ, compounds, compound.names = compound.names, xlim = xlim, ylim = ylim))

          }
)
