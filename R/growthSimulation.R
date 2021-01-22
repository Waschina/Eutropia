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
setGeneric(name="init.simulation",
           def=function(universePolygon,
                        gridHexSize = 1,
                        gridHexHeight = 10,
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
                                gridHexSize     = 1,
                                gridHexHeight   = 10,
                                deltaTime       = 1/6,
                                diffusion.alpha = 6/7,
                                diffusion.niter = floor(deltaTime * 360 * 1/gridHexSize),
                                rMotion         = deltaTime, ...) {


            # init growth environment
            environ <- new("growthEnvironment",
                           polygon.coords = universePolygon,
                           hexagon.size = gridHexSize,
                           diffusion.alpha = diffusion.alpha,
                           hex.height = gridHexHeight)

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
setGeneric(name="add.organism",
           def=function(object,
                        model,
                        name,
                        ncells,
                        coords = NA,
                        distribution.method = "random_centroid",
                        distribution.center = NULL,
                        distribution.radius = NULL,
                        cellDiameter = 1,
                        cellMassInit = 0.28,
                        cellMassAtDivision = 0.55,
                        cellShape = "coccus",
                        vmax = 1,
                        scavengeDist = 1,
                        pFBAcoeff = 1e-6,
                        rm.deadends = T, ...
                        )
           {
             standardGeneric("add.organism")
           }
)

#'
#'
#' Defaults from: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100008 and http://book.bionumbers.org/how-big-is-an-e-coli-cell-and-what-is-its-mass/
setMethod(f          = "add.organism",
          signature  = signature(object = "growthSimulation",
                                 model  = "modelorg",
                                 name   = "character",
                                 ncells = "numeric"),
          definition = function(object,
                                model,
                                name,
                                ncells,
                                coords = NA,
                                distribution.method = "random_centroid",
                                distribution.center = NULL,
                                distribution.radius = NULL,
                                cellDiameter = (3 * 1 / (4 * pi))^(1/3) * 2, # dimeter of sphere with 1 µm^3
                                cellMassInit = 0.28,
                                cellMassAtDivision = 0.56,
                                cellShape = "coccus",
                                vmax = 1,
                                scavengeDist = cellDiameter/2,
                                pFBAcoeff = 1e-6,
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
                                         pFBAcoeff = pFBAcoeff,
                                         rm.deadends = rm.deadends)

            #if not provided -> assign cell positions:
            if(is.na(coords)) {
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


            return(object)
          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# show method for small summary #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
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

            show(object@environ)
          }
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot cell distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
setGeneric(name="plot.cells",
           def=function(object, compounds, xlim = NULL, ylim = NULL, ...)
           {
             standardGeneric("plot.cells")
           }
)

#' Plot spatial distribution of cells
#'
#' @param object An object of class \code{growthEnvrionment} or \code{growthSimulation}
#' @param xlim Limits for the x axis
#' @param ylim Limits for the y axis
#' @param iter Positive integer number of the simulation step/iteration to plot the distribution. Works only if \code{object} is of
#' class \code{growthSimulation} and if the respective metabolite concentrations were recorded (sie \code{link{run.simulation}}).
setMethod(f = "plot.cells",
          signature = signature(object = "growthSimulation"),
          definition = function(object, xlim = NULL, ylim = NULL, iter = NULL) {

            cellposDT <- data.table()
            # current
            if(is.null(iter)) {
              cellposDT <- object@cellDT
            }

            # from history
            if(!is.null(iter)) {
              cellposDT <- object@history[[iter]]$cells
            }

            p <- ggplot(cellposDT, aes(x0 = x,y0 = y, r = size/2, fill = type, col = type)) +
              geom_circle(alpha = 0.3, show.legend = T) +
              coord_equal(xlim = xlim, ylim = ylim) +
              scale_color_brewer(palette = "Set1") +
              theme_bw() +
              theme(plot.background = element_rect(fill = "black", color = NA),
                    axis.title.x = element_text(color = "white"),
                    axis.title.y = element_text(color = "white"),
                    axis.ticks = element_line(color = "white"),
                    axis.text = element_text(color = "white"),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    legend.background = element_blank(),
                    legend.text = element_text(color = "white"),
                    legend.box.background = element_blank(),
                    legend.key = element_blank(),
                    legend.title = element_text(color = "white"),
                    panel.background = element_blank(),
                    plot.margin=unit(c(0,0,0,0), "mm")
              )

            return(p)

          }
)

