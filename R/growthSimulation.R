
#' Structure of the S4 class "growthSimulation"
#'
#' @description Structure of the S4 class \link{growthSimulation} as framework
#' for the growth environment and container for agent-based flux balance
#' analysis.
#'
#' @aliases growthSimulation
#'
#' @exportClass growthSimulation
#'
#' @slot n_rounds integer for the number of simulation rounds that were already
#' performed for this \link{growthSimulation} object.
#' @slot deltaTime double for length of each time step for the simulation in
#' hours.
#' @slot rMotion double. Maximum x- and y distance a cell can travel by means of
#' random movement per simulation round. Unit: \eqn{\mu}m per minute.
#' @slot models List for \link{Organism} objects to represent the different
#' strains in the simulation.
#' @slot history list with recordings of simulation status information at each
#' simulation round.
#' @slot cellDT Data table with individual cell information (e.g. size,
#' position, velocity, type)
#' @slot universePolygon Matrix specifying the corners of the polygon that
#' defines the growth environment boundaries. 2-dimensional: x and y.
#' @slot environ Object of S4-class \link{growthEnvironment}, that specifies the
#' environment mesh layout, compounds, and their concentrations.
#' @slot recordDir Directory name, in which intermediate compound concentrations
#' are recorded if turned on in \link{run_simulation}. Files in this directory
#' are meant as internal resource and not for direct analysis outside of this
#' package.
#' @slot rcdt data.table that stores the last reduced cost values of each cell's
#' exchange reactions. This information could indicate growth limiting
#' compounds.
setClass("growthSimulation",

         slots = c(

           n_rounds = "numeric",

           # simulation slots
           deltaTime      = "numeric",
           rMotion        = "numeric",
           models         = "list",
           history        = "list",

           # ind. cell info
           cellDT = "data.table",
           rcdt   = "data.table",

           # Environment slots
           universePolygon = "matrix",
           environ         = "growthEnvironment",

           # file for compound recordings
           recordDir = "character"
         )
)

# small helper for sanity checks
is.growthSimulation <- function(x) inherits(x, "growthSimulation")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#  Initialize Simulation Object #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' Initialize a growth simulation.
#'
#' Method to initialize a \link{growthSimulation} object
#'
#' @param universePolygon A two column matrix specifying the x and y
#' coordinates of the polygon corners, that describe the growth environment
#' boundaries. Alternatively, a character indicating one of the polygon presets
#' can be provided (see details).
#' @param gridFieldSize double. Distance between neighboring environments 3D
#' mesh field elements (rhombic dodecahedrons) in \eqn{\mu}m.
#' @param gridFieldLayers integer. z-dimension (height) as the number of layers
#' of field elements.
#' @param deltaTime double specifying the length of each time step for the
#' simulation in hours.
#' @param rMotion double. Maximum distance a cell can travel by means of
#' Brownian motion in \eqn{\mu} per minute. Default: 0.1 \eqn{\mu}m
#'
#' @return Object of class \link{growthSimulation}.
#'
#' @details
#' Available universe polygon presets:
#' \itemize{
#'   \item{"Petri_<R>"}{ is a Petri dish-like object (actually a 99-corner polygon),
#'   where `<R>` should be replaced with an integer, indicating the radius of the
#'   dish in \eqn{\mu}m.}
#'   \item{"Rectangle_<X>_<Y>"}{ is a, *surprise*, rectangle. `<X>` and `<Y>` should be
#'   integers specifying the width and height in \eqn{\mu}m, respectively.}
#'   \item{"Kiel_<L>"}{ let microbes thrive within Kiel's city limits. Use `<L>`
#'   to specify the latitude dimension in \eqn{\mu}m (integer). The longitude is automatically
#'   scaled accordingly.}
#' }
#'
#' @examples
#' # Construction a square environment of dimensions 100\eqn{\mu}m x 120\eqn{\mu}m x 3\eqn{\mu}m
#' sim <- init_simulation(cbind(c(-50, -50, 50, 50),
#'                              c(-60, 60, 60, -60)),
#'                        gridFieldSize = 1, gridFieldLayers = 3)
#' sim <- init_simulation("rectangle_100_120", gridFieldSize = 1,
#'                        gridFieldLayers = 3)
#'
#' # Construct a Petri dish-like simulation environment (radius: 75 \eqn{\mu}m)
#' sim <- init_simulation("Petri_75", gridFieldSize = 1,
#'                        gridFieldLayers = 3)
#'
#' @import data.table
#' @importFrom methods new
#'
#' @export
init_simulation <- function(universePolygon,
                            gridFieldSize   = 1,
                            gridFieldLayers = 3,
                            deltaTime       = 1/6,
                            rMotion         = 0.1) {

  # generate poygon matrix if preset is chosen
  if(is.character(universePolygon))
    universePolygon <- universe_polygon_preset(universePolygon)

  # last polygon coordinate must be identical to first one
  if(any(universePolygon[1,] != universePolygon[nrow(universePolygon),])) {
    universePolygon <- rbind(universePolygon,
                             universePolygon[1,])
  }


  # init growth environment
  environ <- new("growthEnvironment",
                 polygon.coords  = universePolygon,
                 field.size      = gridFieldSize,
                 field.layers    = gridFieldLayers)

  # construct empty cell info table
  cellDT <- data.table(cell = numeric(),
                       type = character(),
                       x = numeric(), y = numeric(),
                       x.vel = numeric(), y.vel = numeric(),
                       mass = numeric(),
                       size = numeric(),
                       parent = numeric())
  rcdt   <- data.table()

  # init actual growth simulation object
  simobj <- new("growthSimulation",
                n_rounds = 0,
                deltaTime = deltaTime,
                rMotion = rMotion,
                models = list(),
                history = list(),
                cellDT = cellDT,
                rcdt = rcdt,
                universePolygon = universePolygon,
                environ = environ,
                recordDir = NA_character_
  )

  return(simobj)
}

universe_polygon_preset <- function(universePolygon) {

  universePolygon <- universePolygon[1]

  # Petri dish (99 corner polygon)
  if(grepl("^[P|p]etri_[0-9]+$",universePolygon)) {
    p_r <- as.numeric(sub("^[P|p]etri_(\\S+)$", "\\1", universePolygon))

    if(p_r < 2.5)
      stop("Petri dish radius too small. (min 2.5 \u03BCm)")

    petri <- matrix(0, ncol = 2, nrow = 100)

    petri[,1] <- sin(seq(0, 2*pi, length.out = 100))
    petri[,2] <- cos(seq(0, 2*pi, length.out = 100))
    petri[100,] <- petri[1,]

    petri <- petri * p_r

    return(petri)
  }

  # Square dish (99 corner polygon)
  if(grepl("^[R|r]ectangle_[0-9]+_[0-9]+$",universePolygon)) {
    x <- as.numeric(sub("^[R|r]ectangle_(\\S+)_[0-9]+$", "\\1", universePolygon)) / 2
    y <- as.numeric(sub("^[R|r]ectangle_[0-9]+_(\\S+)$", "\\1", universePolygon)) / 2

    if(x < 5 | y < 5)
      stop("Rectangle dimesions too small. (x,y >= 5 \u03BCm)")

    upg <- cbind(c(x,x,-x,-x),
                 c(y,-y,-y,y))

    return(upg)
  }

  # Kiel city limits (proof of principle)
  if(grepl("^[K|k]iel_[0-9]+$",universePolygon)) {
    Kiel <- as.matrix(fread(system.file("extdata","kiel.csv", package = "Eutropia")))
    Kiel_dim <- max(Kiel[,1])
    Kiel <- Kiel / Kiel_dim

    x <- as.numeric(sub("^[K|k]iel_(\\S+)$", "\\1", universePolygon)) / 2
    Kiel <- Kiel * x/2

    if(x < 10)
      stop("Kiel latitude dimension too small (min 10 \u03BCm).")

    return(Kiel)
  }

  stop(paste0("The polygon preset \"",universePolygon,
              "\" does not exist / is not yet supported."))
}


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#  Add organism cells to simulation object  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' Add a model/organism to simulation
#'
#' Adds an organism an its genome-scale metabolic network model to the growth simulation object.
#'
#' @param object S4-object of type \link{growthSimulation}.
#' @param model The organisms metabolic model of S4-type \link[sybil]{modelorg}
#' @param name Character for the name of the model, that will also be used for
#' plotting.
#' @param ncells integer. Number of initial cells to be added to the growth
#' simulation.
#' @param coords (optional) A two column numerical matrix specifying the
#' coordinates (1st column x, 2nd column y) of the initial cells. If provided,
#' the number of rows should be equal to \code{ncells}. Default: NULL
#' @param distribution.method If `coords` is `NULL`, this parameter specifies
#' the distribution method for initial cells. Default: "random_centroid"
#' @param distribution.center Numeric vector of length 2, which specifies the
#' coordinates of the centre for the `distribution.method`.
#' @param distribution.radius double. Spcifies the radius (in \eqn{\mu}m) in which
#' initial cells are distributed.
#' @param cellDiameter double. Diameter in \eqn{\mu}m of initial cells.
#' @param cellMassInit double. Mass in pg of initial cells. Default is 0.28 pg
#' @param cellMassAtDivision double. Cell mass at which a cell divides into two
#' daughter cells. Default: 0.56 pg
#' @param cellShape character. Shape of cells. Currently only "coccus" is
#' supported.
#' @param vmax double. Maximum velocity of a cell in \eqn{\mu}m per second.
#' @param scavengeDist double. Distance in \eqn{\mu}m a cell can scavenge nutrients from
#' its surrounding/microoenvironment.
#' @param rm.deadends If TRUE, dead-end metabolites and reactions are removed
#' from the `model`, which reduces the computation time for FBA, but has
#' otherwise no effect on the flux distribution solutions.
#' @param chemotaxisCompound Character vector of compound IDs, that are signals
#' for directed movement of the organism.
#' @param chemotaxisStrength Numeric vector that indicates the strength of
#' chemotaxis. Positive value for attraction; Negative for repelling effect. A
#' value of 1 indicates that in case of a maximum gradient (concentration-weighted
#' center in cell's scavenge area is at the edge of the area) the cell moves
#' with its maximum speed (vmax) in the direction of the gradient. Default: 0.01
#' @param chemotaxisHillKA Numeric vector for K_A value (unit: mM) in Hill
#' equation in chemotactic metabolite sensing. Default: 0.1 mM
#' @param chemotaxisHillCoef Numeric vector for the Hill coefficient (unitless)
#' in metabolite sensing. Default: 1.2
#' @param open.bounds Numeric value that is used to reset lower bounds of
#' exchange reactions, which have a current lower bound of 0. See Details.
#' @param color Color of organism in visualizations.
#'
#' @details
#' Genome-scale metabolic models usually come pre-constraint, which means that
#' lower bounds for exchange reactions (= max. uptake rates) are set to represent, both,
#' (a) a specific growth environment and (b) the physiological limit of nutrient uptake.
#' Yet, lower bounds that have a value of 0 might also be utilizable by the
#' organism if the compound is present in the environment.
#' If the option `open.bounds` is used, those 0-lower bounds are replaced with a
#' new lower bound to enable the potential uptake in the agent-based simulation.
#' Please note that the value should by convention be negative; however
#' this package changes the value to it's negative counterpart if a positive value
#' is provided.
#'
#' The default cell diameter (\eqn{(3 * 1 / (4 * pi))^(1/3) * 2}) is that of a
#' sphere with 1 \eqn{\mu}m^3 volume.
#'
#' 'chemotaxisHillKA' and 'chemotaxisHillCoef' are metabolite sensing sensitivity
#' parameters, which is modeled as a Hill equation. Default values correspond to
#' numbers estimated by Sourjik and Berg (2001, PNAS) for \emph{Escherichia coli}.
#'
#' @return Object of class \link{growthSimulation}.
#'
#' @references
#'  \url{https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100008} \cr
#'  \url{http://book.bionumbers.org/how-big-is-an-e-coli-cell-and-what-is-its-mass/} \cr
#'  \url{https://bionumbers.hms.harvard.edu/bionumber.aspx?id=115616&ver=0&trm=speed+e.+coli&org=} \cr
#'  Victor Sourjik and Howard C. Berg. (2001). Receptor sensitivity in bacterial
#'  chemotaxis. \emph{PNAS} \strong{99}, 123-127.
#'
#' @examples
#' # add two bacterial models (Eubacterium rectale, Bifidobacterium longum)
#' # to the environment; each with 15 initial cells
#' models <- list()
#' models[['eure']] <- readRDS(system.file("extdata", "eure.RDS",
#'                             package="Eutropia"))
#' models[['bilo']] <- readRDS(system.file("extdata", "bilo.RDS",
#'                             package="Eutropia"))
#'
#' sim <- init_simulation(cbind(c(-100, -100, 100, 100),
#'                              c(-100, 100, 100, -100)),
#'                        gridFieldSize = 1.75, gridFieldLayers = 3)
#'
#' sim <- add_organism(sim, model = models[["eure"]], name = "E. rectale",
#'                     ncells = 15, distribution.radius = 30)
#' sim <- add_organism(sim, model = models[["bilo"]], name = "B. longum",
#'                     ncells = 15, distribution.radius = 30)
#'
#' plot_cells(sim, xlim = c(-50,50), ylim= c(-50,50))
#'
#' @import particles
#' @import tidygraph
#' @importFrom methods new
#' @import sf
#'
#' @export
add_organism <- function(object,
                         model,
                         name,
                         ncells,
                         coords = NULL,
                         distribution.method = "random_centroid",
                         distribution.center = NULL,
                         distribution.radius = NULL,
                         cellDiameter = (3 * 1 / (4 * pi))^(1/3) * 2, # diameter of sphere with 1 micro-m^3 volume
                         cellMassInit = 0.28,
                         cellMassAtDivision = 0.56,
                         cellShape = "coccus",
                         vmax = 11,
                         scavengeDist = cellDiameter*2.5,
                         rm.deadends = T,
                         chemotaxisCompound = NULL,
                         chemotaxisStrength = 0.01,
                         chemotaxisHillKA   = 0.1,
                         chemotaxisHillCoef = 1.2,
                         open.bounds = NULL,
                         color = NULL) {
  # sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")
  if(!(inherits(model, "modelorg")))
    stop("'model' not of class 'modelorg'.")
  if(ncells <= 0 | (ncells %% 1) != 0)
    stop("'ncells' should be a positive non-zero integer.")


  if(is.null(chemotaxisCompound)) {
    chemotaxisCompound <- character(0)
    chemotaxisStrength <- double(0)
    chemotaxisHillKA   <- double(0)
    chemotaxisHillCoef <- double(0)
  } else {
    if(!(length(chemotaxisStrength) %in% c(1,length(chemotaxisCompound))))
      stop("'chemotaxisStrength' should be the same length as chemotaxisCompound or 1.")
    if(!(length(chemotaxisHillKA) %in% c(1,length(chemotaxisCompound))))
      stop("'chemotaxisHillKA' should be the same length as chemotaxisCompound or 1.")
    if(!(length(chemotaxisHillCoef) %in% c(1,length(chemotaxisCompound))))
      stop("'chemotaxisHillCoef' should be the same length as chemotaxisCompound or 1.")

    if(any(chemotaxisHillKA <= 0) | any(chemotaxisHillCoef <= 0))
      stop("'chemotaxisHillKA' and 'chemotaxisHillCoef' must be non-zero positive numbers.")

    if(length(chemotaxisStrength) == 1)
      chemotaxisStrength <- rep(chemotaxisStrength, length(chemotaxisCompound))
    if(length(chemotaxisHillKA) == 1)
      chemotaxisHillKA <- rep(chemotaxisHillKA, length(chemotaxisCompound))
    if(length(chemotaxisHillCoef) == 1)
      chemotaxisHillCoef <- rep(chemotaxisHillCoef, length(chemotaxisCompound))
  }

  # set color automatically if NULL
  if(is.null(color)) {
    default.colors <- c("#D55E00","#56B4E9","#F0E442","#009E73","#0072B2",
                        "#E69F00","#CC79A7","#DEDEDE")
    colors.in.use <- unlist(lapply(object@models, function(x) x@color))
    color <- default.colors[!(default.colors %in% colors.in.use)]
    if(length(color) == 0) {
      color <- "#333333"
      warnings("Eutropia ran out of distinct and colorblind-friendly colors. Assigning gray instead.")
    } else {
      color <- color[1]
    }
  }

  # init new organism object
  object@models[[name]] <- new("Organism",
                               cellDiameter = cellDiameter,
                               cellMassInit = cellMassInit,
                               cellMassAtDivision = cellMassAtDivision,
                               vmax = vmax,
                               scavengeDist = scavengeDist,
                               mod = model,
                               rm.deadends = rm.deadends,
                               chemotaxisCompound = chemotaxisCompound,
                               chemotaxisStrength = chemotaxisStrength,
                               chemotaxisHillKA = chemotaxisHillKA,
                               chemotaxisHillCoef = chemotaxisHillCoef,

                               open.bounds = open.bounds,
                               color = color)

  #if not provided -> assign cell positions:
  if(is.null(coords)) {
    if(!(distribution.method %in% c("random","random_centroid")))
      stop("Distribution method not supported. Choose one of \"random\",\"random_centroid\".")

    universePG <- st_polygon(list(object@universePolygon))

    if(distribution.method == "random") {
      #coords <- spsample(universePG, ncells, type = "random")@coords
      coords <- st_sample(universePG, ncells, type = "random")
    }

    if(distribution.method == "random_centroid") {

      PGcent <- distribution.center
      if(is.null(PGcent)) {
        PGcent <- st_centroid(universePG)
      } else {
        PGcent <- st_point(PGcent)
      }

      maxRadius <- distribution.radius

      if(is.null(maxRadius)) {
        maxRadius <- st_distance(PGcent, st_cast(universePG,
                                                     "LINESTRING"))[1,1]
      }
      distr_area <- st_buffer(PGcent, dist = maxRadius)
      distr_area <- st_intersection(universePG, distr_area)

      coords <- st_sample(distr_area, ncells, type = "random")
    }

    coords <- lapply(coords, as.matrix)
    coords <- do.call(rbind, coords)
  }

  newcells <- data.table(cell = (1:ncells)+nrow(object@cellDT),
                         type = name,
                         x = coords[,1], y = coords[,2],
                         x.vel = 0, y.vel = 0,
                         mass = cellMassInit,
                         size = cellDiameter,
                         parent = NA_integer_)

  # making sure cells are placed within the universe polygon (trapping cells)
  # oder: "Schäfchen zurück ins Gehege holen"
  sim_addorg <- tbl_graph(nodes = newcells, directed = F, node_key = "cell") %>%
    simulate(setup = predefined_genesis(x = newcells$x,
                                        y = newcells$y,
                                        x_vel = newcells$x.vel,
                                        y_vel = newcells$y.vel)) %>%
    # wield(collision_force, radius = newcells$size/2, n_iter = 50, strength = 0.7) %>%
    wield(trap_force, polygon = object@universePolygon, strength = 0.7,
          min_dist = 5, distance_falloff = 2) %>%
    impose(velocity_constraint,
           vmax = rep(1, nrow(newcells))) %>%
    impose(polygon_constraint,
           polygon = object@universePolygon)

  keepJittering <- T
  tmp_currentpos <- c(newcells$x, newcells$y)
  while(keepJittering) {
    sim_addorg <- evolve(sim_addorg, steps = 10)
    if(all(tmp_currentpos == c(sim_addorg$position[,1],sim_addorg$position[,2]))) {
      keepJittering <- F
    } else {
      # handbreak... (setting velocity to 0)
      sim_addorg$velocity[,1] <- sim_addorg$velocity[,1]*0
      sim_addorg$velocity[,2] <- sim_addorg$velocity[,2]*0
      tmp_currentpos <- c(sim_addorg$position[,1],sim_addorg$position[,2])
    }
  }

  newcells$x <- sim_addorg$position[,1]
  newcells$y <- sim_addorg$position[,2]

  # add cells to cell info table
  object@cellDT <- rbind(object@cellDT, newcells)

  # add compounds, for which there are exchange reactions in at least one model
  # and which are not yet part of the environment
  mod <- object@models[[name]]@mod

  dt_exr <- data.table(id = mod@react_id,
                       lb = mod@lowbnd,
                       name = mod@react_name)
  dt_exr <- dt_exr[grepl("^EX_", id)]
  dt_exr[, id := gsub("^EX_","",id)]
  dt_exr <- dt_exr[id != "cpd11416_c0"]
  dt_exr[, name := gsub("-e0-e0 Exchange$","",name)]
  dt_exr[, name := gsub("-e0 Exchange$","",name)]
  dt_exr[, name := gsub(" Exchange$","",name)]
  dt_exr[, name := gsub(" exchange$","",name)]

  dt_exr <- dt_exr[!(id %in% object@environ@compounds)] # only new compounds should be added

  if(nrow(dt_exr) > 0) {
    object <- add_compounds(object,
                            compounds = dt_exr$id,
                            concentrations = rep(0, nrow(dt_exr)),
                            compound.names = dt_exr$name,
                            is.constant = rep(FALSE, nrow(dt_exr)))
  }

  return(object)
}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# show method for small summary #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

#' @title Print a short summary of a growth simulation
#'
#' @description Displays a few numbers to describe the current status of a growth simulation.
#'
#' @param object S4-object of type \link{growthSimulation}.
#'
#' @import data.table
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
            print(object@cellDT[, .(`cells` = .N, `mass in pg` = sum(`mass`)), by = "type"])
            cat("\n")

            # Environment
            cat("Cell growth environment\n")
            cat("    Universe volume (\u03BCm^3):\t\t", round(object@environ@fieldVol * object@environ@nfields, digits = 2) , "\n")
            cat("    Number of rhombic dodecahedrons:\t",object@environ@nfields,"\n")
            cat("    Number of compounds:\t\t",length(object@environ@compounds),
                "(", sum(!object@environ@conc.isConstant),"variable )\n")

          }
)



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Add compounds to environment  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Add compounds to the growth environment
#'
#' @description The function can be used to add substances to the growth
#' environment by providing concentrations in mM.
#'
#' @param object S4-object of type \link{growthSimulation}.
#' @param compounds Character vector with the compound IDs of substances to add
#' to the environment. Compound IDs should correspond to the models' exchange
#' reactions IDs ("EX_[cpdid]"), without the "EX_" prefix.
#' @param concentrations Numeric vector with the concentrations of the compounds
#' from `compounds` in the same order. Values in mM.
#' @param compound.names Character vector with the compound names.
#' @param is.constant Logical vector that indicates if the compound should
#' remain constant over time despite of potential uptake or production by cells.
#' @param compound.D Numeric vector with the compounds' diffusion coefficients
#' in \eqn{\mu}m^2/s. Default: 75
#'
#' @details Compound concentration are equally distributed across the whole
#' growth environment.
#' If the compound is already present, old and new concentrations are added.
#' More options are planned.
#'
#' If no compound names are provided, the current names are kept (if compound
#' is already present) or the compound ID is also used as name (in case the
#' compound is new).
#'
#' You can also define "Inf" for the compound diffusion rates in 'compound.D'.
#' This has the effect, that the compound will be evenly distributed across the
#' whole growth environment again at every iteration during the simulation.
#'
#' @return Return a S4-object of type \link{growthSimulation}.
#'
#' @examples
#' sim <- init_simulation(cbind(c(-100, -100, 100, 100), c(-100, 100, 100, -100)),
#'                        gridFieldSize = 2, gridFieldLayers = 3)
#'
#' sim <- add_compounds(sim,
#'                      compounds = c("cpd00027_e0","cpd00029_e0","cpd00047_e0",
#'                                    "cpd00159_e0","cpd00211_e0"),
#'                      concentrations = c(50,0,0,0,0),
#'                      compound.names = c("D-Glucose","Acetate","Formate",
#'                                         "L-Lactate","Butyrate"),
#'                      is.constant = rep(FALSE, 5))
#'
#' @export
add_compounds <- function(object, compounds, concentrations,
                          compound.names = NULL,
                          is.constant = NULL,
                          compound.D  = NULL) {
  # Sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  default_D <- 75

  if(length(compounds) != length(concentrations))
    stop("Lengths of 'compounds' and 'concentrations' differ.")

  # make compounds non-constant if nothing else is specified
  if(is.null(is.constant)) {
    is.constant <- rep(F, length(compounds))
  }

  # assign each compound the default diffusion coefficient if nothing
  # else is specified
  if(is.null(compound.D)) {
    compound.D <- rep(default_D, length(compounds))
  } else {
    compound.D <- ifelse(is.na(compound.D), default_D, compound.D)
  }

  # if only one diff. coeff. is specified -> use it for all metabolites
  if(length(compound.D) == 1) {
    compound.D <- rep(compound.D, length(compounds))
  }

  # is.constant should be NULL (see above) or of length 1 or n (nr. of compounds)
  if(length(is.constant) != 1 & length(is.constant) != length(compounds)) {
    stop("Length of 'is.constant' is not 1 or the same length as 'compounds'")
  }

  # compound.D should be NULL (see above) or of length 1 or n (nr. of compounds)
  if(length(compound.D) != 1 & length(compound.D) != length(compounds)) {
    stop("Length of 'compound.D' is not 1 or the same length as 'compounds'")
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
      object@environ@compound.D[c_ind]       <- compound.D[i]
      if(!is.na(compound.names[i]))
        object@environ@compound.names[c_ind] <- compound.names[i]

      # Metabolite is new
    } else {
      object@environ@compounds       <- c(object@environ@compounds, compounds[i])
      object@environ@concentrations  <- cbind(object@environ@concentrations,
                                              matrix(concentrations[i], ncol = 1, nrow = object@environ@nfields))
      object@environ@conc.isConstant <- c(object@environ@conc.isConstant, is.constant[i])
      object@environ@compound.D      <- c(object@environ@compound.D, compound.D[i])
      if(!is.na(compound.names[i])) {
        if(compound.names[i] %in% object@environ@compound.names) {
          warning(paste0("Compound with name '",compound.names[i],"' already exists. Replacing with '",compounds[i],"'."))
          compound.names[i] <- compounds[i]
        }
        object@environ@compound.names <- c(object@environ@compound.names, compound.names[i])
      } else {
        object@environ@compound.names <- c(object@environ@compound.names, compounds[i])
      }
    }
  }

  return(object)
}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Add compound gradient         #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Add compounds to the growth environment in a gradient
#'
#' @description The function can be used to add substances to the growth
#' environment where the compound is distribution in a concentration gradient.
#'
#' @param object S4-object of type \link{growthSimulation}.
#' @param compound Character with the compound ID of substance to add
#' to the environment in a gradient. Compound ID should correspond to the
#' models' exchange reaction ID ("EX_[cpdid]"), without the "EX_" prefix.
#' @param p1 Numeric vector of length 2 or 3 defining the xy(z)-coordinates of
#' the first gradient reference point. See details.
#' @param p2 Numeric vector of length 2 or 3 defining the xy(z)-coordinates of
#' the second gradient reference point. See details.
#' @param c1 Numeric value with the concentration of the compound at 'p1'.
#' Values in mM.
#' @param c2 Numeric value with the concentration of the compound at 'p2'.
#' Values in mM.
#' @param gradient.dir Either "radial", "linear", or "linear_mirrored". See
#' details.
#' @param compound.name Character with the compound name.
#' @param is.constant Logical defining if the compound should remain constant
#' over time despite of potential uptake or production by cells.
#' @param compound.D Numeric value with the compound's diffusion coefficient
#' in \eqn{\mu}m^2/s. Default: 75
#'
#' @details If only xy-coordinated are provided the z-coordinate is assumed to
#' be 0.
#'
#' If ''gradient.dir' is set to "linear", any point that is in the
#' opposite direction of the gradient will get the concentration c1. If
#' "linear_mirrored", the concentration gradient is mirrored at the plane that
#' is defined by p1 and the p1-p2 as normal vector.
#'
#' If the compound is already present, old and new concentrations are added.
#'
#' If no compound names are provided, the current names are kept (if compound
#' is already present) or the the compound ID is also used as name (in case the
#' compound is new).
#'
#'
#' @return Return a S4-object of type \link{growthSimulation}.
#'
#' @examples
#' sim <- init_simulation(cbind(c(-70, -70, 70, 70), c(-45, 45, 45, -45)),
#' gridFieldSize = 2, gridFieldLayers = 3)
#' sim <- add_compound_gradient(sim,
#'                              compound = "cpd00027_e0",
#'                              p1 = c(-60,33), p2 = c(60,-40),
#'                              c1 = 25, c2 = 0,
#'                              compound.name = "D-Glucose")
#' sim <- add_compound_gradient(sim,
#'                              compound = "cpd00036_e0",
#'                              p1 = c(60,-20), p2 = c(50,60),
#'                              c1 = 27, c2 = 0,
#'                              gradient.dir = "linear",
#'                              compound.name = "Succinate")
#'
#' plot_environment(sim, c("cpd00027_e0","cpd00036_e0"), incl.timestamp = FALSE)
#' @import sf
#' @export
add_compound_gradient <- function(object, compound, p1, p2, c1, c2,
                                  gradient.dir = "radial",
                                  compound.name = NULL,
                                  is.constant = FALSE,
                                  compound.D  = 75) {
  # Sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  # grad dir
  if(!(gradient.dir %in% c("radial","linear","linear_mirrored")))
    stop("gradient.dir should be \"radial\", \"linear\", or \"linear_mirrored\".")

  # Both coordinates should have numeric length of 2 or 3.
  if(!(length(p1) %in% c(2,3) & length(p2) %in% c(2,3))) {
    stop("Coordinates p1 and p2 should each be vectors of length 2 (xy) or 3 (xyz).")
  }
  if(length(p1) == 2)
    p1 <- c(p1,0)
  if(length(p2) == 2)
    p2 <- c(p2,0)

  # p1 and p2 same?
  if(all(p1 == p2))
    stop("'p1' and 'p2' cannot be identical for defining a gradient.")

  # Calculate concentration gradient
  p1sf <- st_point(p1)
  p2sf <- st_point(p2)
  d <- st_distance(p1sf, p2sf)[1,1]
  if(gradient.dir == "radial") {
    envp <- st_cast(st_sfc(object@environ@field.pts), "POINT")
    env_dist <- st_distance(envp, p1sf)
  }
  if(grepl("^linear",gradient.dir)) {
    envp <- as.matrix(object@environ@field.pts)

    norm_vec <- p2-p1

    env_dist <- norm_vec[1]*(envp[,1] - p1[1]) + norm_vec[2]*(envp[,2] - p1[2]) + norm_vec[3]*(envp[,3] - p1[3])
    #env_dist <- norm_vec[1]*envp[,1] + norm_vec[2]*envp[,2] + norm_vec[3]*envp[,3]
    #env_dist <- env_dist + norm_vec[1]*p1[1] + norm_vec[2]*p1[2] + norm_vec[3]*p1[3]
    env_dist <- env_dist / sqrt(sum(norm_vec^2))

    if(grepl("mirrored$", gradient.dir)) {
      env_dist <- abs(env_dist)
    } else {
      env_dist <- ifelse(env_dist < 0, 0, env_dist)
    }
  }

  # translate distances to concentrations
  grad_func <- function(x) {
    m <- (c2 - c1) / d
    y <- x * m + c1
    y <- ifelse(y < min(c(c1,c2)), min(c(c1,c2)), y)
    y <- ifelse(y > max(c(c1,c2)), max(c(c1,c2)), y)
  }
  concadd <- grad_func(env_dist)

  # Add new concentration to growth environment
  if(compound %in% object@environ@compounds) {
    # Metabolite is already present in environment
    c_ind <- which(object@environ@compounds == compound)
    object@environ@concentrations[, c_ind] <- object@environ@concentrations[, c_ind] + concadd
    object@environ@conc.isConstant[c_ind]  <- is.constant
    object@environ@compound.D[c_ind]       <- compound.D
    if(!is.null(compound.name) && !is.na(compound.name))
      object@environ@compound.names[c_ind] <- compound.name

  } else {
    # Metabolite is new
    object@environ@compounds       <- c(object@environ@compounds, compound)
    object@environ@concentrations  <- cbind(object@environ@concentrations,
                                            matrix(concadd, ncol = 1, nrow = object@environ@nfields))
    object@environ@conc.isConstant <- c(object@environ@conc.isConstant, is.constant)
    object@environ@compound.D      <- c(object@environ@compound.D, compound.D)
    if(!is.na(compound.name)) {
      if(compound.name %in% object@environ@compound.names) {
        warning(paste0("Compound with name '",compound.name,"' already exists. Replacing with '",compound,"'."))
        compound.name <- compound
      }
      object@environ@compound.names <- c(object@environ@compound.names, compound.name)
    } else {
      object@environ@compound.names <- c(object@environ@compound.names, compound)
    }
  }

  return(object)
}


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Dilute compounds            #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Dilute compounds
#'
#' @description Dilutes all or selected compounds with a given dilution factor.
#'
#' @param object A \link{growthSimulation} object
#' @param dilution.factor Numeric within the range of [0,1], by which the
#' compound concentrations are diluted. `1` completely dilutes concentrations to
#' 0 mM, while `0` does not change anything.
#' @param compounds Character of compound IDs that should be diluted. If `NULL`,
#' all compounds are diluted.
#' @param incl.constant Logical specifying whether also constant compounds
#' should be diluted. Default: FALSE
#'
#' @export
dilute_compounds <- function(object, dilution.factor,
                             compounds = NULL, incl.constant = FALSE) {

  # sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")
  if(dilution.factor > 1 | dilution.factor < 0) {
    stop("Dilution factor should be between 0 and 1.")
  }

  if(is.null(compounds))
    compounds <- object@environ@compounds

  available_compounds <- object@environ@compounds
  if(!incl.constant)
    available_compounds <- object@environ@compounds[!object@environ@conc.isConstant]

  compounds <- compounds[compounds %in% available_compounds]

  if(length(compounds) == 0) {
    warning("No valid compounds selected. Returning original simulation object.")
    return(object)
  }

  ind_relcpds <- match(compounds, object@environ@compounds)

  object@environ@concentrations[, ind_relcpds] <- object@environ@concentrations[, ind_relcpds] * (1-dilution.factor)

  return(object)
}


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Summary exchanges           #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Summary of uptake / production by organism type
#'
#' @description Uptake/Production rates in fmol summarized by organism type
#'
#' @param object S4-object of type \link{growthSimulation}.
#' @param iter Positive integer number of the simulation step/iteration to plot
#' the rates.
#'
#' @return A data.table
#'
#' @export
summary_exchanges <- function(object, iter = NULL) {

  # sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  # Sanity checks
  if(object@n_rounds < 1)
    stop("Simulation did not yet run for at least 1 iteration. No exchange rates availabe (yet).")

  if(!is.null(iter)) {
    if(iter > object@n_rounds) {
      warning(paste0("Simulation did not run ",iter," iterations yet. Displaying results for last iteration (",object@n_rounds,")"))
      iter <- object@n_rounds
    }
  } else {
    iter <- object@n_rounds
  }

  # data.table X R CMD check workaround
  type <- compound <- cell <- exflux <- compound.name <- fmol <- NULL

  dt_exch <- merge(object@history[[iter]]$cell.exchanges,
                   object@history[[iter]]$cells[,.(cell, type)],
                   by = "cell")
  dt_exch <- dt_exch[, .(fmol = sum(exflux)), by = .(type, compound)]
  dt_exch[, compound := gsub("^EX_","",compound)]
  dt_exch <- dt_exch[compound != "cpd11416_c0"] # Biomass...
  dt_exch$compound.name <- object@environ@compound.names[match(dt_exch$compound,
                                                               object@environ@compounds)]

  return(dt_exch[,.(type, compound, compound.name, fmol)])

}
