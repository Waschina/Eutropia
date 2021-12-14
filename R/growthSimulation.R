
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
#' random movement per simulation round. Unit: µm per minute.
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
#' coordinates of the polygon, that describes the growth environment boundaries.
#' Alternatively, a character indicating one of the polygon presets can be
#' provides (see details).
#' @param gridFieldSize double. Distance between neighboring environments 3D
#' mesh field elements (rhombic dodecahedrons) in µm.
#' @param gridFieldLayers integer. z-dimension (height) as the number of layers
#' of field elements.
#' @param deltaTime double specifying the length of each time step for the
#' simulation in hours.
#' @param rMotion double. Maximum distance a cell can travel by means of
#' Brownian motion per minute Default: 0.1 µm
#'
#' @return Object of class \link{growthSimulation}.
#'
#' @details
#' Available universe polygon presets:
#' \itemize{
#'   \item{"Petri_<R>"}{ is a Petri dish-like object (actually a 99-corner polygon),
#'   where `<R>` should be replaced with an integer, indicating the radius of the
#'   dish in µm.}
#'   \item{"Rectangle_<X>_<Y>"}{ is a, *surprise*, rectangle. `<X>` and `<Y>` should be
#'   integers specifying the width and height in µm, respectively.}
#'   \item{"Kiel_<L>"}{ let microbes thrive within Kiel's city limits. Use `<L>`
#'   to specify the latitude dimension in µm (integer). The longitude is automatically
#'   scaled accordingly.}
#' }
#'
#' @examples
#' # Construction a square environment of dimensions 100µm x 100µm x 5µm
#' sim <- init_simulation(cbind(c(-50, -50, 50, 50),
#'                              c(-60, 60, 60, -60)),
#'                        gridFieldSize = 1, gridFieldLayers = 3)
#' sim <- init_simulation("rectangle_100_120", gridFieldSize = 1,
#'                        gridFieldLayers = 5)
#'
#' # Construct a Petri dish-like simulation environment (radius: 100 µm)
#' sim <- init_simulation("Petri_100", gridFieldSize = 1,
#'                        gridFieldLayers = 10)
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
      stop("Petri dish radius too small. (min 2.5 µm)")

    petri <- matrix(0, ncol = 2, nrow = 100)

    petri[,1] <- sin(seq(0, 2*pi, length.out = 100))
    petri[,2] <- cos(seq(0, 2*pi, length.out = 100))

    petri <- petri * p_r

    return(petri)
  }

  # Square dish (99 corner polygon)
  if(grepl("^[R|r]ectangle_[0-9]+_[0-9]+$",universePolygon)) {
    x <- as.numeric(sub("^[R|r]ectangle_(\\S+)_[0-9]+$", "\\1", universePolygon)) / 2
    y <- as.numeric(sub("^[R|r]ectangle_[0-9]+_(\\S+)$", "\\1", universePolygon)) / 2

    if(x < 5 | y < 5)
      stop("Rectangle dimesions too small. (x,y >= 5 µm)")

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
      stop("Kiel latitude dimension too small (min 10 µm).")

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
#' @param distribution.radius double. Spcifies the radius (in µm) in which
#' initial cells are distributed.
#' @param cellDiameter double. Diameter in µm of initial cells.
#' @param cellMassInit double. Mass in pg of initial cells. Default is 0.28 pg
#' @param cellMassAtDivision double. Cell mass at which a cell divides into two
#' daughter cells. Default: 0.56 pg
#' @param cellShape character. Shape of cells. Currently only "coccus" is
#' supported.
#' @param vmax double. Maximum velocity of a cell in µm per second.
#' @param scavengeDist double. Distance in µm a cell can scavenge nutrients from
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
#' sphere with 1 µm^3 volume.
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
#' sim <- init_simulation(cbind(c(-150, -150, 150, 150),
#'                              c(-150, 150, 150, -150)),
#'                        gridFieldSize = 1.75, gridFieldLayers = 5)
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
#' @import rgeos
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
                         cellDiameter = (3 * 1 / (4 * pi))^(1/3) * 2, # diameter of sphere with 1 µm^3 volume
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
                         open.bounds = NULL) {
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
      stop("'chemotaxisStrength' should be the samle length as chemotaxisCompound or 1.")
    if(!(length(chemotaxisHillKA) %in% c(1,length(chemotaxisCompound))))
      stop("'chemotaxisHillKA' should be the samle length as chemotaxisCompound or 1.")
    if(!(length(chemotaxisHillCoef) %in% c(1,length(chemotaxisCompound))))
      stop("'chemotaxisHillCoef' should be the samle length as chemotaxisCompound or 1.")

    if(any(chemotaxisHillKA <= 0) | any(chemotaxisHillCoef <= 0))
      stop("'chemotaxisHillKA' and 'chemotaxisHillCoef' must be non-zero positive numbers.")

    if(length(chemotaxisStrength) == 1)
      chemotaxisStrength <- rep(chemotaxisStrength, length(chemotaxisCompound))
    if(length(chemotaxisHillKA) == 1)
      chemotaxisHillKA <- rep(chemotaxisHillKA, length(chemotaxisCompound))
    if(length(chemotaxisHillCoef) == 1)
      chemotaxisHillCoef <- rep(chemotaxisHillCoef, length(chemotaxisCompound))
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

                               open.bounds = open.bounds)

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
      if(is.null(PGcent)) {
        PGcent <- SpatialPoints(matrix(universePG@labpt,ncol = 2))
      } else {
        PGcent <- SpatialPoints(matrix(PGcent,ncol = 2))
      }

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

  object <- add_compounds(object,
                          compounds = dt_exr$id,
                          concentrations = rep(0, nrow(dt_exr)),
                          compound.names = dt_exr$name,
                          is.constant = rep(FALSE, nrow(dt_exr)))

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
            # cat("    Universe dimensions (µm):\t\t",
            #     abs(round(min(object@environ@field.pts@coords[,1])-max(object@environ@field.pts@coords[,1]), digits = 2))," x ",
            #     abs(round(min(object@environ@field.pts@coords[,2])-max(object@environ@field.pts@coords[,2]), digits = 2))," x ",
            #     abs(round(min(object@environ@field.pts@coords[,3])-max(object@environ@field.pts@coords[,3]), digits = 2)),"\n")
            cat("    Universe volume (µm^3):\t\t", round(object@environ@fieldVol * object@environ@nfields, digits = 2) , "\n")
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
#' in µm^2/s. Default: 75
#'
#' @details Compound concentration are equally distributed across the whole
#' growth environment.
#' If the compound is already present, old and new concentrations are added.
#' More options are planned.
#'
#' If no compound names are provided, the current names are kept (if compound
#' is already present) or the the compound ID is also used as name (in case the
#' compound is new).
#'
#' @return Return a S4-object of type \link{growthSimulation}.
#'
#' @examples
#' sim <- init_simulation(cbind(c(-100, -100, 100, 100), c(-100, 100, 100, -100)),
#'                        gridFieldSize = 2, gridFieldLayers = 5)
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
        }
        object@environ@compound.names <- c(object@environ@compound.names, compound.names[i])
      } else {
        object@environ@compound.names <- c(object@environ@compound.names, compounds[i])
      }
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
