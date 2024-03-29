
#' Structure of the S4 class "Organism"
#'
#' @aliases Organism
#'
#' @exportClass Organism
#'
#' @slot cellDiameter Numeric for the initial Diameter of cell spheres in \eqn{\mu}m
#' @slot cellMassInit Numeric for the initial cell mass in pg
#' @slot cellMassAtDivision Limit of a cell's mass before it divided into two
#' daughter cells. Unit: pg
#' @slot cellShape Character for the cell shape type. Currently, only coccus/sphere
#' shapes are supported.
#' @slot vmax Numeric for the maximum speed a cell can move in \eqn{\mu}m/s
#' @slot scavengeDist Numeric indicating the maximum distance (from cell surface)
#' a cell can scavenge nutrients from its surrounding. Unit: \eqn{\mu}m
#' @slot chemotaxisCompound Character vector with the compound IDs that influence
#' the cells chemotaxis behavior.
#' @slot chemotaxisStrength Numeric vector that indicates the strength of
#' chemotaxis. Positive value for attraction; Negative for repelling effect. A
#' value of 1 indicates that in case of a maximum gradient (concentration-weighted
#' center in cell's scavenge area is at the edge of the area) the cell moves
#' with its maximum speed (vmax) in the direction of the gradient.
#' @slot chemotaxisHillKA Numeric vector for K_A value in Hill equation in
#' chemotactic metabolite sensing. Unit: mM
#' @slot chemotaxisHillCoef Numeric vector for the Hill coefficient. Unitless
#' @slot mod Object of S4-class \link[sybil]{modelorg} for the organisms metabolic
#' network model.
#' @slot exoenzymes Character vector of IDs of the organism's exoenzymes
#' @slot exoenzymes.prod Numeric vector of the production rates of exoenzymes.
#' Unit: nmol / gDW / hr (nmol Enzyme per gDW cells per hr)
#' @slot color Color of organism in visualizations.
setClass("Organism",

         slots = c(
           # biological parameters
           cellDiameter       = "numeric", # micro-m
           cellMassInit       = "numeric", # pg
           cellMassAtDivision = "numeric", # pg
           cellShape          = "character", # currently only coccus
           vmax               = "numeric", # in micro-m/s
           scavengeDist       = "numeric", # in micro-m

           # Chemotaxis
           chemotaxisCompound = "character",
           chemotaxisStrength = "numeric", # positive (attracting) and negative (repelling)
           chemotaxisHillKA   = "numeric",
           chemotaxisHillCoef = "numeric",

           # Genome-scale model for FBA / pFBA
           mod      = "modelorg",

           # Exoenzymes
           exoenzymes = "character", # vector with Exoenzyme IDs
           exoenzymes.prod = "numeric", # Vector for production rates in [nmol/gDW/hr]

           # Aestetics
           color = "character"
         )
)

#' @import data.table
#' @import sybil
setMethod("initialize", "Organism",
          function(.Object,
                   cellDiameter,
                   cellMassInit,
                   cellMassAtDivision,
                   vmax,
                   scavengeDist,
                   mod,
                   rm.deadends,
                   chemotaxisCompound,
                   chemotaxisStrength,
                   chemotaxisHillKA,
                   chemotaxisHillCoef,
                   open.bounds,
                   color,
                   ...) {
            .Object <- callNextMethod(.Object, ...)

            if(cellMassAtDivision <= cellMassInit)
              stop("Initial cell mass should be less than mass at binary fission.")

            # Organism-specific parameters
            .Object@cellDiameter <- cellDiameter
            .Object@cellMassInit <- cellMassInit
            .Object@cellMassAtDivision <- cellMassAtDivision
            .Object@cellShape <- "coccus" # (currently only coccus possible)
            .Object@vmax <- vmax
            .Object@scavengeDist <- scavengeDist
            .Object@chemotaxisCompound <- chemotaxisCompound
            .Object@chemotaxisStrength <- chemotaxisStrength
            .Object@chemotaxisHillKA   <- chemotaxisHillKA
            .Object@chemotaxisHillCoef <- chemotaxisHillCoef
            .Object@color <- color

            .Object@exoenzymes <- character(0)
            .Object@exoenzymes.prod <- double(0)

            # Rm exchange reaction for D-lactate if L-lactate is also present
            # (works only with gapseq models)
            if(all(c("EX_cpd00221_e0","EX_cpd00159_e0") %in% mod@react_id))
              mod <- rmReact(mod, react = "EX_cpd00221_e0")

            # open bounds if wanted
            if(!is.null(open.bounds)) {
              if(open.bounds > 0)
                open.bounds <- open.bounds * -1
              dt_ex <- data.table(rxn = mod@react_id, lb = mod@lowbnd)
              dt_ex <- dt_ex[lb == 0 & grepl("^EX_", rxn)]
              mod <- changeBounds(mod, react = dt_ex$rxn,
                                  lb = rep(open.bounds, nrow(dt_ex)))
            }

            # Network
            if(rm.deadends) {
              der <- deadEndMetabolites(mod)$der
              mod <- rmReact(mod, react = der)
            }
            .Object@mod <- mod

            return(.Object)

          }
)
