setClass("Organism",

         slots = c(
           # biological parameters
           cellDiameter       = "numeric", # µm
           cellMassInit       = "numeric", # pg
           cellMassAtDivision = "numeric", # pg
           cellShape          = "character", # currently only coccus
           vmax               = "numeric", # ?
           scavengeDist       = "numeric", # in µm

           # Genome-scale model for FBA / pFBA
           mod      = "modelorg"
         )
)

setMethod("initialize", "Organism",
          function(.Object,
                   cellDiameter,
                   cellMassInit,
                   cellMassAtDivision,
                   vmax,
                   scavengeDist,
                   mod,
                   rm.deadends,
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

            # Rm exchange reaction for D-lactate if L-lactate is also present
            if(all(c("EX_cpd00221_e0","EX_cpd00159_e0") %in% mod@react_id))
              mod <- rmReact(mod, react = "EX_cpd00221_e0")

            # Network
            if(rm.deadends) {
              der <- deadEndMetabolites(mod)$der
              mod <- rmReact(mod, react = der)
            }
            .Object@mod <- mod

            return(.Object)

          }
)
