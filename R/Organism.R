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
           mod      = "modelorg",
           mod.warm = "sysBiolAlg"
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
                   pFBAcoeff,
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

            # Network
            if(rm.deadends) {
              der <- deadEndMetabolites(mod)$der
              mod <- rmReact(mod, react = der)
            }
            .Object@mod <- mod
            .Object@mod.warm <- sysBiolAlg(mod,
                                           algorithm = "mtf2",
                                           pFBAcoeff = pFBAcoeff)

            return(.Object)

          }
)
