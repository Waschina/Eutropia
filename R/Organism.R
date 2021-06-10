setClass("Organism",

         slots = c(
           # biological parameters
           cellDiameter       = "numeric", # µm
           cellMassInit       = "numeric", # pg
           cellMassAtDivision = "numeric", # pg
           cellShape          = "character", # currently only coccus
           vmax               = "numeric", # in µm/s
           scavengeDist       = "numeric", # in µm
           chemotaxisCompound = "character",
           chemotaxisStrength = "numeric", # between 1 (attracting) and -1 (repelling)

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
                   chemotaxisCompound,
                   chemotaxisStrength,
                   open.bounds,
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

            # Rm exchange reaction for D-lactate if L-lactate is also present
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
