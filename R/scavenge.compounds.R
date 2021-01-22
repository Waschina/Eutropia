# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Scavenge compounds   #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="scavenge.compounds",
           def=function(object, scavengeDist, cellSize, localEnv)
           {
             standardGeneric("scavenge.compounds")
           }
)

#' Scavenge available compounds from local environment
#'
#' @param object A \code{growthSimulation} object.
#' @param localEnv List for local environment. The two expected list elements are: \code{hex.id} is a 
#' vector if hexagon indeces of the environment grid and \code{hex.dist} the distance to each of the respective hexagons
#' within the local environment.
#' @return Lists with two elements: \code{compounds}: character vector of compound ids, and \code{fmol}: numeric
#' vector with the absolute metabolite availability in fmol.
setMethod(f          = "scavenge.compounds",
          signature  = signature(object         = "growthSimulation", 
                                 scavengeDist   = "numeric",
                                 cellSize       = "numeric",
                                 localEnv       = "list"),
          definition = function(object, scavengeDist, cellSize, localEnv) {
            
            # get accessible portion of environment grids based on distance
            accPortion <-  - 1 / scavengeDist * (localEnv$hex.dist - (cellSize/2 + scavengeDist))

            # -  (localEnv$hex.dist -  object@par_scavengeRadius)*(1 / (object@par_scavengeRadius - object@par_cellDiameter/2))
            accPortion <- ifelse(accPortion > 1, 1, accPortion)
            
            # retrieve absolute amount of accessible metabolites
            # TODO: Add constant compounds
            accCpd <- list()
            tmp_met <- object@environ@concentrations[localEnv$hex.id,] / 1000 * object@environ@hexVol # fmol in each hexagon (div. by 1000 because conc. are stored in mM)
            tmp_met <- tmp_met * accPortion # fmol accessible in each hexagon
            
            accCpd$fmolPerHex <- tmp_met
            
            tmp_met <- apply(tmp_met,2,sum) # total fmol accible to cell
              
            accCpd$compounds <- object@environ@compounds
            accCpd$fmol      <- tmp_met
            
            accCpd$hex.id   <- localEnv$hex.id
            accCpd$hex.dist <- localEnv$hex.dist
            
            
            return(accCpd)
          }
)