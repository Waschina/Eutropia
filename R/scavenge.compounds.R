#' # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
#' # Scavenge compounds   #
#' # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
#' setGeneric(name="scavenge.compounds",
#'            def=function(object, localEnv)
#'            {
#'              standardGeneric("scavenge.compounds")
#'            }
#' )
#'
#' #' Scavenge available compounds from local environment
#' #'
#' #' @param object A \code{growthSimulation} object.
#' #' @param localEnv data.table for local environment. The expected columns elements are: \code{field.id} is a
#' #' vector if field indices of the environment grid; \code{field.dist} the distance to each of the respective fields
#' #' within the local environment; \code{acc.prop} specifies the accessible proportion of the field to the cell.
#' #' @return Lists with the following elements: \code{compounds}: character vector of compound ids, and \code{fmol}: numeric
#' #' vector with the absolute metabolite availability in fmol.
#' setMethod(f          = "scavenge.compounds",
#'           signature  = signature(object         = "growthSimulation",
#'                                  localEnv       = "list"),
#'           definition = function(object, localEnv) {
#'
#'             # retrieve absolute amount of accessible metabolites
#'             # TODO: Add constant compounds
#'             accCpd <- list()
#'             tmp_met <- object@environ@concentrations[localEnv$hex.id,] / 1000 * object@environ@hexVol # fmol in each hexagon (div. by 1000 because conc. are stored in mM and hex volume in µm^3)
#'             tmp_met <- tmp_met * localEnv$acc.prop # fmol accessible in each hexagon
#'
#'             accCpd$fmolPerHex <- tmp_met
#'
#'             tmp_met <- apply(tmp_met,2,sum) # total fmol accessible to cell
#'
#'             accCpd$compounds <- object@environ@compounds
#'             accCpd$fmol      <- tmp_met
#'
#'             accCpd$hex.id   <- localEnv$hex.id
#'             accCpd$hex.dist <- localEnv$hex.dist
#'
#'             return(accCpd)
#'           }
#' )
