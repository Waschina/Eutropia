# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Scavenge compounds   #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

# Scavenge available compounds from local environment
#
# @param object A \code{growthSimulation} object.
# @param localEnv data.table for local environment. The expected columns elements are: \code{field.id} is a
# vector if field indeces of the environment grid; \code{field.dist} the distance to each of the respective fields
# within the local environment; \code{acc.prop} specifies the accessible proportion of the field to the cell.
# @return Lists with the following elements: \code{compounds}: character vector of compound ids, and \code{fmol}: numeric
scavenge_compounds <- function(localEnv, env_conc, env_cpds, env_fieldVol) {

  # retrieve absolute amount of accessible metabolites
  # TODO: Add constant compounds
  accCpd <- list()
  tmp_met <- env_conc / 1000 * env_fieldVol # fmol in each field (div. by 1000 because conc. are stored in mM and field volume in micro-m^3)
  tmp_met <- tmp_met * localEnv$acc.prop # fmol accessible in each field

  accCpd$fmolPerField <- tmp_met

  tmp_met <- apply(tmp_met,2,sum) # total fmol accessible to cell

  accCpd$compounds <- env_cpds
  accCpd$fmol      <- tmp_met

  accCpd$field.id   <- localEnv$field.id
  accCpd$field.dist <- localEnv$field.dist

  return(accCpd)
}
