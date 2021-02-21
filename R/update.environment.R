# # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# # Scavenge compounds   #
# # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# setGeneric(name="update.environment",
#            def=function(object, model, cMass, accCpdFMOL, sol)
#            {
#              standardGeneric("update.environment")
#            }
# )
#
# setMethod(f          = "update.environment",
#           signature  = signature(object         = "growthSimulation",
#                                  model          = "character",
#                                  cMass          = "numeric",
#                                  accCpdFMOL     = "list",
#                                  sol            = "list"),
#           definition = function(object, model, cMass, accCpdFMOL, sol) {
#             if(sol$stat != ok | sol$obj <= 0) {
#               return(object@environ)
#             }
#
#             # accessibility matrix
#             accmat <- as.matrix(accCpdFMOL$fmolPerHex)
#             colnames(accmat) <- accCpdFMOL$compounds
#             accmat.row2hex   <- accCpdFMOL$hex.id
#             accmat <- t(t(accmat)/colSums(accmat)) # relative accessibilities
#
#             # get exchange reactions with non-zero flux
#             ex.ind <- grep("^EX_",object@models[[model]]@mod@react_id)
#             ex.flx <- sol$fluxes[ex.ind]
#             names(ex.flx) <- object@models[[model]]@mod@react_id[ex.ind]
#             ex.flx <- ex.flx[abs(ex.flx) > 0]
#
#             # normalise to time and cell mass
#             ex.flx <- ex.flx * cMass * object@deltaTime # uptake / production in abolute fmol in this time step and by this cell of its specific mass
#
#
#             # -
#             # Uptake
#             # -
#             ex.upt <- ex.flx[ex.flx < 0]
#             names(ex.upt) <- gsub("^EX_","", names(ex.upt))
#             ex.upt <- ex.upt[names(ex.upt) %in% colnames(accmat)]
#             if(length(ex.upt) > 0) {
#
#               accmat.upt <- accmat[, names(ex.upt), drop = F]
#               accmat.upt <- t(t(accmat.upt) * ex.upt) # fmol taken up in from each hexagon
#               accmat.upt <- (accmat.upt / 1e12)  / (object@environ@hexVol / 1e15 ) # mM taken up from each hexagon
#
#               concMatInd <- match(colnames(accmat.upt), object@environ@compounds)
#
#               I <- matrix(c(rep(accmat.row2hex, ncol(accmat.upt)),
#                             rep(concMatInd, each =nrow(accmat.upt)),
#                             as.numeric(accmat.upt)),
#                           ncol = 3, byrow = F)
#
#               object@environ@concentrations[I[,1:2]] <- object@environ@concentrations[I[,1:2]] + I[,3]
#               object@environ@concentrations[I[,1:2]] <- ifelse(object@environ@concentrations[I[,1:2]] < 1e-7, 0, object@environ@concentrations[I[,1:2]])
#
#             }
#
#             # -
#             # Production
#             # -
#             ex.pro <- ex.flx[ex.flx > 0]
#             names(ex.pro) <- gsub("^EX_","", names(ex.pro))
#             ex.pro <- ex.pro[names(ex.pro) %in% object@environ@compounds]
#             if(length(ex.pro) > 0) {
#               #hex.frac <- 1 - accCpdFMOL$hex.dist # WRONG
#               hex.frac <- max(accCpdFMOL$hex.dist) - accCpdFMOL$hex.dist
#               hex.frac <- hex.frac/sum(hex.frac)
#
#               ex.pro <- (ex.pro / 1e12) / (object@environ@hexVol / 1e15) # mM produced if all given to a single hexagon
#
#               pro.mat <- hex.frac %*% t(ex.pro)
#
#               concMatInd <- match(colnames(pro.mat), object@environ@compounds)
#
#               I <- matrix(c(rep(accmat.row2hex, ncol(pro.mat)),
#                             rep(concMatInd, each =nrow(pro.mat)),
#                             as.numeric(pro.mat)),
#                           ncol = 3, byrow = F)
#
#               object@environ@concentrations[I[,1:2]] <- object@environ@concentrations[I[,1:2]] + I[,3]
#             }
#
#
#             return(object@environ)
#           }
# )
