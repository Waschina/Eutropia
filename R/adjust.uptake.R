# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Update models LB                 #
# based on compound accessibility  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="adjust.uptake",
           def=function(object, model, cMass, accCpdFMOL)
           {
             standardGeneric("adjust.uptake")
           }
)

#' Adjust a model's lower bounds for exchange reactions based on compound accessibility
#'
setMethod(f          = "adjust.uptake",
          signature  = signature(object         = "growthSimulation", 
                                 model          = "character",
                                 cMass          = "numeric",
                                 accCpdFMOL     = "list"),
          definition = function(object, model, cMass, accCpdFMOL) {
            
            # get overlap of models exchange reactions with lb < 0 and metablites in environment
            accMets <- paste0("EX_",accCpdFMOL$compounds)
            
            rel.mets.inds <- which(accMets %in% object@models[[model]]@mod@react_id)
            if(length(rel.mets.inds) == 0)
              return(list(ex.react.ind = NULL,
                          ex.react.lb  = NULL))
            
            rel.mets <- accMets[rel.mets.inds]
            model.ex.ind <- which(object@models[[model]]@mod@react_id %in% rel.mets)
            names(model.ex.ind) <- object@models[[model]]@mod@react_id[model.ex.ind]
            model.ex.ind <- model.ex.ind[rel.mets]
            
            lb.mat <- matrix(0, nrow = length(rel.mets), ncol = 2)
            
            lb.mat[,1] <- accCpdFMOL$fmol[rel.mets.inds] * (1 / object@deltaTime) / cMass
            lb.mat[,2] <- -object@models[[model]]@mod@lowbnd[model.ex.ind]
            
            ex.react.lb <- apply(lb.mat,1,min)
            
            return(list(ex.react.ind = model.ex.ind,
                        ex.react.lb  = -ex.react.lb))
          }
)