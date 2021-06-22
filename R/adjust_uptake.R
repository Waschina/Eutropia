# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Update models LB                  #
# based on compound accessibility & #
# availability                      #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# Adjust a model's lower bounds for exchange reactions based on compound accessibility
#
adjust_uptake <- function(model, cMass, accCpdFMOL, deltaTime) {

  # get overlap of models exchange reactions with lb < 0 and metabolites in environment
  accMets <- paste0("EX_",accCpdFMOL$compounds)

  rel.mets.inds <- which(accMets %in% model@react_id)
  if(length(rel.mets.inds) == 0)
    return(list(ex.react.ind = NULL,
                ex.react.lb  = NULL))

  rel.mets <- accMets[rel.mets.inds]
  model.ex.ind <- which(model@react_id %in% rel.mets)
  names(model.ex.ind) <- model@react_id[model.ex.ind]
  model.ex.ind <- model.ex.ind[rel.mets]

  lb.mat <- matrix(0, nrow = length(rel.mets), ncol = 2)

  lb.mat[,1] <- accCpdFMOL$fmol[rel.mets.inds] * (1 / deltaTime) / cMass # fmol/(pg * hr) = mmol/ (g * hr)
  lb.mat[,2] <- -model@lowbnd[model.ex.ind]

  ex.react.lb <- apply(lb.mat,1,min)

  return(list(ex.react.ind = model.ex.ind,
              ex.react.lb  = -ex.react.lb))
}
