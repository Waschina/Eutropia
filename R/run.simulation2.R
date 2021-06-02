# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Run simulation       #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

#' @title Run simulation
#'
#' @description Run the agent- and FBA-based simulation.
#'
#' @param object An object of class \code{growthSimulation}
#' @param niter Number of rounds to simulate
#' @param verbose Control the 'chattiness' of the simulation logs. 0 - no logs, 1 - only main logs, 2 - more lines per round.
#' @param lim_cells Simulation terminates of total number of cells exceed this value.
#' @param lim_time Simulation terminated after the first iteration that finished after this time limit (in minutes).
#' @param record Character vector that indicated which simulation variables should be recorded after each simulation round. See Details.
#' @param n.cores Number of CPUs to use for parallelisation. If NULL (default), it will use the number of detectable cores minus 1 and in maximum 10 cores.
#'
#' @details TODO.
#'
#' @export
setGeneric(name="run.simulation",
           def=function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300,
                        record = c("cells","global_compounds"), n.cores = NULL, ...)
           {
             standardGeneric("run.simulation")
           },
           signature = c("object","niter")
)

setMethod(f          = "run.simulation",
          signature  = signature(object    = "growthSimulation",
                                 niter     = "numeric"),
          definition = function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300,
                                record = c("cells","global_compounds"), n.cores = NULL) {

            # initialise multi core processing (with a copy of each model in warm for each parallel fork)
            cmad <- unlist(lapply(object@models, function(x) x@cellMassAtDivision))
            #ok <- 1
            if(is.null(n.cores))
              n.cores <- min(10, detectCores()-1)
            #fork_models       <- object@models
            #fork_deltaTime    <- object@deltaTime
            #fork_fieldVol     <- object@environ@fieldVol
            #fork_envCompounds <- object@environ@compounds
            cat("Initalising simulations using",n.cores,"CPU cores...\n")
            cl <- makeCluster(n.cores)

            pre_mod_list <- lapply(1:n.cores, function(x) { return(object@models) })
            fork_ids <- clusterApply(cl, pre_mod_list, fun = init_warm_mods)
            rm(pre_mod_list)

            # keeping an eye on time
            t_start <- Sys.time()
            j <- 1

            # get grid field positions as data.table
            gridDT <- as.data.table(object@environ@field.pts)
            gridDT[, field := 1:.N]
            setkeyv(gridDT, c("z","x","y"))

            # get organims' scvenge radius
            get_scv_radius <- function(ctype) {
              unlist(lapply(ctype, function(x) object@models[[x]]@scavengeDist))
            }

            #ram_usage <<- c(mem_used())
            # while loop that checks termination criteria before starting a new round
            while(j <= niter & difftime(Sys.time(), t_start, units = "mins") < lim_time & nrow(object@cellDT) < lim_cells) {
              simRound <- object@n_rounds + 1

              # # give workers the current metabolite concentrations
              # fork_conc <- object@environ@concentrations
              # updateConc <- clusterApply(cl, 1:n.cores, function(fork_conc) {fieldConc <<- fork_conc; return(NULL)})
              # rm(updateConc)

              # get the cells' info
              #cellDT <- copy(object@cellDT)
              ncells <- nrow(object@cellDT)
              smass  <- sum(object@cellDT$mass)

              # get elapsed time
              elapT <- round(difftime(Sys.time(), t_start, units = "secs"))
              elapT <- as.character(as_hms(elapT))

              if(verbose > 0)
                cat("[",elapT,"] Simulation round ",simRound," \t(",ncells," cells, ",format(round(smass, digits = 2), nsmall = 2)," pg dBM)\n", sep ='')


              # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
              # (1) get each cell's neighboring environment grid fields           #
              # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... getting cells' local environment\n", sep ='')

              # the following gets a list for each cell, that has these elements:
              # -> First element 'field.id' is an vector of ids of neighboring fields
              # -> Second element 'field.dist' is the distance of the field centers to the cell's center
              # -> Third element: 'acc.prop' is the proportion of the respective field, claimed/accessible by the cell

              pre_celFieldEnv <- lapply(1:ncells, FUN = function(i) {
                res <- list()
                res[["x"]]      <- object@cellDT[i, x]
                res[["y"]]      <- object@cellDT[i, y]
                res[["size"]]   <- object@cellDT[i, size]
                res[["scvr"]]   <- get_scv_radius(object@cellDT[i, type])
                # TODO: Theres probably a smater way of doing the next command - but works for now
                res[["gridLC"]] <- gridDT[abs(x - res[["x"]]) <= (res[["size"]] /2 + res[["scvr"]]) &
                                            abs(y - res[["y"]]) <= (res[["size"]] /2 + res[["scvr"]]) &
                                            abs(z - 0                  ) <= (res[["size"]] /2 + res[["scvr"]])]
                return(res)

              })

              cellFieldEnv <- parLapply(cl, pre_celFieldEnv, get_cellFieldEnv_ex)

              # close cells may claim in sum more than 100% of the resources from a field
              # Proportianally scale down accessibility in those cases:
              cellFieldEnv <- rbindlist(cellFieldEnv, idcol = T)
              cellFieldEnv[, field_claim_sum := sum(acc.prop), by = field.id]
              cellFieldEnv[field_claim_sum > 1, acc.prop.tmp := acc.prop^exp(1)]
              cellFieldEnv[field_claim_sum > 1, acc.prop.tmp := acc.prop.tmp / sum(acc.prop.tmp) * field_claim_sum, by = field.id]
              cellFieldEnv[field_claim_sum > 1, acc.prop := acc.prop.tmp / field_claim_sum]

              cellFieldEnv <- lapply(1:ncells, FUN = function(x) {
                res <- list()
                res[["field.id"]]   <- cellFieldEnv[.id == x, field.id]
                res[["field.dist"]] <- cellFieldEnv[.id == x, field.dist]
                res[["acc.prop"]]   <- cellFieldEnv[.id == x, acc.prop]

                res[["field.conc"]] <- object@environ@concentrations[res[["field.id"]],]
                res[["field.cpds"]] <- object@environ@compounds
                res[["fieldVol"]]   <- object@environ@fieldVol
                res[["model"]]      <- object@models[[object@cellDT[x, type]]]
                res[["deltaTime"]]  <- object@deltaTime

                res[["type"]]       <- object@cellDT[x, type]
                res[["cMass"]]      <- object@cellDT[x, mass]
                res[["size"]]       <- object@cellDT[x, size]

                return(res)
              })

              # - - - - - - - - - - - - #
              # (2) parallel agentFBA   #
              # - - - - - - - - - - - - #K
              if(verbose > 1)
                cat("... performing cell-agent-FBA\n", sep ='')

              agFBA_results <- parLapply(cl, cellFieldEnv, agentFBA_ex)

              # update the environment
              if(verbose > 1)
                cat("... update environment\n", sep ='')
              I <- rbindlist(lapply(agFBA_results, function(x) rbindlist(list(x$I.up, x$I.pd))))
              I <- I[, .(concChange = sum(concChange)), by = c("field.id", "concMatInd")]
              ind_var <- which(!object@environ@conc.isConstant)
              I <- I[concMatInd %in% ind_var]
              I <- as.matrix(I)

              I[,3] <- ifelse(is.na(I[,3]), 0, I[,3]) # TODO: Check where NA originate come from
              #print(I)

              object@environ@concentrations[I[,1:2]] <- object@environ@concentrations[I[,1:2]] + I[,3]
              object@environ@concentrations[I[,1:2]] <- ifelse(object@environ@concentrations[I[,1:2]] < 1e-7, 0, object@environ@concentrations[I[,1:2]])

              # - - - - - - - - - - - - #
              # (3) compound diffusion  #
              # - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... diffusion of compounds\n", sep ='')

              object@environ <- diffuse.compounds(object@environ,
                                                  n_iter = object@diffusionNIter,
                                                  cl = cl,
                                                  n.cores = n.cores)


              # - - - - - - - - - - - - - -#
              # (4) grow and devide cells  #
              # - - - - - - - - - - - - - -#
              if(verbose > 1)
                cat("... binary fission of large cells\n", sep ='')

              # grow cells
              object@cellDT$mass <- unlist(lapply(agFBA_results, function(x) x$cMass_new))
              object@cellDT$size <- unlist(lapply(agFBA_results, function(x) x$cSize_new))


              object@cellDT[, cellMassAtDivision := cmad[type]]

              gind <- which(object@cellDT$mass >= object@cellDT$cellMassAtDivision)
              if(length(gind) > 0) {
                newCells <- copy(object@cellDT[gind])
                newCells[, parent := cell] # save parent information
                newCells[, cell := 1:.N]
                newCells[, cell := cell + ncells]

                # update cell mass and size of divided cells
                object@cellDT <- rbind(object@cellDT, newCells)
                object@cellDT[mass >= cellMassAtDivision, size := size * 0.5^(1/3)]
                object@cellDT[mass >= cellMassAtDivision, mass := mass / 2]
                rm(newCells)

                # place daughter cell next to parent cell (random angle/direction)
                object@cellDT[cell > ncells, x := x + size * cos(runif(.N, max = 2*pi))]
                object@cellDT[cell > ncells, y := y + size * sin(runif(.N, max = 2*pi))]

                # invert velocity of daughter cells
                object@cellDT[cell > ncells, x.vel := x.vel * (-1)]
                object@cellDT[cell > ncells, y.vel := y.vel * (-1)]
              }

              object@cellDT[, cellMassAtDivision := NULL]

              ncells_new <- nrow(object@cellDT)

              cvmax <- unlist(lapply(object@models, function(x) x@vmax))

              # set up particle simulation
              sim_new <- tbl_graph(nodes = object@cellDT, directed = F, node_key = "cell") %>%
                simulate(setup = predefined_genesis(x = object@cellDT$x,
                                                    y = object@cellDT$y,
                                                    x_vel = object@cellDT$x.vel,
                                                    y_vel = object@cellDT$y.vel)) %>%
                wield(random_force,
                      xmin = -object@rMotion,
                      xmax =  object@rMotion,
                      ymin = -object@rMotion,
                      ymax =  object@rMotion) %>%
                wield(collision_force, radius = object@cellDT$size/2, n_iter = 50, strength = 0.7) %>%
                impose(velocity_constraint,
                       vmax = cvmax[object@cellDT$type]) %>%
                impose(polygon_constraint,
                       polygon = object@universePolygon)

              #saveRDS(sim_new, file ="sim_new.RDS")

              # - - - - - - - - - - - - - - - - - - - - - - - - #
              # (5) Move & collide / re-organise cells in space #
              # - - - - - - - - - - - - - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... sliding and colliding cells\n", sep ='')
              sim_new <- evolve(sim_new, steps = round(object@deltaTime * 60))

              object@cellDT$x <- sim_new$position[,1]
              object@cellDT$y <- sim_new$position[,2]
              object@cellDT$x.vel <- sim_new$velocity[,1]
              object@cellDT$y.vel <- sim_new$velocity[,2]

              # - - - - - - - - #
              # (6) Record data #
              # - - - - - - - - #
              if(!is.null(record)) {
                if(verbose > 1)
                  cat("... recording simulation data\n", sep ='')

                object@history[[simRound]] <- list()

                # cell positions
                if("cells" %in% record) {
                  object@history[[simRound]]$cells <- copy(object@cellDT)
                }

                # global metabolite concentrations (only variable)
                if("global_compounds" %in% record) {
                  ind_var <- !object@environ@conc.isConstant
                  dt_conc_tmp <- data.table(cpd.id = object@environ@compounds[ind_var],
                                            cpd.name = object@environ@compound.names[ind_var],
                                            global_concentration = apply(object@environ@concentrations[,ind_var],2,mean))
                  object@history[[simRound]]$global_compounds <- dt_conc_tmp
                }

                # specific compound concentrations per grids
                if(any(grepl("^compound_", record)) | ("compounds" %in% record)) {
                  COI <- record[grepl("^compound_", record)] # metabolites of interest
                  COI <- gsub("^compound_","",COI)
                  COI <- COI[COI %in% object@environ@compounds]
                  if("compounds" %in% record)
                    COI <- object@environ@compounds

                  if(length(COI > 0)) {
                    COI.ind <- match(COI, object@environ@compounds)

                    COI.rec <- object@environ@concentrations[,COI.ind]
                    colnames(COI.rec) <- object@environ@compounds[COI.ind]

                    object@history[[simRound]]$compound <- COI.rec
                  }
                }

                #gc()
                # TODO: Record individual fluxes
              }

              # small cell summary
              if(verbose > 0) {
                cellsum <- copy(object@cellDT[,.(mass = round(sum(mass), digits = 2)), by = type])
                cellsum[, tmp.sum := paste0(type,"(",mass,")")]
                cellsum <- paste(cellsum$tmp.sum, collapse = " ")
                cat("           ",cellsum,"\n", sep ='')
              }


              j <- j + 1

              object@n_rounds <- object@n_rounds + 1

            }

            # delete cplex problem object to free memory
            fork_ids <- clusterApply(cl, 1:n.cores, close_clusters)

            # stop cluster forks
            stopCluster(cl)

            return(object)
          }
)

#
#
#
#
#
init_warm_mods <- function(x) {
  require(EcoAgents)
  SYBIL_SETTINGS("SOLVER","cplexAPI")
  #sybil::SYBIL_SETTINGS("METHOD", "hybbaropt") # Experimental
  #tid <<- x
  fork_mods <- list()
  for(mi in names(x)) {
    fork_mods[[mi]] <- sysBiolAlg(x[[mi]]@mod,
                                  algorithm = "mtf2",
                                  pFBAcoeff = 1e-6)
  }
  fork_mods <<- fork_mods
  ok <<- 1
  #mod_env <<- environment()
  #return(fork_mods)
  return(NULL)
}

#
#
#
#
#
get_cellFieldEnv_ex <- function(x) {
  icell_x    <- x$x
  icell_y    <- x$y
  icell_size <- x$size
  scv_dist   <- x$scvr

  gsp  <- SpatialPoints(x$gridLC[,.(x,y,z)])

  csp <- SpatialPoints(matrix(c(icell_x, icell_y, 0), ncol = 3))

  # qry.env  <- which(gWithinDistance(csp,
  #                                   gsp,
  #                                   dist = icell_size/2 + scv_dist, byid = T),
  #                   useNames = F)
  q.c.dist <- spDists(gsp, csp)[,1]

  qry.env <- which(q.c.dist <= (icell_size/2 + scv_dist))

  #qry.dist <- gDistance(gsp[qry.env,], csp, byid = T)[1,]
  qry.dist <- q.c.dist[qry.env]

  # get accessible portion of environment grids based on distance
  accPortion <-  - 1 / scv_dist * (qry.dist - (icell_size/2 + scv_dist))
  accPortion <- ifelse(accPortion > 1, 1, accPortion)

  return(data.table(field.id   = x$gridLC[qry.env, field],
                    field.dist = qry.dist,
                    acc.prop   = accPortion))
}

#
#
#
#
#
agentFBA_ex <- function(x) {
  # (2.0) Get local accessible compounds
  accCompounds <- scavenge.compounds(localEnv     = x[c("acc.prop","field.id","field.dist")],
                                     env_conc     = x$field.conc,
                                     env_cpds     = x$field.cpds,
                                     env_fieldVol = x$fieldVol)
  updatedEX    <- adjust.uptake(model      = x$model@mod,
                                cMass      = x$cMass,
                                accCpdFMOL = accCompounds,
                                deltaTime  = x$deltaTime)

  # (2.1) agentFBA(model, envGrids, curMass) for independent cells
  # constrain model to lcoal environment
  ccbnds <- changeColsBnds(problem(fork_mods[[x$type]]), updatedEX$ex.react.ind,
                           lb = updatedEX$ex.react.lb, ub = x$model@mod@uppbnd[updatedEX$ex.react.ind])

  sol.fba <- optimizeProb(fork_mods[[x$type]])
  mu <- sol.fba$obj * x$deltaTime
  stat <- sol.fba$stat

  # restore former bounds
  ccbnds <- changeColsBnds(problem(fork_mods[[x$type]]), updatedEX$ex.react.ind,
                           lb = x$model@mod@lowbnd[updatedEX$ex.react.ind],
                           ub = x$model@mod@uppbnd[updatedEX$ex.react.ind])

  # update cell mass and size
  if(sol.fba$obj > 0 & sol.fba$stat == ok) {
    cMass_new <- x$cMass * (1+mu)
    cSize_new <- x$size * (1+mu)^(1/3) # sphere
  } else {
    cMass_new <- x$cMass
    cSize_new <- x$size
  }

  # get fluxes
  flx <- sol.fba$fluxes[1:x$model@mod@react_num]

  # get exchange reactions with non-zero flux
  # adjusting them to time and cell mass
  ex.ind <- grep("^EX_",x$model@mod@react_id)
  ex.flx <- sol.fba$fluxes[ex.ind]
  names(ex.flx) <- x$model@mod@react_id[ex.ind]
  ex.flx <- ex.flx[abs(ex.flx) > 0]
  # normalise to time and cell mass
  ex.flx <- ex.flx * x$cMass * x$deltaTime # uptake / production in absolute fmol in this time step and by this cell of its specific mass

  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  # contribution to environment (change in metabolite concentrations) #
  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
  # accessibility matrix:
  accmat           <- as.matrix(accCompounds$fmolPerField)
  colnames(accmat) <- accCompounds$compounds
  accmat.row2field <- accCompounds$field.id
  accmat           <- t(t(accmat)/colSums(accmat)) # relative accessibilities

  #--------#
  # Uptake #
  #--------#
  ex.upt <- ex.flx[ex.flx < 0]
  names(ex.upt) <- gsub("^EX_","", names(ex.upt))
  ex.upt <- ex.upt[names(ex.upt) %in% colnames(accmat)]
  if(length(ex.upt) > 0) {

    accmat.upt <- accmat[, names(ex.upt), drop = F]
    accmat.upt <- t(t(accmat.upt) * ex.upt) # fmol taken up in from each field
    accmat.upt <- (accmat.upt / 1e12)  / (x$fieldVol / 1e15 ) # mM taken up from each field

    concMatInd <- match(colnames(accmat.upt), x$field.cpds)

    I.up <- data.table(field.id = rep(accmat.row2field, ncol(accmat.upt)),
                       concMatInd = rep(concMatInd, each =nrow(accmat.upt)),
                       concChange = as.numeric(accmat.upt))

  } else {
    I.up <- data.table(field.id   = integer(0),
                       concMatInd = integer(0),
                       concChange = double(0))
  }

  #------------#
  # Production #
  #------------#
  ex.pro <- ex.flx[ex.flx > 0]
  names(ex.pro) <- gsub("^EX_","", names(ex.pro))
  ex.pro <- ex.pro[names(ex.pro) %in% x$field.cpds]
  if(length(ex.pro) > 0) {
    field.frac <- max(accCompounds$field.dist) - accCompounds$field.dist
    field.frac <- field.frac/sum(field.frac)

    ex.pro <- (ex.pro / 1e12) / (x$fieldVol / 1e15) # mM produced if all given to a single field

    pro.mat <- field.frac %*% t(ex.pro)

    concMatInd <- match(colnames(pro.mat), x$field.cpds)

    I.pd <- data.table(field.id   = rep(accmat.row2field, ncol(pro.mat)),
                       concMatInd = rep(concMatInd, each =nrow(pro.mat)),
                       concChange = as.numeric(pro.mat))

  } else {
    I.pd <- data.table(field.id   = integer(0),
                       concMatInd = integer(0),
                       concChange = double(0))
  }

  rm(ccbnds)
  rm(sol.fba)

  return(list(mu        = mu,
              fba.stat  = stat,
              cMass_new = cMass_new,
              cSize_new = cSize_new,
              fluxes    = flx,
              ex.fluxes = ex.flx,
              accmat    = accmat,
              I.up      = I.up,
              I.pd      = I.pd
  ))
}

close_clusters <- function(x){
  for(mi in names(fork_mods)) {
    delProbCPLEX(fork_mods[[mi]]@problem@oobj@env, fork_mods[[mi]]@problem@oobj@lp)
    closeEnvCPLEX(fork_mods[[mi]]@problem@oobj@env)
  }
}
