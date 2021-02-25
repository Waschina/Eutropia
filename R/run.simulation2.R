# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Run simulation       #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="run.simulation",
           def=function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300, record = NULL, n.cores = NULL, ...)
           {
             standardGeneric("run.simulation")
           },
           signature = c("object","niter")
)

#' Run simulation
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
setMethod(f          = "run.simulation",
          signature  = signature(object    = "growthSimulation",
                                 niter     = "numeric"),
          definition = function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300, record = NULL, n.cores = NULL) {

            # initialise multi core processing (with a copy of each model for each parellel fork)
            cmad <- unlist(lapply(object@models, function(x) x@cellMassAtDivision))
            ok <- 1
            if(is.null(n.cores))
              n.cores <- min(10, detectCores()-1)
            fork_models       <- object@models
            fork_deltaTime    <- object@deltaTime
            fork_fieldVol     <- object@environ@fieldVol
            fork_envCompounds <- object@environ@compounds
            cat("Initalising simulations using",n.cores,"CPU cores...\n")
            cl <- makeCluster(n.cores)

            init_forks <- function(x){
              require(EcoAgents)
              SYBIL_SETTINGS("SOLVER","cplexAPI")
              sybil::SYBIL_SETTINGS("METHOD", "hybbaropt") # Experimental
              #tid <<- x
              fork_mods <- list()
              for(mi in names(fork_models)) {
                fork_mods[[mi]] <- sysBiolAlg(fork_models[[mi]]@mod,
                                               algorithm = "mtf2",
                                               pFBAcoeff = 1e-6)
              }
              fork_mods <<- fork_mods
              #mod_env <<- environment()
              #return(fork_mods)
              return(NULL)
            }
            fork_ids <- clusterApply(cl, 1:n.cores, init_forks)
            rm(fork_ids)

            # keeping an eye on time
            t_start <- Sys.time()
            j <- 1

            # get grid field positions as data.table
            gridDT <- as.data.table(object@environ@field.pts)
            gridDT[, field := 1:.N]

            # get organims' scvenge radius
            get_scv_radius <- function(ctype) {
              unlist(lapply(ctype, function(x) object@models[[x]]@scavengeDist))
            }

            #ram_usage <<- c(mem_used())
            # while loop that checks termination criteria before starting a new round
            while(j <= niter & difftime(Sys.time(), t_start, units = "mins") < lim_time & nrow(object@cellDT) < lim_cells) {
              simRound <- object@n_rounds + 1

              # # give workers the current metabolite concentrations
              fork_conc <- object@environ@concentrations
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
              object@cellDT[, tmp_scv_radius := get_scv_radius(type)]

              get_cellFieldEnv <- function(x, y, size, tmp_scv_radius, gridLC) {
                icell_x    <- x
                icell_y    <- y
                icell_size <- size
                scv_dist   <- tmp_scv_radius

                # pre-filter grid fields
                # gtmp <- gridDT[abs(x - icell_x) <= (icell_size/2 + scv_dist) &
                #                  abs(y - icell_y) <= (icell_size/2 + scv_dist)]
                # gsp  <- SpatialPoints(gtmp[,.(x,y)])
                gsp  <- SpatialPoints(gridLC[,.(x,y,z)])

                csp <- SpatialPoints(matrix(c(icell_x, icell_y, 0), ncol = 3))

                qry.env  <- which(gWithinDistance(csp,
                                                  gsp,
                                                  dist = icell_size/2 + scv_dist, byid = T),
                                  useNames = F)
                qry.dist <- gDistance(gsp[qry.env,], csp, byid = T)[1,]

                # get accessible portion of environment grids based on distance
                accPortion <-  - 1 / scv_dist * (qry.dist - (icell_size/2 + scv_dist))
                accPortion <- ifelse(accPortion > 1, 1, accPortion)

                return(data.table(field.id   = gridLC[qry.env, field],
                                  field.dist = qry.dist,
                                  acc.prop   = accPortion))
              }

              fork_cellDT <- copy(object@cellDT)
              # calculate the local environmental conditions. Use parallel version only if more than 500 cell
              if(ncells < 500) {
                cellFieldEnv <- lapply(1:ncells, function(i) get_cellFieldEnv(x = fork_cellDT[i, x],
                                                                              y = fork_cellDT[i, y],
                                                                              size = fork_cellDT[i, size],
                                                                              tmp_scv_radius = fork_cellDT[i, tmp_scv_radius],
                                                                              gridLC = gridDT[abs(x - fork_cellDT[i, x]) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius]) &
                                                                                              abs(y - fork_cellDT[i, y]) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius]) &
                                                                                              abs(z - 0                  ) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius])]))

              } else {
                cellFieldEnv <- parLapply(cl, 1:ncells, function(i) get_cellFieldEnv(x = fork_cellDT[i, x],
                                                                                     y = fork_cellDT[i, y],
                                                                                     size = fork_cellDT[i, size],
                                                                                     tmp_scv_radius = fork_cellDT[i, tmp_scv_radius],
                                                                                     gridLC = gridDT[abs(x - fork_cellDT[i, x]) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius]) &
                                                                                                     abs(y - fork_cellDT[i, y]) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius]) &
                                                                                                     abs(z - 0                  ) <= (fork_cellDT[i, size]/2 + fork_cellDT[i, tmp_scv_radius])]))
              }

              object@cellDT[, tmp_scv_radius := NULL] # a clean up

              # close cells may claim in sum more than 100% of the resources from a field
              # Proportianally scale down accessibility in those cases:
              cellFieldEnv <- rbindlist(cellFieldEnv, idcol = T)
              cellFieldEnv[, field_claim_sum := sum(acc.prop), by = field.id]
              cellFieldEnv[ field_claim_sum > 1, acc.prop := acc.prop / field_claim_sum]

              # - - - - - - - - - - - - #
              # (2) parallel agentFBA   #
              # - - - - - - - - - - - - #K
              if(verbose > 1)
                cat("... performing cell-agent-FBA\n", sep ='')

              # the core part. Lets do this
              agentFBA <- function(localEnv, cMass, size, type, env_conc, env_cpds) {

                # (2.0) Get local accessible compounds
                accCompounds <- scavenge.compounds(localEnv     = localEnv,
                                                   env_conc     = env_conc,
                                                   env_cpds     = env_cpds,
                                                   env_fieldVol = fork_fieldVol)
                updatedEX    <- adjust.uptake(model      = fork_models[[type]]@mod,
                                              cMass      = cMass,
                                              accCpdFMOL = accCompounds,
                                              deltaTime  = fork_deltaTime)

                # (2.1) agentFBA(model, envGrids, curMass) for independent cells
                # constrain model to lcoal environment
                ccbnds <- changeColsBnds(problem(fork_mods[[type]]), updatedEX$ex.react.ind,
                                         lb = updatedEX$ex.react.lb, ub = fork_models[[type]]@mod@uppbnd[updatedEX$ex.react.ind])

                sol.fba <- optimizeProb(fork_mods[[type]])
                mu <- sol.fba$obj * fork_deltaTime
                stat <- sol.fba$stat

                # restore former bounds
                ccbnds <- changeColsBnds(problem(fork_mods[[type]]), updatedEX$ex.react.ind,
                                         lb = fork_models[[type]]@mod@lowbnd[updatedEX$ex.react.ind],
                                         ub = fork_models[[type]]@mod@uppbnd[updatedEX$ex.react.ind])

                # update cell mass and size
                if(sol.fba$obj > 0 & sol.fba$stat == ok) {
                  cMass_new <- cMass * (1+mu)
                  cSize_new <- size * (1+mu)^(1/3) # sphere
                } else {
                  cMass_new <- cMass
                  cSize_new <- size
                }

                # get fluxes
                flx <- sol.fba$fluxes[1:fork_models[[type]]@mod@react_num]

                # get exchange reactions with non-zero flux
                # adjusting them to time and cell mass
                ex.ind <- grep("^EX_",fork_models[[type]]@mod@react_id)
                ex.flx <- sol.fba$fluxes[ex.ind]
                names(ex.flx) <- fork_models[[type]]@mod@react_id[ex.ind]
                ex.flx <- ex.flx[abs(ex.flx) > 0]
                # normalise to time and cell mass
                ex.flx <- ex.flx * cMass * fork_deltaTime # uptake / production in abolute fmol in this time step and by this cell of its specific mass

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
                  accmat.upt <- (accmat.upt / 1e12)  / (fork_fieldVol / 1e15 ) # mM taken up from each field

                  concMatInd <- match(colnames(accmat.upt), fork_envCompounds)

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
                ex.pro <- ex.pro[names(ex.pro) %in% fork_envCompounds]
                if(length(ex.pro) > 0) {
                  field.frac <- max(accCompounds$field.dist) - accCompounds$field.dist
                  field.frac <- field.frac/sum(field.frac)

                  ex.pro <- (ex.pro / 1e12) / (fork_fieldVol / 1e15) # mM produced if all given to a single field

                  pro.mat <- field.frac %*% t(ex.pro)

                  concMatInd <- match(colnames(pro.mat), fork_envCompounds)

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
                #gc() # das bremst, passt aber auf, dass kein RAM leakt
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

              agFBA_results <- parLapply(cl, 1:ncells, function(x) agentFBA(cellFieldEnv[.id == x],
                                                                            fork_cellDT[x, mass],
                                                                            fork_cellDT[x, size],
                                                                            fork_cellDT[x, type],
                                                                            env_conc = fork_conc[cellFieldEnv[.id == x, field.id],],
                                                                            env_cpds = fork_envCompounds))

              # update the environment
              if(verbose > 1)
                cat("... update environment\n", sep ='')
              I <- rbindlist(lapply(agFBA_results, function(x) rbindlist(list(x$I.up, x$I.pd))))
              I <- I[, .(concChange = sum(concChange)), by = c("field.id", "concMatInd")]
              I <- as.matrix(I)

              object@environ@concentrations[I[,1:2]] <- object@environ@concentrations[I[,1:2]] + I[,3]
              object@environ@concentrations[I[,1:2]] <- ifelse(object@environ@concentrations[I[,1:2]] < 1e-7, 0, object@environ@concentrations[I[,1:2]])

              # - - - - - - - - - - - - #
              # (3) compound diffusion  #
              # - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... diffusion of compounds\n", sep ='')

              object@environ <- diffuse.compounds(object@environ, n_iter = object@diffusionNIter)


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

                # specific compound contentrations per grids
                if(any(grepl("^compound_", record)) | ("compounds" %in% record)) {
                  COI <- record[grepl("^compound_", record)] # metabolies of interest
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


              j <- j + 1
              #ram_usage <<- c(ram_usage, mem_used())
              object@n_rounds <- object@n_rounds + 1
              #rm(agFBA_results)
              #gc()
              #clusterApply(cl, 1:n.cores, gc)
            }

            # delete cplex problem object to free memory
            fork_ids <- clusterApply(cl, 1:n.cores, function(x){
              for(mi in names(fork_mods)) {
                delProbCPLEX(fork_mods[[mi]]@problem@oobj@env, fork_mods[[mi]]@problem@oobj@lp)
                closeEnvCPLEX(fork_mods[[mi]]@problem@oobj@env)
              }
            })

            # stop cluster forks
            stopCluster(cl)

            return(object)
          }
)
