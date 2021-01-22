# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Run simulation       #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="run.simulation",
           def=function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300, record = NULL, ...)
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
#'
#' @details TODO.
setMethod(f          = "run.simulation",
          signature  = signature(object    = "growthSimulation",
                                 niter     = "numeric"),
          definition = function(object, niter, verbose = 1, lim_cells = 1e5, lim_time = 300, record = NULL) {

            t_start <- Sys.time()
            j <- 1

            while(j <= niter & difftime(Sys.time(), t_start, units = "mins") < lim_time & nrow(object@cellDT) < lim_cells) {
              simRound <- object@n_rounds + 1

              # get the cells' info
              cellDT <- copy(object@cellDT)
              ncells <- nrow(cellDT)
              smass  <- sum(cellDT$mass)

              # get elapsed time
              elapT <- round(difftime(Sys.time(), t_start, units = "secs"))
              elapT <- as.character(as_hms(elapT))

              if(verbose > 0)
                cat("[",elapT,"] Simulation round ",simRound," \t(",ncells," cells, ",format(round(smass, digits = 2), nsmall = 2)," pg dBM)\n", sep ='')


              # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
              # (1) get each cell's neighboring environment hexagonal grid cells  #
              # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... getting cells' local environment\n", sep ='')

              # get grid hexagon positions
              gridDT <- as.data.table(object@environ@hex.pts)
              gridDT[, hex := 1:.N]

              cellHexEnv <- apply(cellDT, 1, function(cDT) {
                gtmp <- gridDT[abs(x - as.numeric(cDT["x"])) <= (as.numeric(cDT["size"])/2 + object@models[[cDT["type"]]]@scavengeDist) &
                               abs(y - as.numeric(cDT["y"])) <= (as.numeric(cDT["size"])/2 + object@models[[cDT["type"]]]@scavengeDist) ]
                gsp <- SpatialPoints(gtmp[,.(x,y)])

                csp <- SpatialPoints(matrix(c(as.numeric(cDT["x"]),as.numeric(cDT["y"])), ncol = 2))

                qry.env  <- which(gWithinDistance(csp,
                                                  gsp,
                                                  dist = as.numeric(cDT["size"])/2 + object@models[[cDT["type"]]]@scavengeDist, byid = T),
                                  useNames = F)
                qry.dist <- gDistance(gsp[qry.env,], csp, byid = T)[1,]

                return(list(hex.id   = gtmp[qry.env, hex],
                            hex.dist = qry.dist))
              })


              # - - - - - - - - - - - - #
              # (2) iterative agentFBA  #
              # - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... performing cell-agent-FBA\n", sep ='')

              # random order to iterate through cell list
              cell_order <- sample(1:nrow(cellDT), nrow(cellDT))

              for(i in cell_order) {
                if(any(object@environ@concentrations < 0)) {
                  saveRDS(object, "object.RDS")
                  stop("neg. conc")
                }
                # (2.0) Get local accessible compounds
                cMass <- cellDT[i, mass]
                accCompounds <- scavenge.compounds(object,
                                                   scavengeDist = object@models[[cellDT[i, type]]]@scavengeDist,
                                                   cellSize = cellDT[i, size],
                                                   localEnv = cellHexEnv[[i]])
                updatedEX    <- adjust.uptake(object,
                                              model = cellDT[i, type],
                                              cMass = cMass,
                                              accCpdFMOL = accCompounds)

                # (2.1) agentFBA(model, envGrids, curMass) for independent cells

                # constrain model to lcoal environment
                changeColsBnds(problem(object@models[[cellDT[i, type]]]@mod.warm), updatedEX$ex.react.ind,
                               lb = updatedEX$ex.react.lb, ub = object@models[[cellDT[i, type]]]@mod@uppbnd[updatedEX$ex.react.ind])

                sol.fba <- optimizeProb(object@models[[cellDT[i, type]]]@mod.warm)
                mu <- sol.fba$obj * object@deltaTime
                #
                # if(sol.fba$stat == 3) {
                #   print(cellDT[i])
                #   print(cellHexEnv[[i]])
                #   print(accCompounds)
                #   print(updatedEX)
                #   stop("nanu")
                # }

                # update cell mass and size
                if(sol.fba$obj > 0 & sol.fba$stat == ok) {
                  cellDT[cell == i, mass := mass * (1+mu)]
                  cellDT[cell == i, size := size * (1+mu)^(1/3)] # sphere
                }

                # restore former bounds
                changeColsBnds(problem(object@models[[cellDT[i, type]]]@mod.warm), updatedEX$ex.react.ind,
                               lb = object@models[[cellDT[i, type]]]@mod@lowbnd[updatedEX$ex.react.ind],
                               ub = object@models[[cellDT[i, type]]]@mod@uppbnd[updatedEX$ex.react.ind])


                # (2.2) Update environment
                if(sol.fba$obj > 0 & sol.fba$stat == ok) {
                  object@environ <- update.environment(object, cellDT[i, type], cMass, accCompounds, sol.fba)
                }

                # (2.3) Back to (2.1) with next cell

              }

              # - - - - - - - - - - - - #
              # (3) compound diffusion  #
              # - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... diffusion of compounds\n", sep ='')

              object@environ <- diffuse.compounds(object@environ, n_iter = object@diffusionNIter)

              # - - - - - - - - - #
              # (4) devide cells  #
              # - - - - - - - - - #
              if(verbose > 1)
                cat("... binary fission of large cells\n", sep ='')

              cmad <- unlist(lapply(object@models, function(x) x@cellMassAtDivision))
              cellDT[, cellMassAtDivision := cmad[type]]

              gind <- which(cellDT$mass >= cellDT$cellMassAtDivision)
              if(length(gind) > 0) {
                newCells <- copy(cellDT[gind])
                newCells[, parent := cell] # save parent information
                newCells[, cell := 1:.N]
                newCells[, cell := cell + ncells]

                # update cell mass and size of divided cells
                cellDT <- rbind(cellDT, newCells)
                cellDT[mass >= cellMassAtDivision, size := size * 0.5^(1/3)]
                cellDT[mass >= cellMassAtDivision, mass := mass / 2]

                # place daughter cell next to parent cell (random angle/direction)
                cellDT[cell > ncells, x := x + size * cos(runif(.N, max = 2*pi))]
                cellDT[cell > ncells, y := y + size * sin(runif(.N, max = 2*pi))]

                # invert velocity of daughter cells
                cellDT[cell > ncells, x.vel := x.vel * (-1)]
                cellDT[cell > ncells, y.vel := y.vel * (-1)]
              }

              cellDT[, cellMassAtDivision := NULL]

              object@cellDT <- copy(cellDT)

              ncells_new <- nrow(cellDT)

              cvmax <- unlist(lapply(object@models, function(x) x@vmax))

              # set up particle simulation
              sim_new <- tbl_graph(nodes = cellDT, directed = F, node_key = "cell") %>%
                simulate(setup = predefined_genesis(x = cellDT$x,
                                                    y = cellDT$y,
                                                    x_vel = cellDT$x.vel,
                                                    y_vel = cellDT$y.vel)) %>%
                wield(random_force,
                      xmin = -object@rMotion,
                      xmax =  object@rMotion,
                      ymin = -object@rMotion,
                      ymax =  object@rMotion) %>%
                wield(collision_force, radius = cellDT$size/2, n_iter = 50, strength = 0.7) %>%
                impose(velocity_constraint,
                       vmax = cvmax[cellDT$type]) %>%
                impose(polygon_constraint,
                       polygon = object@universePolygon)

              saveRDS(sim_new, file ="sim_new.RDS")
              # - - - - - - - - - - - - - - - - - - - - - - - - #
              # (5) Move & collide / re-organise cells in space #
              # - - - - - - - - - - - - - - - - - - - - - - - - #
              if(verbose > 1)
                cat("... sliding and colliding cells\n", sep ='')
              sim_new <- evolve(sim_new, steps = 1)

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
                    #COI.rec <- cbind(data.table(gridID = 1:object@environ@nhex), COI.rec)

                    object@history[[simRound]]$compound <- COI.rec
                  }
                }

                # TODO: Record individual fluxes
              }


              j <- j + 1
              object@n_rounds <- object@n_rounds + 1
            }

            return(object)
          }
)
