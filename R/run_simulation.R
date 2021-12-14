# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Run simulation       #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

#' @title Run simulation
#'
#' @description Run the agent- and FBA-based simulation.
#'
#' @param object An object of class \link{growthSimulation}
#' @param niter Number of rounds to simulate
#' @param verbose Control the 'chattiness' of the simulation logs. 0 - no logs,
#' 1 - only main logs, 2 - chaffinch.
#' @param lim_cells Simulation terminates of total number of cells exceed this
#' value.
#' @param lim_time Simulation terminated after the first iteration that finished
#' after this time limit (in minutes).
#' @param convergence.e Numeric indicating when community growth (pg) is
#' considered to have reached convergence
#' @param record Character vector that indicates, which simulation variables
#' should be recorded after each simulation iteration. See Details.
#' @param n.cores Number of CPUs to use for parallelisation. If NULL (default),
#' it will use the number of detectable cores minus 1 and in maximum 10 cores.
#' @param on.iteration An optional function that is performed at the end of each
#' iteration and returns again an object of class \link{growthSimulation}.
#' @param live.plot Logical. If TRUE, positions of cells are plotted at each
#' iteration. Caution: May reduce performance of simulation. Default: FALSE
#' @param ... Additional arguments to `on.iteration`
#'
#' @details
#' Recording: The cells' positions, masses, sizes, and metabolite exchanges are
#' always recorded. Also recorded are the global metabolite concentrations. Due
#' to memory considerations, local metabolite concentrations are not recorded
#' by default. However, users can specify which concentrations are tracked during
#' the simulation using the `record` option. E.g.
#'
#' \code{record = c("compound_cpd00029_e0","compound_cpd00211_e0")}
#'
#' records the concentration fo the to two metabolites `cpd00029_e0` and
#' `cpd00211_e0`. All compound concentration can be recorded with
#'
#' \code{record = "compounds"} ,
#'
#' yet, is not advised as this could be storage- and time-consuming. Compound
#' concentrations are not stored in memory but on the hard drive in a file
#' within the working directory. Tracked concentrations can also be plotted with
#' \link{plot_environment}.
#' Exoenzyme concentrations cannot be recorded in the current version.\cr
#'
#' Convergence is checked by calculating the ratio:
#' \deqn{c := | a / min(a_{i-1},...,a_{i-5}) - 1 |}
#' \eqn{a_i} is the total biomass at iteration \eqn{i}. The simulation
#' terminates if \eqn{c} is below `convergence.e`. Thus, one can expect
#' longer simulations when reducing `convergence.e`.
#'
#' @import lamW
#' @import parallel
#' @import hms
#' @import data.table
#' @import tidygraph
#' @import particles
#' @importFrom stats weighted.mean runif rnorm sd
#' @importFrom stringi stri_rand_strings
#'
#' @export
run_simulation <- function(object, niter, verbose = 1, lim_cells = 1e5,
                           lim_time = 300, convergence.e = 1e-4,
                           record = NULL,
                           n.cores = NULL,
                           live.plot = F,
                           on.iteration = NULL, ...) {

  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  # If something big is going to be recorded create a directory for recording files
  if(!is.null(record) & is.na(object@recordDir)) {
    object@recordDir <- paste0("recordings_",stringi::stri_rand_strings(1,6))
    if(!dir.exists(object@recordDir))
      dir.create(object@recordDir)
  }

  # save initial simulation status
  if(object@n_rounds == 0) {
    object@history[[1]] <- list()

    # cell positions
    object@history[[1]]$cells <- copy(object@cellDT)

    # global metabolite concentrations (only variable)
    ind_var <- !object@environ@conc.isConstant
    dt_conc_tmp <- data.table(cpd.id = object@environ@compounds[ind_var],
                              cpd.name = object@environ@compound.names[ind_var],
                              global_concentration = apply(object@environ@concentrations[,ind_var],2,mean))
    object@history[[1]]$global_compounds <- dt_conc_tmp

    # specific compound concentrations per grids
    if(!is.null(record) && (any(grepl("^compound_", record)) | ("compounds" %in% record))) {
      COI <- record[grepl("^compound_", record)] # metabolites of interest
      COI <- gsub("^compound_","",COI)
      COI <- COI[COI %in% object@environ@compounds]
      if("compounds" %in% record)
        COI <- object@environ@compounds

      if(length(COI) > 0) {
        rec_file <- paste0(object@recordDir,"/cpdrec_",1,"_0.RDS")
        k_tmp <- 1
        while(file.exists(rec_file)) {
          rec_file <- paste0(object@recordDir,"/cpdrec_",1,"_",k_tmp,".RDS")
          k_tmp <- k_tmp + 1
        }

        COI.ind <- match(COI, object@environ@compounds)
        COI.rec <- object@environ@concentrations[,COI.ind, drop = FALSE]
        COI.rec <- data.table(round(COI.rec, digits = 7))

        # save compound recording in recording file
        saveRDS(COI.rec, file = rec_file, compress = T)

        # save order of compounds
        object@history[[1]]$compounds <- object@environ@compounds[COI.ind]

        # save path to recording file
        object@history[[1]]$compounds.record <- rec_file
      }
    }
  }

  # devtools::check() does not recognize column names as variables and
  # throws warnings. This is prevented by assigned these names to NULL
  # https://github.com/Rdatatable/data.table/issues/850


  # initialize multi core processing (with a copy of each model in warm for each parallel fork)
  cmad <- unlist(lapply(object@models, function(x) x@cellMassAtDivision))

  if(is.null(n.cores))
    n.cores <- min(10, detectCores()-1)

  if(verbose >= 1)
    cat("Initalising simulations using",n.cores,"CPU cores...\n")
  cl <- makeCluster(n.cores)

  # Check if cplexAPI is installed
  lpsolver <- "glpkAPI"
  okcode   <- 5
  if("cplexAPI" %in% rownames(utils::installed.packages())) {
    lpsolver <- "cplexAPI"
    okcode   <- 1
  }
  # lpsolver <- "glpkAPI"; okcode <- 5 # DEBUG LINE -> force use of glpk
  if(verbose >= 1)
    cat(paste0("LP-solver: ", lpsolver, "\n"))

  pre_mod_list <- lapply(1:n.cores, function(x) { return(object@models) })
  fork_ids <- clusterApply(cl, pre_mod_list, fun = init_warm_mods,
                           lpsolver = lpsolver, okcode = okcode)
  rm(pre_mod_list)

  # keeping an eye on time
  t_start <- Sys.time()
  j <- 1

  # get grid field positions as data.table
  gridDT <- as.data.table(object@environ@field.pts)
  gridDT[, "field" := 1:.N]
  setkeyv(gridDT, c("z","x","y"))

  # get organims' scavenge radius
  get_scv_radius <- function(ctype) {
    unlist(lapply(ctype, function(x) object@models[[x]]@scavengeDist))
  }

  # track convergence of growth
  n_conv <- 5
  BM_converged <- F
  last_smass  <- rep(NA_real_, n_conv)

  #ram_usage <<- c(mem_used())
  # while loop that checks termination criteria before starting a new round
  while(j <= niter & difftime(Sys.time(), t_start, units = "mins") < lim_time & nrow(object@cellDT) < lim_cells & !BM_converged) {
    simRound <- object@n_rounds + 1

    # # give workers the current metabolite concentrations
    # fork_conc <- object@environ@concentrations
    # updateConc <- clusterApply(cl, 1:n.cores, function(fork_conc) {fieldConc <<- fork_conc; return(NULL)})
    # rm(updateConc)

    # get the cells' info
    #cellDT <- copy(object@cellDT)
    ncells <- nrow(object@cellDT)
    smass  <- sum(object@cellDT$mass)

    # check convergence
    last_smass <- c(last_smass[2:n_conv], smass)
    if(!any(is.na(last_smass))) {
      dBM <- abs(last_smass[n_conv] / min(last_smass[-n_conv]) - 1)
      #print(dBM)
      if(dBM < convergence.e)
        BM_converged <- TRUE
    }

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
      res[["x"]]      <- object@cellDT[i, get("x")]
      res[["y"]]      <- object@cellDT[i, get("y")]
      res[["size"]]   <- object@cellDT[i, get("size")]
      res[["scvr"]]   <- get_scv_radius(object@cellDT[i, get("type")]) # scavenge radius
      # TODO: There's probably a smarter way of doing the next command - but works for now
      # it gets all grid field in the cube surrounding the cell
      res[["gridLC"]] <- gridDT[abs(get("x") - res[["x"]]) <= (res[["size"]] /2 + res[["scvr"]]) &
                                  abs(get("y") - res[["y"]]) <= (res[["size"]] /2 + res[["scvr"]]) &
                                  abs(get("z") - 0         ) <= (res[["size"]] /2 + res[["scvr"]])]
      return(res)

    })

    cellFieldEnv <- parLapply(cl, pre_celFieldEnv, get_cellFieldEnv_ex)

    # close cells may claim in sum more than 100% of the resources from a field
    # Proportionally scale down accessibility in those cases:
    field_claim_sum = acc.prop = acc.prop.tmp = NULL
    cellFieldEnv <- rbindlist(cellFieldEnv, idcol = T)
    cellFieldEnv[, field_claim_sum := sum(acc.prop), by = "field.id"]
    cellFieldEnv[field_claim_sum > 1, acc.prop.tmp := acc.prop^exp(1)]
    cellFieldEnv[field_claim_sum > 1, acc.prop.tmp := acc.prop.tmp / sum(acc.prop.tmp) * field_claim_sum, by = "field.id"]
    cellFieldEnv[field_claim_sum > 1, acc.prop := acc.prop.tmp / field_claim_sum]

    cellFieldEnv <- lapply(1:ncells, FUN = function(x) {
      res <- list()
      res[["field.id"]]    <- cellFieldEnv[get(".id") == x, get("field.id")]
      res[["field.dist"]]  <- cellFieldEnv[get(".id") == x, get("field.dist")]
      res[["acc.prop"]]    <- cellFieldEnv[get(".id") == x, get("acc.prop")]

      res[["field.conc"]]  <- object@environ@concentrations[res[["field.id"]],]
      res[["field.cpds"]]  <- object@environ@compounds
      res[["field.execs"]] <- names(object@environ@exoenzymes)
      res[["fieldVol"]]    <- object@environ@fieldVol
      res[["model"]]       <- object@models[[object@cellDT[x, get("type")]]]
      res[["deltaTime"]]   <- object@deltaTime

      res[["type"]]        <- object@cellDT[x, get("type")]
      res[["cMass"]]       <- object@cellDT[x, get("mass")]
      res[["size"]]        <- object@cellDT[x, get("size")]

      # Chemotaxis
      res[["CT.x"]]       <- 0
      res[["CT.y"]]       <- 0
      ic_x <- object@cellDT[x, get("x")]
      ic_y <- object@cellDT[x, get("y")]
      for(icpd in object@models[[object@cellDT[x, get("type")]]]@chemotaxisCompound) {
        ind_ct <- which(object@models[[object@cellDT[x, get("type")]]]@chemotaxisCompound == icpd)
        ind <- which(res[["field.cpds"]] == icpd)
        delta_field_x <- object@environ@field.pts[res[["field.id"]]]$x - ic_x
        delta_field_y <- object@environ@field.pts[res[["field.id"]]]$y - ic_y
        # normalize
        delta_field_x <- delta_field_x / (res[["size"]]/2 + object@models[[object@cellDT[x, get("type")]]]@scavengeDist)
        delta_field_y <- delta_field_y / (res[["size"]]/2 + object@models[[object@cellDT[x, get("type")]]]@scavengeDist)

        f_conc <- res[["field.conc"]][,ind]
        f_conc_orig <- f_conc
        if(sum(f_conc) == 0)
          f_conc <- rep(1,length(f_conc))
        f_conc <- f_conc/sum(f_conc)

        # calculate sensing activity/strength
        centr_x    <- weighted.mean(delta_field_x, f_conc)
        centr_y    <- weighted.mean(delta_field_y, f_conc)
        Hill_ka    <- object@models[[object@cellDT[x, get("type")]]]@chemotaxisHillKA[ind_ct]
        Hill_coef  <- object@models[[object@cellDT[x, get("type")]]]@chemotaxisHillCoef[ind_ct]

        # calculate sensing strength based on hill equation and gradient steepness
        sensing_strength <- 1 / (1 + (Hill_ka / weighted.mean(f_conc_orig, f_conc))^Hill_coef)

        res[["CT.x"]] <- res[["CT.x"]] + centr_x * sensing_strength * object@models[[object@cellDT[x, get("type")]]]@chemotaxisStrength[ind_ct] * object@models[[object@cellDT[x, get("type")]]]@vmax * 60 * object@deltaTime * 60
        res[["CT.y"]] <- res[["CT.y"]] + centr_y * sensing_strength * object@models[[object@cellDT[x, get("type")]]]@chemotaxisStrength[ind_ct] * object@models[[object@cellDT[x, get("type")]]]@vmax * 60 * object@deltaTime * 60
      }

      return(res)
    })

    # - - - - - - - - - - - - #
    # (2) parallel agentFBA   #
    # - - - - - - - - - - - - #
    if(verbose > 1)
      cat("... performing cell-agent-FBA\n", sep ='')

    agFBA_results <- parLapply(cl, cellFieldEnv, agentFBA_ex)

    #print(table(unlist(lapply(agFBA_results, function(x) x$fba.stat))))

    # update the environment
    if(verbose > 1)
      cat("... update environment\n", sep ='')
    I <- rbindlist(lapply(agFBA_results, function(x) rbindlist(list(x$I.up, x$I.pd))))
    I <- I[, list("concChange" = sum(get("concChange"))), by = c("field.id", "concMatInd")]
    ind_var <- which(!object@environ@conc.isConstant)
    I <- I[get("concMatInd") %in% ind_var]
    I <- as.matrix(I)

    I[,3] <- ifelse(is.na(I[,3]), 0, I[,3]) # TODO: Check where NA originate from
    #print(I)

    object@environ@concentrations[I[,1:2]] <- object@environ@concentrations[I[,1:2]] + I[,3]
    object@environ@concentrations[I[,1:2]] <- ifelse(object@environ@concentrations[I[,1:2]] < 1e-12, 0, object@environ@concentrations[I[,1:2]])

    # extract reduced costs
    object@rcdt <- rbindlist(lapply(agFBA_results, function(x) x$rc), idcol = "cell")

    # exoenzymes
    I.exec <- rbindlist(lapply(agFBA_results, function(x) x$I.exec.pd))
    I.exec <- I.exec[, list("concChange" = sum(get("concChange"))), by = c("field.id", "concMatInd")]
    I.exec <- as.matrix(I.exec)

    I.exec[,3] <- ifelse(is.na(I.exec[,3]), 0, I.exec[,3]) # TODO: Check where NA originate come from
    object@environ@exoenzymes.conc[I.exec[,1:2]] <- object@environ@exoenzymes.conc[I.exec[,1:2]] + I.exec[,3]
    object@environ@exoenzymes.conc[I.exec[,1:2]] <- ifelse(object@environ@exoenzymes.conc[I.exec[,1:2]] < 1e-12, 0, object@environ@exoenzymes.conc[I.exec[,1:2]])
    #print(I.exec)

    # - - - - - - - - - - - - - #
    # (2.5) Exoenzyme activity  #
    # - - - - - - - - - - - - - #
    if(length(object@environ@exoenzymes) > 0) {
      if(verbose > 1)
        cat("... exoenzyme activity\n", sep ='')

      # Catalysis
      k <- 1
      for(exec in object@environ@exoenzymes) {
        exec.vmax <- exec@Kcat * object@environ@exoenzymes.conc[,k] / 1e6
        met_inds <- match(exec@mets, object@environ@compounds)

        S_0 <- object@environ@concentrations[,met_inds[1]] # current Substrate concentration

        # DOI: 10.1016/j.bej.2012.01.010
        S_t <- exec@Km * lambertW0(S_0/exec@Km * exp((S_0 - exec.vmax * object@deltaTime * 60 * 60)/exec@Km))

        dS <- S_0 - S_t

        # update concentrations
        #met_inds <- met_inds[!object@environ@conc.isConstant[met_inds]] # skip constant compounds
        for(i in 1:length(met_inds)) {
          # skip constant compounds
          if(object@environ@conc.isConstant[met_inds[i]] == F) {
            object@environ@concentrations[,met_inds[i]] <-
              object@environ@concentrations[,met_inds[i]] + exec@stoich[i] * dS
          }
        }

        k <- k + 1
      }

      # Exoenzyme decay
      all_lambda <- unlist(lapply(object@environ@exoenzymes, function(x) x@lambda))
      for(i in 1:length(all_lambda))
        object@environ@exoenzymes.conc[,i] <- object@environ@exoenzymes.conc[,i] * exp(-all_lambda[i] * object@deltaTime)
    }


    # - - - - - - - - - - - - #
    # (3) compound diffusion  #
    # - - - - - - - - - - - - #
    if(verbose > 1)
      cat("... diffusion of compounds\n", sep ='')
    #saveRDS(object@environ@concentrations, file = paste0("sim_conc_",object@n_rounds,".RDS"))
    object@environ <- diffuse_compounds(object@environ,
                                        deltaTime = object@deltaTime,
                                        cl = cl,
                                        n.cores = n.cores)

    # - - - - - - - - - #
    # (3.5) Chemotaxis  #
    # - - - - - - - - - #
    CT.x <- unlist(lapply(cellFieldEnv, function(x) x$CT.x))
    CT.y <- unlist(lapply(cellFieldEnv, function(x) x$CT.y))

    # add CT speed to existing speed
    object@cellDT$x.vel <- object@cellDT$x.vel + CT.x
    object@cellDT$y.vel <- object@cellDT$y.vel + CT.y

    # - - - - - - - - - - - - - -#
    # (4) grow and divide cells  #
    # - - - - - - - - - - - - - -#
    if(verbose > 1)
      cat("... binary fission of large cells\n", sep ='')

    # grow cells
    object@cellDT$mass <- unlist(lapply(agFBA_results, function(x) x$cMass_new))
    object@cellDT$size <- unlist(lapply(agFBA_results, function(x) x$cSize_new))


    object@cellDT[, "cellMassAtDivision" := cmad[get("type")]]

    gind <- which(object@cellDT$mass >= object@cellDT$cellMassAtDivision)
    if(length(gind) > 0) {
      newCells <- copy(object@cellDT[gind])
      newCells[, "parent" := get("cell")] # save parent information
      newCells[, "cell" := 1:.N]
      newCells[, "cell" := get("cell") + ncells]

      # update cell mass and size of divided cells
      object@cellDT <- rbind(object@cellDT, newCells)
      object@cellDT[get("mass") >= get("cellMassAtDivision"), "size" := get("size") * 0.5^(1/3)]
      object@cellDT[get("mass") >= get("cellMassAtDivision"), "mass" := get("mass") / 2]
      rm(newCells)

      # place daughter cell next to parent cell (random angle/direction)
      object@cellDT[get("cell") > ncells, "x" := get("x") + get("size") * cos(stats::runif(.N, max = 2*pi))]
      object@cellDT[get("cell") > ncells, "y" := get("y") + get("size") * sin(stats::runif(.N, max = 2*pi))]

      # invert velocity of daughter cells
      object@cellDT[get("cell") > ncells, "x.vel" := get("x.vel") * (-1)]
      object@cellDT[get("cell") > ncells, "y.vel" := get("y.vel") * (-1)]
    }

    object@cellDT[, "cellMassAtDivision" := NULL]

    ncells_new <- nrow(object@cellDT)

    cvmax <- unlist(lapply(object@models, function(x) x@vmax)) * 3600 * object@deltaTime # max speed in µm per deltaTime

    # Brownian motion of cells (assuming normal distribution)
    r_motion <- matrix(stats::rnorm(ncells_new * 2,sd = object@rMotion * 60 * object@deltaTime),
                       ncol = 2)

    # set up particle simulation
    sim_new <- tbl_graph(nodes = object@cellDT, directed = F, node_key = "cell") %>%
      simulate(setup = predefined_genesis(x = object@cellDT$x,
                                          y = object@cellDT$y,
                                          x_vel = object@cellDT$x.vel + r_motion[,1],
                                          y_vel = object@cellDT$y.vel + r_motion[,2])) %>%
      wield(collision_force, radius = object@cellDT$size/2, n_iter = 50, strength = 0.7) %>%
      impose(velocity_constraint,
             vmax = cvmax[object@cellDT$type]) %>%
      impose(polygon_constraint,
             polygon = object@universePolygon)

    #saveRDS(sim_new, file ="sim_new.RDS")

    # - - - - - - - - - - - - - - - - - - - - - - - - #
    # (5) Move & collide / reorganize cells in space  #
    # - - - - - - - - - - - - - - - - - - - - - - - - #
    if(verbose > 1)
      cat("... sliding and colliding cells\n", sep ='')
    sim_new <- evolve(sim_new, steps = round(object@deltaTime * 60))

    object@cellDT$x <- sim_new$position[,1]
    object@cellDT$y <- sim_new$position[,2]
    object@cellDT$x.vel <- sim_new$velocity[,1]
    object@cellDT$y.vel <- sim_new$velocity[,2]

    # Live Plot
    if(live.plot)
      print(plot_cells(object))

    # - - - - - - - - #
    # (6) Record data #
    # - - - - - - - - #
    if(verbose > 1)
      cat("... recording simulation data\n", sep ='')

    object@history[[simRound+1]] <- list()

    # cell positions
    object@history[[simRound+1]]$cells <- copy(object@cellDT)

    # global metabolite concentrations (only variable)
    ind_var <- !object@environ@conc.isConstant
    dt_conc_tmp <- data.table(cpd.id = object@environ@compounds[ind_var],
                              cpd.name = object@environ@compound.names[ind_var],
                              global_concentration = apply(object@environ@concentrations[,ind_var],2,mean))
    object@history[[simRound+1]]$global_compounds <- dt_conc_tmp


    # Cell-individual exchange fluxes in fmol
    cell.ex.fluxes <- lapply(agFBA_results, function(icell) {
      data.table(compound = names(icell$ex.fluxes),
                 exflux   = icell$ex.fluxes)
    })
    cell.ex.fluxes <- rbindlist(cell.ex.fluxes, idcol = "cell")
    object@history[[simRound+1]]$cell.exchanges <- cell.ex.fluxes

    # specific compound concentrations per grids
    if(!is.null(record) && (any(grepl("^compound_", record)) | ("compounds" %in% record))) {
      COI <- record[grepl("^compound_", record)] # metabolites of interest
      COI <- gsub("^compound_","",COI)
      COI <- COI[COI %in% object@environ@compounds]
      if("compounds" %in% record)
        COI <- object@environ@compounds

      if(length(COI) > 0) {
        rec_file <- paste0(object@recordDir,"/cpdrec_",simRound+1,"_0.RDS")
        k_tmp <- 1
        while(file.exists(rec_file)) {
          rec_file <- paste0(object@recordDir,"/cpdrec_",simRound+1,"_",k_tmp,".RDS")
          k_tmp <- k_tmp + 1
        }

        COI.ind <- match(COI, object@environ@compounds)
        COI.rec <- object@environ@concentrations[,COI.ind, drop = FALSE]
        COI.rec <- data.table(round(COI.rec, digits = 7))

        # save compound recording in recording file
        saveRDS(COI.rec, file = rec_file, compress = T)

        # save order of compounds
        object@history[[simRound+1]]$compounds <- object@environ@compounds[COI.ind]

        # save path to recording file
        object@history[[simRound+1]]$compounds.record <- rec_file
      }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - #
    # (7) Apply user function to simulation object  #
    # - - - - - - - - - - - - - - - - - - - - - - - #
    if (!is.null(on.iteration)) {
      tmpobj <- on.iteration(object,
                             ...
      )
      if (class(tmpobj) == "growthSimulation") object <- tmpobj
    }

    # - - - - - - - - - #
    # Finish iteration  #
    # - - - - - - - - - #

    # small cell summary
    if(verbose > 0) {
      cellsum <- copy(object@cellDT[, list("mass" = round(sum(get("mass")), digits = 2)), by = "type"])
      cellsum[, "tmp.sum" := paste0(get("type"),"(",get("mass"),")")]
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

#
#
#
#
#' @import sybil
#'
init_warm_mods <- function(x, lpsolver, okcode) {
  require(Eutropia)

  SYBIL_SETTINGS("SOLVER",lpsolver)
  ok <<- okcode

  fork_mods <- list()
  for(mi in names(x)) {
    fork_mods[[mi]] <- sysBiolAlg(x[[mi]]@mod,
                                  algorithm = "mtf2",
                                  pFBAcoeff = 1e-6)
  }
  fork_mods <<- fork_mods

  return(NULL)
}

#' @import sp
#' @import data.table
get_cellFieldEnv_ex <- function(x) {
  icell_x    <- x$x
  icell_y    <- x$y
  icell_size <- x$size
  scv_dist   <- x$scvr

  gsp  <- SpatialPoints(x$gridLC[,list("x"=get("x"),"y"=get("y"),"z"=get("z"))])

  csp <- SpatialPoints(matrix(c(icell_x, icell_y, 0), ncol = 3))
  q.c.dist <- spDists(gsp, csp)[,1]

  qry.env <- which(q.c.dist <= (icell_size/2 + scv_dist))

  #qry.dist <- gDistance(gsp[qry.env,], csp, byid = T)[1,]
  qry.dist <- q.c.dist[qry.env]

  # get accessible portion of environment grids based on distance
  accPortion <-  - 1 / scv_dist * (qry.dist - (icell_size/2 + scv_dist))
  accPortion <- ifelse(accPortion > 1, 1, accPortion)

  return(data.table(field.id   = x$gridLC[qry.env, get("field")],
                    field.dist = qry.dist,
                    acc.prop   = accPortion))
}



#' @import sybil
#' @import data.table
agentFBA_ex <- function(x) {

  # (2.0) Get local accessible compounds
  accCompounds <- scavenge_compounds(localEnv     = x[c("acc.prop","field.id","field.dist")],
                                     env_conc     = x$field.conc,
                                     env_cpds     = x$field.cpds,
                                     env_fieldVol = x$fieldVol)
  updatedEX    <- adjust_uptake(model      = x$model@mod,
                                cMass      = x$cMass,
                                accCpdFMOL = accCompounds,
                                deltaTime  = x$deltaTime)

  # (2.1) agentFBA(model, envGrids, curMass) for independent cells
  # constrain model to local environment
  ccbnds <- changeColsBnds(problem(fork_mods[[x$type]]), updatedEX$ex.react.ind,
                           lb = updatedEX$ex.react.lb, ub = x$model@mod@uppbnd[updatedEX$ex.react.ind])

  sol.fba <- optimizeProb(fork_mods[[x$type]])
  mu <- sol.fba$obj * x$deltaTime
  stat <- sol.fba$stat

  # extract reduced costs of exchange reactions
  rc <- getRedCosts(fork_mods[[x$type]]@problem)
  rc <- data.table(type      = x$type,
                   exrxn     = x$model@mod@react_id[updatedEX$ex.react.ind],
                   exname    = x$model@mod@react_name[updatedEX$ex.react.ind],
                   flux      = sol.fba$fluxes[updatedEX$ex.react.ind],
                   red.costs = rc[updatedEX$ex.react.ind],
                   lb.org    = x$model@mod@lowbnd[updatedEX$ex.react.ind],
                   lb.env    = updatedEX$ex.react.lb)

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
  if(sol.fba$stat != ok)
    flx <- rep(0, x$model@mod@react_num)

  # get exchange reactions with non-zero flux
  # adjusting them to time and cell mass
  ex.ind <- grep("^EX_",x$model@mod@react_id)
  ex.flx <- sol.fba$fluxes[ex.ind]
  names(ex.flx) <- x$model@mod@react_id[ex.ind]
  ex.flx <- ex.flx[abs(ex.flx) > 0]
  # normalize to time and cell mass
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

  #------------------------#
  # Production of Exozymes #
  #------------------------#
  if(length(x$model@exoenzymes) > 0) {
    exec.prod <- cMass_new * x$deltaTime * x$model@exoenzymes.prod / 1000 # Total Exoenzyme production in absolute fmol in this time step and by this cell of its specific mass
    names(exec.prod) <- x$model@exoenzymes

    field.frac <- max(accCompounds$field.dist) - accCompounds$field.dist
    field.frac <- field.frac/sum(field.frac)

    exec.prod <- (exec.prod / 1e12) / (x$fieldVol / 1e15) * 1e3 # nM produced if all given to a single field

    exec.pro.mat <- field.frac %*% t(exec.prod)

    exec.concMatInd <- match(colnames(exec.pro.mat), x$field.execs)

    I.exec.pd <- data.table(field.id   = rep(accmat.row2field, ncol(exec.pro.mat)),
                            concMatInd = rep(exec.concMatInd, each =nrow(exec.pro.mat)),
                            concChange = as.numeric(exec.pro.mat))
  } else {
    I.exec.pd <- data.table(field.id   = integer(0),
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
              I.pd      = I.pd,
              I.exec.pd = I.exec.pd,
              rc        = rc
  ))
}

#' @import sybil
close_clusters <- function(x){
  for(mi in names(fork_mods)) {
    if(SYBIL_SETTINGS("SOLVER") == "cplexAPI") {
      cplexAPI::delProbCPLEX(fork_mods[[mi]]@problem@oobj@env, fork_mods[[mi]]@problem@oobj@lp)
      cplexAPI::closeEnvCPLEX(fork_mods[[mi]]@problem@oobj@env)
    } else {
      glpkAPI::delProbGLPK(fork_mods[[mi]]@problem@oobj)
    }
  }
  return(NULL)
}
