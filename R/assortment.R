#' @title Calculate spatial segregation and assortment
#'
#' @description This functions calculates cell-wise and mean cell type
#' segregation and assortment values as suggested by Yanni et al. (2019) Current
#' Biology.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param r Numeric. Maximum distance to to other cells to considered as
#' neighbour. Refers to the surface-to-surface distance of cells. Unit: \eqn{\mu}m.
#' @param iter Positive integer number of the simulation step/iteration at which
#' the cell assortment should be calculated. If NULL, current distribution is
#' displayed.
#' @param n_simulations positive integer, specifying the number of permutations
#' to be used to estimate the expected assortment if cell types are
#' distributed randomly.
#'
#' @return A list with three elements: 'cells' is a data table with the
#' cell-wise calculated segregation and assortment values. 'summary_observed' is
#' a summary (mean and sd) of the observed segregation and assortment by cell
#' types. 'summary_expected' is a summary of expected segregation and assortment
#' if cell types are reassigned randomly to cells in the environment while
#' maintaining the cell positions and community-wide cell type counts.
#'
#' @details
#' Cell segregation and assortment is calculated as suggested by Yanni et al.
#' (2019) Current Biology.
#'
#' @export
calc_assortment <- function(object, r = 5, iter = NULL, n_simulations = 10) {
  # sanity checks
  if(object@n_rounds == 0 & !is.null(iter)) {
    if(iter != 0)
      warning("Simulation did not run yet. Showing initial simulation status.")
    iter <- NULL
  }

  if(nrow(object@cellDT) == 0)
    stop("No cells in the environment")

  #--------------------#
  # get cell positions #
  #--------------------#
  cellposDT <- data.table()
  # current
  if(is.null(iter)) {
    cellposDT <- copy(object@cellDT)
  }

  # from history
  if(!is.null(iter)) {
    if(iter > object@n_rounds) {
      warning(paste0("Simulation did not run ",iter," iterations yet. Displaying results for last iteration (",object@n_rounds,")"))
      iter <- object@n_rounds
    }
    i_round <- iter
    cellposDT <- copy(object@history[[iter+1]]$cells)
  }

  # core function for assortment and segregation calculations
  tmp_assortment <- function(dttmp) {
    c_assortment <- rep(NA_real_, nrow(dttmp))
    c_segregation <- rep(NA_real_, nrow(dttmp))

    for(i in 1:nrow(dttmp)) {
      c_x <- dttmp[i,x]
      c_y <- dttmp[i,y]
      c_size <- dttmp[i,size]
      c_type <- dttmp[i,type]
      freq_type <- dttmp[type == c_type, .N]/dttmp[,.N]
      dttmp[, tmp_dist := sqrt((x-c_x)^2 + (y-c_y)^2)]
      dist_dt_tmp <- copy(dttmp[-1])[tmp_dist <= c_size/2 + size/2 + r]
      if(nrow(dist_dt_tmp) > 0) {
        dist_dt_tmp[, tmp_score := -1]
        dist_dt_tmp[type == c_type, tmp_score := 1]
        c_segr_tmp <- 1/((pi * (r + c_size/2)^2) - (pi * (c_size/2)^2)) * sum(dist_dt_tmp$tmp_score)
      } else {
        c_segr_tmp <- 0
      }
      c_segregation[i] <- c_segr_tmp
      c_assortment[i] <- (c_segr_tmp-freq_type) / (1-freq_type)
    }

    return(data.table(cell = 1:nrow(dttmp),
                      segregation = c_segregation,
                      assortment = c_assortment))
  }

  asttmp <- tmp_assortment(cellposDT)

  out_dt <- copy(cellposDT[,.(cell, type, x, y)])
  out_dt <- merge(out_dt, asttmp, by = "cell")

  summary_observed <-  out_dt[, .(mean_segregation = mean(segregation),
                                  sd_segregation = sd(segregation),
                                  mean_assortment = mean(assortment),
                                  sd_assortment = sd(assortment)), by = "type"]


  # get expected assortment
  simulated_assortment <- list()

  for(j in 1:n_simulations) {
    cat("\r",j)
    cellposDT$type <- sample(cellposDT$type)
    dt_tmp <- copy(cellposDT[,.(cell, type, x, y)])
    simulated_assortment[[j]] <- merge(dt_tmp,
                                       tmp_assortment(cellposDT), by = "cell")
  }
  cat("\n")
  simulated_assortment <- rbindlist(simulated_assortment, idcol = "simrun")

  summary_expected <-  simulated_assortment[, .(mean_segregation = mean(segregation),
                                                sd_segregation = sd(segregation),
                                                mean_assortment = mean(assortment),
                                                sd_assortment = sd(assortment)),
                                            by = "type"]

  return(list(cells = out_dt,
              summary_observed = summary_observed,
              summary_expected = summary_expected))
}
