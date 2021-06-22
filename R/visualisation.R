# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot cell distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Plot spatial distribution of cells
#'
#' @description Plots the distribution of cells in a growth simulation.
#'
#' @param object S4-object of type \code{growthSimulation}.
#' @param xlim Numeric vector of length 2, specifying the x-range to be displayed.
#' @param ylim Numeric vector of length 2, specifying the y-range to be displayed.
#' @param iter Positive integer number of the simulation step/iteration to plot the cell distribution. If past simulation are displayed, the cell postitions needed to be recorded when running the simulation before (see \code{link{run.simulation}}). If NULL, current distribution is displayed.
#' @param scalebar.color Color of the scale bar and its annotation. Default: "white".
#'
#' @return A ggplot object.
#'
#' @import ggforce
#' @import ggplot2
#' @import hms
#'
#' @export
plot_cells <- function(object, xlim = NULL, ylim = NULL, iter = NULL,
                       scalebar.color = "white") {

  # sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  # If no limits are defined use polygon universe limits
  if(is.null(xlim)) {
    xlim <- c(min(object@universePolygon[,1]),
              max(object@universePolygon[,1]))
  }
  if(is.null(ylim)) {
    ylim <- c(min(object@universePolygon[,2]),
              max(object@universePolygon[,2]))
  }

  i_round <- object@n_rounds

  cellposDT <- data.table()
  # current
  if(is.null(iter)) {
    cellposDT <- object@cellDT
  }

  # from history
  if(!is.null(iter)) {
    if(iter > object@n_rounds) {
      warning(paste0("Simulation did not run ",iter," iterations yet. Displaying results for last iteration (",object@n_rounds,")"))
      iter <- object@n_rounds
    }
    i_round <- iter
    cellposDT <- object@history[[iter]]$cells
  }

  # get expansion factor so scale bar is not to close to panel margins
  x_exp_fac <- (xlim[2]-xlim[1]) * 0.05
  y_exp_fac <- (ylim[2]-ylim[1]) * 0.05

  # get scale bar length to span approx 10% of x-axis
  x_span <- xlim[2]-xlim[1]
  x_magn <- floor(log(x_span, base = 10))
  bar_wd <- 10^x_magn / 10
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd * 2.5, bar_wd) # e.g. 10 -> 25
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd / 2.5 * 5, bar_wd) # e.g. 10 -> 50

  # get coordinated to draw growth polygon fences
  kidt <- data.table(object@universePolygon)
  colnames(kidt) <- c("x","y")

  n_types <- length(unique(cellposDT$type))
  col_code <- "Set1"
  if(n_types > 8)
    col_code <- "Set3"

  p <- ggplot(cellposDT, aes(x0 = x,y0 = y)) +
    geom_polygon(data = kidt, mapping = aes(x = x, y = y),
                 col = "#333333", fill = "black", linetype = "solid") +
    geom_circle(aes(r = size/2, fill = type, col = type), alpha = 0.3,
                show.legend = T, key_glyph = draw_key_point) +
    coord_equal(xlim = xlim, ylim = ylim) +
    geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac,
                     xend = xlim[2]-x_exp_fac,
                     y = ylim[1]+y_exp_fac,
                     yend = ylim[1]+y_exp_fac), color = scalebar.color) +
    annotate("text",
             x = xlim[2]-bar_wd/2-x_exp_fac,
             y = ylim[1]+y_exp_fac,
             label = paste0(bar_wd," \u00B5m"),
             color = scalebar.color, hjust = 0.5, vjust = -1, size = 2.5)

  # use nicer brewer colors if here are not to many distinct cell types
  # Otherwise: using ggplot defaults
  if(n_types <= 12) {
    p <- p + scale_color_brewer(palette = col_code)
  }

  p <- p +
    scale_y_continuous(sec.axis = sec_axis(~ .)) + scale_x_continuous(sec.axis = sec_axis(~ .)) +
    theme_bw() +
    theme(plot.background = element_blank(),
          axis.line.x.top = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.x.bottom = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.left = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.right = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white"),
          legend.background = element_blank(),
          legend.text = element_text(color = "black", face = "italic"),
          legend.box.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "bottom", legend.direction="vertical"
    ) +
    guides(color = guide_legend(title = "Organism",
                                override.aes = list(alpha = 1, size = 5)),
           fill  = guide_legend(title = "Organism")) +
    labs(subtitle = round_hms(hms(hours = object@deltaTime * i_round),
                              digits = 0))

  return(p)

}


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot compound distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Plot spatial distribution of compounds
#'
#' @description Plots the spatial distribution of compound concentrations in a
#' heatmap.
#'
#' @param object S4-object of type \link{growthSimulation}.
#' @param compounds Character string of compound ids to plot
#' @param compound.names Character string of compound names that should be
#' displayed as facet header instead of compound names from the simulation object.
#' @param xlim Numeric vector of length 2 specifying the limits (left and right)
#' for x axis; i.e. the horizontal dimension.
#' @param ylim Numeric vector of length 2 specifying the limits (top and bottom)
#' for y axis; i.e. the vertical dimension.
#' @param iter Positive integer number of the simulation step/iteration to plot
#' the distribution. Works only if the respective metabolite
#' concentrations were prior recorded (see \link{run_simulation}). If `NULL`,
#' current distribution of compound concentrations.
#' @param scalebar.color Color of the scale bar and its annotation. Default:
#' "white".
#' @param layer Positive integer, specifying which layer (z-dimension) to plot.
#' Default: 0 (base plane).
#' @param gradient.limits Numeric vector of length 2 specifying the
#' concentration range that is spanned by the color gradient.
#' @param gradient.option Character string indicating the colormap to be used
#' for visualizing metabolite concentrations. Please refer to \link[ggplot2]{scale_colour_viridis_d}
#' so see possible options.
#'
#' @return A \link{ggplot}
#'
#' @import ggplot2
#' @import hms
#'
#' @export
plot_environment <- function(object, compounds, compound.names = NULL,
                             xlim = NULL, ylim = NULL, iter = NULL,
                             scalebar.color = "white", layer = 0,
                             gradient.limits = NULL, gradient.option = "viridis") {

  # sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")

  i_round <- object@n_rounds

  # Argument sanity check
  if(!all(compounds %in% object@environ@compounds)) {
    compounds <- compounds[compounds %in% object@environ@compounds]

    if(length(compounds) == 0)
      stop("Selected compound(s) not in list of variable environment compounds.")

    warning("Not all selected compounds are in the list of variable environment compounds. Continuing with the rest...")
  }

  # If no limits are defined use polygon universe limits
  if(is.null(xlim)) {
    xlim <- c(min(object@universePolygon[,1]),
              max(object@universePolygon[,1]))
  }
  if(is.null(ylim)) {
    ylim <- c(min(object@universePolygon[,2]),
              max(object@universePolygon[,2]))
  }

  # If no compound names are defined, use original compound names instead
  cpd_nameDT <- data.table(Compound = compounds,
                           Compound.name = object@environ@compound.names[match(compounds, object@environ@compounds)])
  if(!is.null(compound.names)) {
    if(length(compound.names) != length(compounds)){
      warning("Number of compound names not equal to number of compound IDs. Using original compound names instead.")
    } else if(any(duplicated(compound.names))) {
      warning("Duplicate compound names. Using original compound names instead.")
    } else {
      cpd_nameDT$Compound.name <- compound.names
    }
  }

  # which layer to plot? By default: Baseplane
  zvals <- sort(unique(object@environ@field.pts@coords[,3]))
  zvalOI <- zvals[1]
  if(layer < 0 | layer > (length(zvals)-1) | layer %% 1 != 0) {
    warning("Invalid layer. Plotting data for base plane (0).")
  } else {
    zvalOI <- zvals[layer + 1]
  }

  if(!(length(iter) %in% c(0,1))) {
    iter <- iter[1]
  }

  if(!is.null(iter) && iter > object@n_rounds) {
    warning(paste0("Simulation did not yet run ",iter, " iterations.",
                   " Displaying lastest iteration status. (",
                   object@n_rounds,")"))
    iter <- NULL
  }

  # which time step to plot?
  if(is.null(iter) || iter == object@n_rounds) {
    envDT <- data.table(object@environ@concentrations)
    names(envDT) <- object@environ@compounds
  } else {
    if(is.null(object@history[[iter]]$compounds))
      stop(paste0("No compound concentrations were recorded at iteration ",
                  iter, "."))

    compounds <- compounds[compounds %in% object@history[[iter]]$compounds]
    if(length(compounds) == 0)
      stop(paste0("None of the compounds have recorded concentrations at iteration",
                 iter, "."))

    if(!file.exists(object@history[[iter]]$compounds.record))
      stop(paste0("Compound recording file ",
                  object@history[[iter]]$compounds.record,
                  " not found."))

    i_round <- iter

    envDT <- readRDS(object@history[[iter]]$compounds.record)
    names(envDT) <- object@history[[iter]]$compounds
  }


  envDT <- envDT[, ..compounds]

  envDT <- cbind(object@environ@field.pts@coords, envDT)
  envDT <- envDT[z == zvalOI] # baseplane
  envDT <- melt(envDT, id.vars = c("x","y"), variable.name = "Compound", value.name = "mM")
  envDT <- merge(envDT, cpd_nameDT)

  # get expansion factor so scale bar is not to close to panel margins
  x_exp_fac <- (xlim[2]-xlim[1]) * 0.05
  y_exp_fac <- (ylim[2]-ylim[1]) * 0.05

  # get scale bar length to span approx 10% of x-axis
  x_span <- xlim[2]-xlim[1]
  x_magn <- floor(log(x_span, base = 10))
  bar_wd <- 10^x_magn / 10
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd * 2.5, bar_wd) # e.g. 10 -> 25
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd / 2.5 * 5, bar_wd) # e.g. 10 -> 50

  p <- ggplot(envDT, aes(x,y, fill = mM)) +
    #geom_hex(aes(colour = mM), stat = "identity") +
    geom_bin2d(aes(colour = mM), stat = "identity") +
    coord_equal(xlim = xlim, ylim = ylim) +
    geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac, xend = xlim[2]-x_exp_fac,
                     y = ylim[1]+y_exp_fac, yend = ylim[1]+y_exp_fac), color = scalebar.color) +
    annotate("text",
             x = xlim[2]-bar_wd/2-x_exp_fac,
             y = ylim[1]+y_exp_fac,
             label = paste0(bar_wd," \u00B5m"),
             color = scalebar.color, hjust = 0.5, vjust = -1, size = 2.5) +
    theme_bw() +
    scale_fill_viridis_c(guide = guide_colourbar(direction = "vertical"), limits = gradient.limits,
                         option = gradient.option) +
    scale_color_viridis_c(limits = gradient.limits, option = gradient.option) +
    facet_wrap("Compound.name") +
    scale_y_continuous(sec.axis = sec_axis(~ .)) + scale_x_continuous(sec.axis = sec_axis(~ .)) +
    theme(legend.position = "right",
          axis.line.x.top = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.x.bottom = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.left = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.right = element_line(color = "white", size = 1.5, lineend = "round"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black")) +
    labs(caption = round_hms(hms(hours = object@deltaTime * i_round),
                             digits = 0))


  return(p)

}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot exoenzyme distribution #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Plot spatial distribution of exoenzymes
#'
#' @description Plots the spatial distribution of exoenzyme concentrations in a
#' heatmap.
#'
#' @param object S4-object of type \link{growthSimulation}
#' @param exoenzymes Character string of exoenzyme ids to plot
#' @param exoenzyme.names Character string of exoenzyme names that should be
#' displayed as facet header instead of exoenzyme names from the simulation
#' object.
#' @param xlim Numeric vector of length 2 specifying the limits (left and right)
#' for x axis; i.e. the horizontal dimension.
#' @param ylim Numeric vector of length 2 specifying the limits (top and bottom)
#' for y axis; i.e. the vertical dimension.
#' @param iter Positive integer number of the simulation step/iteration to plot
#' the distribution. Works only if the respective metabolite
#' concentrations were prior recorded (see \link{run_simulation}). If
#' \code{NULL}, current distribution of exoenzyme concentrations.
#' @param scalebar.color Color of the scale bar and its annotation. Default:
#' "white".
#' @param layer Positive integer, specifying which layer (z-dimension) to plot.
#' Default: 0 (base plane).
#' @param gradient.limits Numeric vector of length 2 specifying the
#' concentration range that is spanned by the color gradient.
#' @param gradient.option Character string indicating the colormap to be used
#' for visualizing exoenzyme concentrations. Please refer to \link[ggplot2]{scale_colour_viridis_d}
#' so see possible options.
#'
#' @return A \link{ggplot}
#'
#' @import ggplot2
#'
#' @export
plot_environment_exoenzymes <- function(object, exoenzymes, exoenzyme.names = NULL,
                                        xlim = NULL, ylim = NULL, iter = NULL,
                                        scalebar.color = "white", layer = 0,
                                        gradient.limits = NULL, gradient.option = "viridis") {

  # Argument sanity check
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")
  if(!all(exoenzymes %in% names(object@environ@exoenzymes))) {
    exoenzymes <- exoenzymes[exoenzymes %in% names(object@environ@exoenzymes)]

    if(length(exoenzymes) == 0)
      stop("Selected exoenzyme(s) not part of the environment/simulation.")

    warning("Not all selected exoenzymes are part of the environment/simulation. Continuing with the ones that are...")
  }

  # If no limits are defined use polygon universe limits
  if(is.null(xlim)) {
    xlim <- c(min(object@universePolygon[,1]),
              max(object@universePolygon[,1]))
  }
  if(is.null(ylim)) {
    ylim <- c(min(object@universePolygon[,2]),
              max(object@universePolygon[,2]))
  }

  # If no exoenzyme names are defined, use original exoenzyme names instead
  exec_names <- unlist(lapply(object@environ@exoenzymes,
                              function(x) x@name))
  exec_nameDT <- data.table(Exoenzyme = exoenzymes,
                            Exoenzyme.name = exec_names[match(exoenzymes,
                                                              names(object@environ@exoenzymes))])
  if(!is.null(exoenzyme.names)) {
    if(length(exoenzyme.names) != length(exoenzymes)){
      warning("Number of exoenzyme names not equal to number of exoenzyme IDs. Using original exoenzyme names instead.")
    } else if(any(duplicated(exoenzyme.names))) {
      warning("Duplicate exoenzyme names. Using original exoenzyme names instead.")
    } else {
      exec_nameDT$Exoenzyme.name <- exoenzyme.names
    }
  }

  # which layer to plot? By default: Baseplane
  zvals <- sort(unique(object@environ@field.pts@coords[,3]))
  zvalOI <- zvals[1]
  if(layer < 0 | layer > (length(zvals)-1) | layer %% 1 != 0) {
    warning("Invalid layer. Plotting data for base plane (0).")
  } else {
    zvalOI <- zvals[layer + 1]
  }

  # which time step to plot?
  if(is.null(iter)) {
    envDT <- data.table(object@environ@exoenzymes.conc)
  } else {
    # TODO: When a former simulation step is selected
  }

  names(envDT) <- names(object@environ@exoenzymes)
  envDT <- envDT[, ..exoenzymes]

  envDT <- cbind(object@environ@field.pts@coords, envDT)
  envDT <- envDT[z == zvalOI] # baseplane
  envDT <- melt(envDT, id.vars = c("x","y"), variable.name = "Exoenzyme",
                value.name = "nM")
  envDT <- merge(envDT, exec_nameDT)

  # get expansion factor so scale bar is not to close to panel margins
  x_exp_fac <- (xlim[2]-xlim[1]) * 0.05
  y_exp_fac <- (ylim[2]-ylim[1]) * 0.05

  # get scale bar length to span approx 10% of x-axis
  x_span <- xlim[2]-xlim[1]
  x_magn <- floor(log(x_span, base = 10))
  bar_wd <- 10^x_magn / 10
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd * 2.5, bar_wd) # e.g. 10 -> 25
  bar_wd <- ifelse(bar_wd/x_span < 0.05, bar_wd <- bar_wd / 2.5 * 5, bar_wd) # e.g. 10 -> 50

  p <- ggplot(envDT, aes(x,y, fill = nM)) +
    #geom_hex(aes(colour = mM), stat = "identity") +
    geom_bin2d(aes(colour = nM), stat = "identity") +
    coord_equal(xlim = xlim, ylim = ylim) +
    geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac,
                     xend = xlim[2]-x_exp_fac,
                     y = ylim[1]+y_exp_fac,
                     yend = ylim[1]+y_exp_fac),
                 color = scalebar.color) +
    annotate("text", x = xlim[2]-bar_wd/2-x_exp_fac,
             y = ylim[1]+y_exp_fac,
             label = paste0(bar_wd," \u00B5m"),
             color = scalebar.color, hjust = 0.5, vjust = -1, size = 2.5) +
    theme_bw() +
    scale_fill_viridis_c(guide = guide_colourbar(direction = "vertical"), limits = gradient.limits,
                         option = gradient.option) +
    scale_color_viridis_c(limits = gradient.limits, option = gradient.option) +
    facet_wrap("Exoenzyme.name") +
    scale_y_continuous(sec.axis = sec_axis(~ .)) + scale_x_continuous(sec.axis = sec_axis(~ .)) +
    theme(legend.position = "right",
          axis.line.x.top = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.x.bottom = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.left = element_line(color = "white", size = 1.5, lineend = "round"),
          axis.line.y.right = element_line(color = "white", size = 1.5, lineend = "round"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          strip.background = element_blank(), strip.text = element_text(face = "bold", color = "black"))


  return(p)

}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot growth curves      #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Plot growth curves
#'
#' @description Plot population growth over time
#'
#' @param object S4-object of class \link{growthSimulation}.
#' @param tlim Numeric vector of length 2, specifying the x-range (Time) to be displayed.
#' @param ylim Numeric vector of length 2, specifying the y-range (Mass) to be displayed.
#'
#' @return A \link{ggplot}
#'
#' @import ggplot2
#'
#' @export
plot_growth <- function(object, tlim = NULL, ylim = NULL) {

  # Sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")
  if(object@n_rounds < 2)
    stop("Simulation did not yet run for at least 2 iterations. No growth can be plotted.")
  if(is.null(object@history[[1]]$cells))
    stop("No cell/population data has been recorded.")

  dt_cellsg <- rbindlist(lapply(object@history, function(x) x$cells), idcol = "iteration")
  dt_cellsg <- dt_cellsg[, .(mass = sum(mass)), by = .(iteration, type)]
  dt_cellsg[, time := iteration * object@deltaTime]

  p_growth <- ggplot(dt_cellsg, aes(time, mass, col = type, group = type)) +
    geom_line() +
    #geom_point(cex = 0.5) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Time [hr]", y = "Mass [pg]",col = "Organism") +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(xlim = tlim, ylim = ylim) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          legend.text = element_text(color = "black", face = "italic"))

  return(p_growth)
}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot compound time curves #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#' @title Plot compound time curves
#'
#' @description Plot concentrations of metabolites over time
#'
#' @param object S4-object of class \link{growthSimulation}.
#' @param compounds Vector of compound IDs whose concentration is plotted. If no IDs are provided,
#' the top variable (SD) compounds are plotted (max. 8 compounds).
#' @param tlim Numeric vector of length 2, specifying the x-range (Time) to be displayed.
#' @param ylim Numeric vector of length 2, specifying the y-range (mM) to be displayed.
#'
#' @return A \link{ggplot}.
#'
#' @import ggplot2
#'
#' @export
plot_compounds <- function(object, compounds = NULL, tlim = NULL, ylim = NULL) {

  # Sanity checks
  if(!is.growthSimulation(object))
    stop("'Object' not of class 'growthSimulation'.")
  if(object@n_rounds < 2)
    stop("Simulation did not yet run for at least 2 iterations. No growth can be plotted.")
  if(is.null(object@history[[1]]$global_compounds))
    stop("No compound concentration data has been recorded.")

  dt_cpdg <- rbindlist(lapply(object@history, function(x) x$global_compounds),
                       idcol = "iteration")
  dt_cpdg[, time := iteration * object@deltaTime]


  if(is.null(compounds)) {
    # identify most variable compounds
    rel_cpds <- dt_cpdg[, sd(global_concentration), by = cpd.id][V1 > 1e-4][order(-V1), cpd.id]
    rel_cpds <- rel_cpds[1:min(c(8, length(rel_cpds)))]
  } else {
    # check if slected compounds are part of models/simulation
    compounds <- compounds[compounds %in% object@environ@compounds]
    if(length(compounds) == 0)
      stop("Selected compounds not part of the simulation.")
    rel_cpds <- compounds
  }

  # filter to relevant cpds
  dt_cpdg <- dt_cpdg[cpd.id %in% rel_cpds]

  p_cpdth <- ggplot(dt_cpdg, aes(time, global_concentration,
                                 col = cpd.name, group = cpd.name)) +
    geom_line() +
    #geom_point(cex = 0.5) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Time [hr]", y = "Concentration [mM]",col = "Compound") +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(xlim = tlim, ylim = ylim) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          legend.text = element_text(color = "black"))

  return(p_cpdth)

}
