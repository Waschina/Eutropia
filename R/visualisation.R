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
#' concentrations were prior recorded (see \link{run.simulation}). If
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
#' @return A ggplot object.
#'
#' @export
setGeneric(name="plot.environment.exoenzymes",
           def=function(object, exoenzymes, exoenzyme.names = NULL,
                        xlim = NULL, ylim = NULL, iter = NULL,
                        scalebar.color = "white", layer = 0,
                        gradient.limits = NULL, gradient.option = "viridis", ...)
           {
             standardGeneric("plot.environment.exoenzymes")
           }
)

setMethod(f = "plot.environment.exoenzymes",
          signature = signature(object = "growthSimulation",
                                exoenzymes = "character"),
          definition = function(object, exoenzymes, exoenzyme.names = NULL,
                                xlim = NULL, ylim = NULL, iter = NULL,
                                scalebar.color = "white", layer = 0,
                                gradient.limits = NULL, gradient.option = "viridis") {

            # Argument sanity check
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
              geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac, xend = xlim[2]-x_exp_fac,
                               y = ylim[1]+y_exp_fac, yend = ylim[1]+y_exp_fac), color = scalebar.color) +
              annotate("text", x = xlim[2]-bar_wd/2-x_exp_fac, y = ylim[1]+y_exp_fac, label = paste0(bar_wd," µm"),
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
)
