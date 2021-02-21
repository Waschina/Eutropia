setClass("growthEnvironment",

         slots = c(
           hex.pts         = "SpatialPoints",
           compounds      = "character",
           concentrations = "matrix",
           hexHeight      = "numeric",
           nhex           = "integer",
           hexSize        = "numeric",
           hexVol         = "numeric",
           DCM            = "Matrix"
         )
)

#' Initialising a cell growth environment
#'
#' @param polygon.coords 2-column numeric matrix with coordinates; first point (row) should equal last coordinates (row)
#' @param hexagon.size Hexagon grid cell size (edge to opoosite edge). in µm
#' @param hexHeight Height of hexagonal prisms. in µm
#' @param expand Length of hexagon grid exceeding the \code{polygon.cords}. Should be higher than 2x\code{hexagon.size}.
#'
#' @example
#' # a sqaure as boundary polygon:
#' build.growthEnvironment(polygon.coords = cbind(c(-100, -100, 100, 100), c(-100, 100, 100, -100)),
#'                         hexagon.size = 0.5)
setMethod("initialize", "growthEnvironment",
          function(.Object,
                   polygon.coords, hexagon.size, diffusion.alpha, hex.height = 10, expand = 1, ...) {
            .Object <- callNextMethod(.Object, ...)

            area <- Polygon(polygon.coords)
            area <- Polygons(list(area), "s1")

            area.sp <- SpatialPolygons(list(area))

            hexGrid <- gBuffer(area.sp, width = expand)

            hex.pts <- spsample(hexGrid,
                                type="hexagonal",
                                cellsize=hexagon.size)

            nhex <- nrow(hex.pts@coords)

            .Object@hex.pts   <- hex.pts
            .Object@compounds <- character(0)
            .Object@concentrations <- matrix(ncol = 0, nrow = nhex)
            .Object@hexHeight <- hex.height
            .Object@nhex      <- nhex
            .Object@hexSize   <- hexagon.size

            # calculate a hexagons volume (in µm^3 or Litre*1e-15)
            hex.a <- tan(pi/3) * hexagon.size
            .Object@hexVol    <- (3 * sqrt(3) / 2) * hex.a^2 * hex.height

            .Object@DCM <- build.DCM(.Object, alpha = diffusion.alpha) # Diffusion coefficient matrix

            return(.Object)
          }
)



# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# show method for small summary #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
setMethod(f          = "show",
          signature  = signature(object = "growthEnvironment"),
          definition = function(object) {
            cat("Cell growth environment\n")
            cat("  Number of hexagons:\t\t",object@nhex,"\n")
            cat("  Number of variable compounds:\t",length(object@compounds),"\n")
            cat("  Hexgonal prism height:\t", object@hexHeight,"\n")
          }
)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Adding compounds     #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="add.compounds",
           def=function(object, compounds, concentrations)
           {
             standardGeneric("add.compounds")
           }
)


#' Add compounds to cell growth environment.
#'
#' @param object Either a 'growthEnvironment' or \code{growthSimulation} object.
#' @param compound Character vector of compound id.
#' @param concentrations Numeric vector of same length as \code{compound} with the compounds concentration in mM
#' @return \code{growthEnvironment} object with updated compound concentrations.
setMethod(f          = "add.compounds",
          signature  = signature(object         = "growthEnvironment",
                                 compounds      = "character",
                                 concentrations = "numeric"),
          definition = function(object, compounds, concentrations) {

            if(length(compounds) != length(concentrations))
              stop("Lengths of 'compounds' and 'concentrations' differ.")

            for(i in 1:length(compounds)) {
              if(compounds[i] %in% object@compounds) {
                c_ind <- which(object@compounds == compounds[i])
                object@concentrations[, c_ind] <- object@concentrations[, c_ind] + concentrations[i]
              } else {
                object@compounds      <- c(object@compounds, compounds[i])
                object@concentrations <- cbind(object@concentrations,
                                               matrix(concentrations[i], ncol = 1, nrow = object@nhex))
              }
            }

            return(object)
          }
)

setMethod(f          = "add.compounds",
          signature  = signature(object         = "growthSimulation",
                                 compounds      = "character",
                                 concentrations = "numeric"),
          definition = function(object, compounds, concentrations) {

            object@environ <- add.compounds(object@environ, compounds, concentrations)

            return(object)
          }
)

setGeneric(name="build.DCM",
           def=function(object, alpha)
           {
             standardGeneric("build.DCM")
           }
)

# Build diffusion coefficient matirx
setMethod(f          = "build.DCM",
          signature  = signature(object         = "growthEnvironment",
                                 alpha          = "numeric"),
          definition = function(object, alpha) {

            hs <- object@hexSize

            lis <- as.data.table(object@hex.pts)
            lis[, id := 1:.N]

            # identify neigbors
            lis_L  <- copy(lis[,.(x = x + cos(pi*6/3)*hs, y = y + sin(pi*6/3)*hs, id.L = 1:.N)])
            lis_LO <- copy(lis[,.(x = x + cos(pi*5/3)*hs, y = y + sin(pi*5/3)*hs, id.LO = 1:.N)])
            lis_RO <- copy(lis[,.(x = x + cos(pi*4/3)*hs, y = y + sin(pi*4/3)*hs, id.RO = 1:.N)])
            lis_R  <- copy(lis[,.(x = x + cos(pi*3/3)*hs, y = y + sin(pi*3/3)*hs, id.R = 1:.N)])
            lis_RU <- copy(lis[,.(x = x + cos(pi*2/3)*hs, y = y + sin(pi*2/3)*hs, id.RU = 1:.N)])
            lis_LU <- copy(lis[,.(x = x + cos(pi*1/3)*hs, y = y + sin(pi*1/3)*hs, id.LU = 1:.N)])

            lis[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]

            lis_L[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]
            lis_LO[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]
            lis_RO[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]
            lis_R[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]
            lis_RU[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]
            lis_LU[, `:=`(x = round(x, digits = 3), y = round(y, digits = 3))]

            lis <- merge(lis, lis_L, by = c("x","y"), all.x = T)
            lis <- merge(lis, lis_LO, by = c("x","y"), all.x = T)
            lis <- merge(lis, lis_RO, by = c("x","y"), all.x = T)
            lis <- merge(lis, lis_R, by = c("x","y"), all.x = T)
            lis <- merge(lis, lis_RU, by = c("x","y"), all.x = T)
            lis <- merge(lis, lis_LU, by = c("x","y"), all.x = T)

            lis <- lis[order(id)]

            # build DCM
            lis <- melt(lis, id.vars = "id", measure.vars = c("id.L","id.LO","id.RO","id.R","id.RU","id.LU"))
            lis <- lis[!is.na(value), .(id, value)][order(id)]

            DCM <- Matrix(0, ncol = max(lis$id), nrow = max(lis$id), sparse = T)
            DCM[as.matrix(lis[,1:2])] <- alpha * 1/6
            diag(DCM) <- 1 - Matrix::rowSums(DCM)
            DCM <- Matrix::t(DCM)

            return(DCM) # Diffusion coefficient matrix
          }
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Compound diffusion   #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
setGeneric(name="diffuse.compounds",
           def=function(object, n_iter,  ...)
           {
             standardGeneric("diffuse.compounds")
           }
)


setMethod(f          = "diffuse.compounds",
          signature  = signature(object         = "growthEnvironment",
                                 n_iter         = "numeric"),
          definition = function(object, n_iter) {

            conctmp <- object@concentrations
            for(i in 1:n_iter)
              conctmp <- object@DCM %*% conctmp

            object@concentrations <- as.matrix(conctmp)

            # low_ind <- which(object@concentrations > 0 & object@concentrations < 1e-8, arr.ind = T) # TODO specify this in an argument / global setting
            #
            # if(nrow > 0)
            #   object@concentrations[low_ind] <- 0

            return(object)
          }
)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# Plot compound distribution  #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
setGeneric(name="plot.environment",
           def=function(object, compounds, compound.names = NULL, xlim = NULL, ylim = NULL, ...)
           {
             standardGeneric("plot.environment")
           }
)

#' Plot spatial distribution of compounds
#'
#' @param object An object of class \code{growthEnvrionment} or \code{growthSimulation}
#' @param compounds Character string of compound name to plot
#' @param compound.names Character string of compound names that should be desplayed as facet header instead of compound IDs from \code{compounds}.
#' @param xlim Numeric vector of length 2 specifying the limits (left and right) for x axis; i.e. the horizontal dimension.
#' @param ylim Numeric vector of length 2 specifying the limits (top and bottom) for y axis; i.e. the vertical dimension.
#' @param iter Positive integer number of the simulation step/iteration to plot the distribution. Works only if \code{object} is of
#' class \code{growthSimulation} and if the respective metabolite concentrations were recorded (see \code{link{run.simulation}}).
setMethod(f = "plot.environment",
          signature = signature(object = "growthEnvironment",
                                compounds = "character"),
          definition = function(object, compounds, compound.names = NULL, xlim = NULL, ylim = NULL, ...) {

            # Argument sanity check
            if(!all(compounds %in% object@compounds)) {
              compounds <- compounds[compounds %in% object@compounds]

              if(length(compounds) == 0)
                stop("Selected compound(s) not in list of variable environment compounds.")

              warning("Not all selected compounds are in the list of variable environment compounds. Continuing with the rest...")
            }

            # If no compound names are defined, use compound IDs instead
            cpd_nameDT <- data.table(Compound = compounds, Compound.name = compounds)
            if(!is.null(compound.names)) {
              if(length(compound.names) != length(compounds)){
                warning("Number of compound names not equal to number of compound IDs. Using compound IDs instead.")
              } else if(any(duplicated(compound.names))) {
                warning("Duplicate compound names. Using compound IDs instead.")
              } else {
                cpd_nameDT$Compound.name <- compound.names
              }
            }

            envDT <- data.table(object@concentrations)
            names(envDT) <- object@compounds
            envDT <- envDT[, ..compounds]

            envDT <- cbind(object@hex.pts@coords, envDT)
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
              geom_hex(aes(colour = mM), stat = "identity") +
              coord_equal(xlim = xlim, ylim = ylim) +
              geom_segment(aes(x = xlim[2]-bar_wd-x_exp_fac, xend = xlim[2]-x_exp_fac,
                               y = ylim[1]+y_exp_fac, yend = ylim[1]+y_exp_fac), color = "white") +
              annotate("text", x = xlim[2]-bar_wd/2-x_exp_fac, y = ylim[1]+y_exp_fac, label = paste0(bar_wd," µm"),
                       color = "white", hjust = 0.5, vjust = -1, size = 2.5) +
              theme_bw() +
              scale_fill_viridis_c() + scale_color_viridis_c() + facet_wrap("Compound.name") +
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

setMethod(f = "plot.environment",
          signature = signature(object = "growthSimulation",
                                compounds = "character"),
          definition = function(object, compounds, compound.names = NULL, xlim = NULL, ylim = NULL, iter = NULL) {

            # If no limits are defined use polygon universe limits
            if(is.null(xlim)) {
              xlim <- c(min(object@universePolygon[,1]),
                        max(object@universePolygon[,1]))
            }
            if(is.null(ylim)) {
              ylim <- c(min(object@universePolygon[,2]),
                        max(object@universePolygon[,2]))
            }

            if(is.null(iter))
              return(plot.environment(object@environ, compounds, compound.names = compound.names, xlim = xlim, ylim = ylim))

            # in case a specific simulation round is selected
            # TODO

          }
)
