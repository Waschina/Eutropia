setClass("growthEnvironment",

         slots = c(
           field.pts      = "SpatialPoints",
           compounds      = "character",
           concentrations = "matrix",
           fieldLayers    = "integer",
           nfields        = "integer",
           fieldSize      = "numeric",
           fieldVol       = "numeric",
           DCM            = "Matrix"
         )
)

#' Initialising a cell growth environment
#'
#' @param polygon.coords 2-column numeric matrix with coordinates; first point (row) should equal last coordinates (row)
#' @param field.size Is the diameter of the circumscribed sphere. Its relation to a is: z = 4*a/sqrt(3)
#' @param expand Length of field grid exceeding the \code{polygon.cords}. Should be higher than 2x\code{field.size}. Really?
#'
#' @example TODO
setMethod("initialize", "growthEnvironment",
          function(.Object,
                   polygon.coords, field.size, diffusion.alpha, field.layers = 3, expand = field.size, ...) {
            .Object <- callNextMethod(.Object, ...)

            field.layers <- as.integer(field.layers)
            # Important: Distance of field centers of neighboring fields: 2*H = 1*R_i = 2 * a * sqrt(2/3)
            # equal double the distance from center to plane faces

            # Make base layer of field (Rhombic dodecahedrons)
            area <- Polygon(polygon.coords)
            area <- Polygons(list(area), "s1")

            area.sp <- SpatialPolygons(list(area))

            fieldGrid <- gBuffer(area.sp, width = expand)

            field.pts_base <- spsample(fieldGrid,
                                       type="regular",
                                       cellsize=field.size)
            field.pts_base <- as.matrix(field.pts_base@coords)
            field.pts_base <- cbind(field.pts_base, matrix(0, ncol = 1, nrow = nrow(field.pts_base)))
            colnames(field.pts_base) <- c("x","y","z")

            #field.pts_base <- SpatialPoints(field.pts_base)

            # Make additional layers
            #fieldHeight    <- sqrt(6)/3 * field.size  # z-axis difference between field centers of neighboring layers
            field.a        <- field.size * 3 / (2 * sqrt(6)) # edge length (main metric of the Rhombic dodecahedrons)
            fieldVol       <- 16/9 * field.a^3 * sqrt(3) # volume calculation for Rhombic dodecahedrons

            field.pts <- field.pts_base
            if(field.layers > 1) {
              i <- 2
              while(i <= field.layers) {
                lay_i <- field.pts_base
                lay_i[,"z"] <- field.size * (i - 1)
                if(i %% 2 == 0) {
                  lay_i[,"x"] <- lay_i[,"x"] - field.size/2
                  lay_i[,"y"] <- lay_i[,"y"] - field.size/2
                }
                field.pts <- rbind(field.pts,lay_i)
                i <- i + 1
              }
            }
            field.pts <- SpatialPoints(field.pts)

            nfields <- nrow(field.pts@coords)

            .Object@field.pts      <- field.pts
            .Object@compounds      <- character(0)
            .Object@concentrations <- matrix(ncol = 0, nrow = nfields)
            .Object@fieldLayers    <- field.layers
            .Object@nfields        <- nfields
            .Object@fieldSize      <- field.size
            .Object@fieldVol       <- fieldVol

            .Object@DCM <- build.DCM(field.pts, field.size,  alpha = diffusion.alpha) # Diffusion coefficient matrix

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
            cat("    Universe dimensions (µm):\t\t",
                abs(round(min(object@field.pts@coords[,1])-max(object@field.pts@coords[,1]), digits = 2))," x ",
                abs(round(min(object@field.pts@coords[,2])-max(object@field.pts@coords[,2]), digits = 2))," x ",
                abs(round(min(object@field.pts@coords[,3])-max(object@field.pts@coords[,3]), digits = 2)),"\n")
            cat("    Universe volume (µm^3):\t\t", round(object@fieldSize * object@nfields, digits = 2) , "\n")
            cat("    Number of rhombic dodecahedrons:\t",object@nfields,"\n")
            cat("    Number of variable compounds:\t",length(object@compounds),"\n")

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
                                               matrix(concentrations[i], ncol = 1, nrow = object@nfields))
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
           def=function(field.pts, field.size, alpha)
           {
             standardGeneric("build.DCM")
           }
)

# Build diffusion coefficient matirx
setMethod(f          = "build.DCM",
          signature  = signature(field.pts      = "SpatialPoints",
                                 field.size     = "numeric",
                                 alpha          = "numeric"),
          definition = function(field.pts, field.size, alpha) {

            fs <- field.size

            lis <- as.data.table(field.pts)
            lis[, id := 1:.N]

            # correct for rounding issues
            all_x <- sort(unique(lis$x)); all_x <- data.table(s = all_x, val = all_x)
            all_y <- sort(unique(lis$y)); all_y <- data.table(s = all_y, val = all_y)
            all_z <- sort(unique(lis$z)); all_z <- data.table(s = all_z, val = all_z)
            setattr(all_x, "sorted", "s"); setkey(all_x, s)
            setattr(all_y, "sorted", "s"); setkey(all_y, s)
            setattr(all_z, "sorted", "s"); setkey(all_z, s)

            # M : mid layer, U upper layer, L - lower layer
            # l - links, o - oben, u - unten, r - rechts
            lis_C_o  <- copy(lis[,.(x = x        , y = y + fs  , z = z     , id.Co  = 1:.N)])[order(x,y,z)]
            lis_C_r  <- copy(lis[,.(x = x + fs   , y = y       , z = z     , id.Cr  = 1:.N)])[order(x,y,z)]
            lis_C_l  <- copy(lis[,.(x = x - fs   , y = y       , z = z     , id.Cl  = 1:.N)])[order(x,y,z)]
            lis_C_u  <- copy(lis[,.(x = x        , y = y - fs  , z = z     , id.Cu  = 1:.N)])[order(x,y,z)]

            lis_U_ol  <- copy(lis[,.(x = x - fs/2, y = y + fs/2, z = z + fs, id.Uol  = 1:.N)])[order(x,y,z)]
            lis_U_or  <- copy(lis[,.(x = x + fs/2, y = y + fs/2, z = z + fs, id.Uor  = 1:.N)])[order(x,y,z)]
            lis_U_ur  <- copy(lis[,.(x = x + fs/2, y = y - fs/2, z = z + fs, id.Uur  = 1:.N)])[order(x,y,z)]
            lis_U_ul  <- copy(lis[,.(x = x - fs/2, y = y - fs/2, z = z + fs, id.Uul  = 1:.N)])[order(x,y,z)]

            lis_L_ol  <- copy(lis[,.(x = x - fs/2, y = y + fs/2, z = z - fs, id.Lol  = 1:.N)])[order(x,y,z)]
            lis_L_or  <- copy(lis[,.(x = x + fs/2, y = y + fs/2, z = z - fs, id.Lor  = 1:.N)])[order(x,y,z)]
            lis_L_ur  <- copy(lis[,.(x = x + fs/2, y = y - fs/2, z = z - fs, id.Lur  = 1:.N)])[order(x,y,z)]
            lis_L_ul  <- copy(lis[,.(x = x - fs/2, y = y - fs/2, z = z - fs, id.Lul  = 1:.N)])[order(x,y,z)]

            correct_numbers <- function(dt) {
              # x
              qn <- dt$x
              new_n <- all_x[J(qn), roll = -Inf]$val
              new_n <- ifelse(abs(new_n-qn) > field.size/10, NA, new_n)
              dt$x <- new_n
              # y
              qn <- dt$y
              new_n <- all_y[J(qn), roll = -Inf]$val
              new_n <- ifelse(abs(new_n-qn) > field.size/10, NA, new_n)
              dt$y <- new_n
              # z
              qn <- dt$z
              new_n <- all_z[J(qn), roll = -Inf]$val
              new_n <- ifelse(abs(new_n-qn) > field.size/10, NA, new_n)
              dt$z <- new_n

              return(dt)
            }

            lis_C_o <- correct_numbers(lis_C_o)
            lis_C_r <- correct_numbers(lis_C_r)
            lis_C_l <- correct_numbers(lis_C_l)
            lis_C_u <- correct_numbers(lis_C_u)

            lis_U_ol <- correct_numbers(lis_U_ol)
            lis_U_or <- correct_numbers(lis_U_or)
            lis_U_ur <- correct_numbers(lis_U_ur)
            lis_U_ul <- correct_numbers(lis_U_ul)

            lis_L_ol <- correct_numbers(lis_L_ol)
            lis_L_or <- correct_numbers(lis_L_or)
            lis_L_ur <- correct_numbers(lis_L_ur)
            lis_L_ul <- correct_numbers(lis_L_ul)

            #lis <- correct_numbers(lis)

            #setkeyv(lis, keycols)
            #setkeyv(lis_C_lo, keycols)

            lis <- merge(lis, lis_C_o, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_C_r, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_C_l, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_C_u, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_U_ol, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_U_or, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_U_ur, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_U_ul, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_L_ol, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_L_or, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_L_ur, by = c("x","y","z"), all.x = T)
            lis <- merge(lis, lis_L_ul, by = c("x","y","z"), all.x = T)

            lis <- lis[order(id)]

            table(apply(lis[,5:16], 1, function(x) sum(!is.na(x))))

            # build DCM
            lis <- melt(lis, id.vars = "id", measure.vars = c("id.Co","id.Cr","id.Cu","id.Cl",
                                                              "id.Uol","id.Uor","id.Uur","id.Uul",
                                                              "id.Lol","id.Lor","id.Lur","id.Lul"))
            lis <- lis[!is.na(value), .(id, value)][order(id)]

            DCM <- Matrix(0, ncol = max(lis$id), nrow = max(lis$id), sparse = T)
            DCM[as.matrix(lis[,1:2])] <- alpha * 1/12
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

            envDT <- cbind(object@field.pts@coords, envDT)
            envDT <- envDT[z == 0]
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
              theme(legend.position = "bottom",
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
