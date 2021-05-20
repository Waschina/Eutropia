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

# Initialising a cell growth environment
#
# @param polygon.coords 2-column numeric matrix with coordinates; first point (row) should equal last coordinates (row)
# @param field.size Is the diameter of the circumscribed sphere. Its relation to a is: z = 4*a/sqrt(3)
# @param expand Length of field grid exceeding the \code{polygon.cords}. Should be higher than 2x\code{field.size}. Really?
#
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

            # # diffusion for n iterations
            # for(i in 1:120)
            #   DCM <- DCM %*% Matrix::t(DCM)

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


            # with RccpArmadillo
            #conctmp <- Matrix(object@concentrations, sparse = T)

            #conctmp <- diffuse_arma(object@DCM, conctmp, n_iter)
            #conctmp <- diffuse_eigen(object@DCM, conctmp, n_iter)

            object@concentrations <- as.matrix(conctmp)

            return(object)
          }
)
