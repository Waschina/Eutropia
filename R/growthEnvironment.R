
#' Structure of the S4 class "growthEnvironment"
#'
#' @aliases growthEnvironment
#'
#' @exportClass growthEnvironment
#'
#' @import Matrix
#' @import sf
#'
#' @slot field.pts Object of class 'XYZ'. Coordinates of rhombic dodecahedron
#' centroids.
#' @slot compounds Character vector with compound IDs
#' @slot compound.names Character vector with compound names
#' @slot compound.D Numeric vector of the diffusion coefficient values for each
#' compound. Unit: \eqn{\mu}m^2/s
#' @slot concentrations (n x m) numeric matrix with n columns representing compounds
#' and m grid field (field.pts). Units of matrix entries: mM
#' @slot conc.isConstant Logical vector indicating if the respective compound is
#' fixed (i.e. buffered) in its concentration.
#' @slot fieldLayers Integer specifying how many layers (z-dimension) of rhomic
#' dodecahedra fields are in the environment representation.
#' @slot nfields Integer indicating the number of fields in the environment
#' representation.
#' @slot fieldSize Size of each rhomic dodecahedron in \eqn{\mu}m as the distance between
#' opposite faces.
#' @slot fieldVol Volume of each field. Unit: \eqn{\mu}m^3
#' @slot mat.in A two column numeric matrix as a representation of the grid environment
#' as a directed graph. First column (From) and second column (To) are denoting
#' field's indices in `field.pts`
#' @slot mat.out A two column matrix. First column: indices of each field in `field.pts`;
#' second column: Number of neighboring fields
#' @slot exoenzymes A list with the \link{Exoenzyme} S4 objects.
#' @slot exoenzymes.conc Same as `concentrations`, but for the exoenzymes and in
#' nM.
#'
setClass("growthEnvironment",

         slots = c(
           field.pts       = "XYZ",
           compounds       = "character",
           compound.names  = "character",
           compound.D      = "numeric",
           concentrations  = "matrix",
           conc.isConstant = "logical",
           fieldLayers     = "integer",
           nfields         = "integer",
           fieldSize       = "numeric",
           fieldVol        = "numeric",
           mat.in          = "matrix",
           mat.out         = "matrix",

           # Exoenzymes

           exoenzymes      = "list",
           exoenzymes.conc = "matrix" # Concentrations of exoenzymes per field [mg/L]
         )
)

# Initialising a cell growth environment
#
# @param polygon.coords 2-column numeric matrix with coordinates; first point (row) should equal last coordinates (row)
# @param field.size Is the diameter of the circumscribed sphere. Its relation to a is: z = 4*a/sqrt(3)
# @param expand Length of field grid exceeding the \code{polygon.cords}. Should be higher than 2x\code{field.size}. Really?
#
#' @import sf
#' @import data.table
setMethod("initialize", "growthEnvironment",
          function(.Object,
                   polygon.coords, field.size, field.layers = 3, expand = field.size, ...) {
            .Object <- callNextMethod(.Object, ...)

            field.layers <- as.integer(field.layers)
            # Important: Distance of field centers of neighboring fields: 2*H = 1*R_i = 2 * a * sqrt(2/3)
            # equal double the distance from center to plane faces

            # Make base layer of fields (Rhombic dodecahedrons)
            area <- st_polygon(list(polygon.coords))
            #plot(area, col = "blue")

            # expand the polygon boundaries a bit (bleeding)
            fieldGrid <- st_buffer(area, expand)
            #plot(fieldGrid, add = TRUE, alpha = 0.5)

            # create raster grid
            field.pts_base <- st_make_grid(fieldGrid, cellsize = field.size,
                                           what = "centers")
            field.pts_base <- field.pts_base[fieldGrid]
            #plot(field.pts_base, add = TRUE)


            field.pts_base <- lapply(field.pts_base, as.matrix)
            field.pts_base <- do.call(rbind, field.pts_base)
            field.pts_base <- cbind(field.pts_base,
                                    matrix(0, ncol = 1,
                                           nrow = nrow(field.pts_base)))
            colnames(field.pts_base) <- c("x","y","z")

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
            field.pts <- st_multipoint(field.pts)

            nfields <- nrow(field.pts)

            .Object@field.pts       <- field.pts
            .Object@compounds       <- character(0)
            .Object@compound.names  <- character(0)
            .Object@compound.D      <- double(0)
            .Object@concentrations  <- matrix(ncol = 0, nrow = nfields)
            .Object@conc.isConstant <- logical(0)
            .Object@fieldLayers     <- field.layers
            .Object@nfields         <- nfields
            .Object@fieldSize       <- field.size
            .Object@fieldVol        <- fieldVol

            # Exoenzymes
            .Object@exoenzymes      <- list()
            .Object@exoenzymes.conc <- matrix(ncol = 0, nrow = nfields)

            # Diffusion kernel
            DCM <- build.DCM(field.pts, field.size) # Diffusion coefficient matrix

            dt <- data.table(which(DCM != 0, arr.ind = T))
            setnames(dt, new = c("from","to"))
            setkey(dt, "to")
            dt[, n_neighbors := .N, by = from]

            mat.in  <- as.matrix(dt[from != to, .(from, to)])
            mat.out <- as.matrix(dt[from == to, .(from, n_neighbors)])

            .Object@mat.in  <- mat.in
            .Object@mat.out <- mat.out

            return(.Object)
          }
)


setGeneric(name="build.DCM",
           def=function(field.pts, field.size)
           {
             standardGeneric("build.DCM")
           }
)

# Build diffusion coefficient matirx
setMethod(f          = "build.DCM",
          signature  = signature(field.pts      = "XYZ",
                                 field.size     = "numeric"),
          definition = function(field.pts, field.size) {

            fs <- field.size

            lis <- as.data.table(field.pts[,])
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
            DCM[as.matrix(lis[,1:2])] <- 12/13 * 1/12
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
#' @import parallel
diffuse_compounds <- function(object, deltaTime, cl, n.cores) {

  # - - - - - #
  # Compounds #
  # - - - - - #
  n.chunks <- n.cores
  ind_variable <- which(apply(object@concentrations,2,sd) > 0 & object@compound.D > 0)



  # no variable compounds -> no need for diffusion
  if(length(ind_variable) > 0) {

    if(length(ind_variable) / n.chunks < 3)
      n.chunks <- ceiling(n.chunks/2)

    # split the metabolites, which are subject to diffusion in (nearly) equally
    # sized chunks, so chunks can be computed in parallel.
    if(n.chunks > 1) {
      ind_var_chunks <- split(ind_variable,
                              cut(seq_along(ind_variable),
                                  n.chunks, labels = F))
    } else {
      ind_var_chunks <- list(`1` = ind_variable)
    }

    # VIA RccpArmadillo
    conc.list.tmp <- lapply(ind_var_chunks, function(x) {

      diff_shere_surface_area <- object@compound.D[x] * 60 * 60 * deltaTime # surface area in micro-m^2 after one iteration step
      diff_shere_radius       <- sqrt(diff_shere_surface_area / (4*pi))
      diffusion.niter         <- ceiling(diff_shere_radius/object@fieldSize)

      list(mat.in  = object@mat.in,
           mat.out = object@mat.out,
           conc    = as.matrix(object@concentrations[,x]),
           n_iter  = diffusion.niter)
    })
    #print(ind_var_chunks[[1]])
    conc.list.tmp <- parLapply(cl, conc.list.tmp, indDiff_worker)

    for(k in 1:length(ind_var_chunks)) {
      for(i in 1:ncol(conc.list.tmp[[k]])) {
        object@concentrations[,ind_var_chunks[[k]][i]] <- conc.list.tmp[[k]][,i]
      }
    }
  }

  # - - - - - - #
  # Exoenzymes  #
  # - - - - - - #
  if(length(object@exoenzymes) > 0) {
    exec.D <- unlist(lapply(object@exoenzymes, function(x) x@D))
    ind_variable <- which(apply(object@exoenzymes.conc,2,sd) > 0 & exec.D > 0)

    if(length(ind_variable) > 0) {
      diff_shere_surface_area <- exec.D[ind_variable] * 60 * 60 * deltaTime # surface area in micro-m^2 after one iteration step
      diff_shere_radius       <- sqrt(diff_shere_surface_area / (4*pi))
      diffusion.niter         <- ceiling(diff_shere_radius/object@fieldSize)

      exec.conc.tmp <- list(mat.in  = object@mat.in,
                            mat.out = object@mat.out,
                            conc    = as.matrix(object@exoenzymes.conc[,ind_variable, drop = FALSE]),
                            n_iter  = diffusion.niter)

      exec.conc.tmp <- indDiff_worker(exec.conc.tmp)

      for(i in 1:length(ind_variable)) {
        object@exoenzymes.conc[,ind_variable[i]] <- exec.conc.tmp[,i]
      }

    }
  }

  return(object)
}

indDiff_worker <- function(x) {
  #res <- diffChangeVec(x$mat.in, x$mat.out, x$conc, x$n_iter)
  res <- diffChange(x$mat.in, x$mat.out, x$conc, x$n_iter)
  return(res)
}

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Compound diffusion   #
# with openMP in c++11 #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
diffuse_compoundsPar <- function(object, deltaTime, cl, n.cores) {

  # - - - - - #
  # Compounds #
  # - - - - - #
  ind_variable <- which(apply(object@concentrations,2,sd) > 0 & object@compound.D > 0 & object@conc.isConstant == FALSE)
  ind_infD     <- ind_variable[is.infinite(object@compound.D[ind_variable])]
  ind_variable <- ind_variable[!(ind_variable %in% ind_infD)]

  # no variable compounds -> no need for diffusion
  if(length(ind_variable) > 0) {

    # calculate number of steps required per compound
    diff_shere_surface_area <- object@compound.D[ind_variable] * 60 * 60 * deltaTime # surface area in micro-m^2 after one iteration step
    diff_shere_radius       <- sqrt(diff_shere_surface_area / (4*pi))
    diffusion.niter         <- ceiling(diff_shere_radius/object@fieldSize)

    # VIA RccpArmadillo + OpenMP
    object@concentrations <- diffChangePar(object@mat.in,
                                           object@mat.out,
                                           as.matrix(object@concentrations),
                                           diffusion.niter,
                                           ind_variable,
                                           n.cores)

  }
  if(length(ind_infD) > 0) {
    if(length(ind_infD) == 1) {
      object@concentrations[,ind_infD] <- mean(object@concentrations[,ind_infD])
    } else {
      object@concentrations[,ind_infD] <- apply(object@concentrations[,ind_infD],
                                                2,
                                                function(x) rep(mean(x), length(x)))
    }
  }

  # - - - - - - #
  # Exoenzymes  #
  # - - - - - - #
  if(length(object@exoenzymes) > 0) {
    exec.D <- unlist(lapply(object@exoenzymes, function(x) x@D))
    ind_variable <- which(apply(object@exoenzymes.conc,2,sd) > 0 & exec.D > 0)
    ind_infD     <- ind_variable[is.infinite(exec.D[ind_variable])]
    ind_variable <- ind_variable[!(ind_variable %in% ind_infD)]

    if(length(ind_variable) > 0) {
      diff_shere_surface_area <- exec.D[ind_variable] * 60 * 60 * deltaTime # surface area in micro-m^2 after one iteration step
      diff_shere_radius       <- sqrt(diff_shere_surface_area / (4*pi))
      diffusion.niter         <- ceiling(diff_shere_radius/object@fieldSize)

      object@exoenzymes.conc <- diffChangePar(object@mat.in, object@mat.out,
                                              as.matrix(object@exoenzymes.conc),
                                              diffusion.niter,
                                              ind_variable,
                                              min(n.cores,
                                                  length(ind_variable)))

    }
    if(length(ind_infD) > 0) {
      if(length(ind_infD) == 1) {
        object@exoenzymes.conc[,ind_infD] <- mean(object@exoenzymes.conc[,ind_infD])
      } else {
        object@exoenzymes.conc[,ind_infD] <- apply(object@exoenzymes.conc[,ind_infD],
                                                   2,
                                                   function(x) rep(mean(x), length(x)))
      }
    }
  }

  return(object)
}
