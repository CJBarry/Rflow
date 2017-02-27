# Rflow package - plotting


# plot matrix images ------------------------------------------------------

#' Plot MODFLOW variable matrix
#'
#' @param mtx
#' matrix;
#' variable to plot
#' @param gccs
#' numeric \code{[NCOL + 1]};
#' column divider \eqn{x} co-ordinates, in ascending order
#' @param grcs
#' numeric \code{[NROW + 1]};
#' row divider \eqn{y} co-ordinates, in ascending order
#' @param zlim
#' numeric \code{[2]} or character string;
#' min. and max. values to plot; alternatively use \code{breaks} argument
#'  as in \code{graphics::image}; additionally, three character string
#'  options are given: \code{"auto"} for automatic range calculation
#'  (pretty), \code{"sym"} for automatic range that is symmetric about \eqn{0}
#'  and "angle" for a range between \eqn{-pi} and \eqn{pi}
#' @param col
#' character \code{[]};
#' colours
#' @param show.range
#' logical \code{[1]};
#' print computed \code{zlim} to screen? (in case that
#'  \code{zlim = "angle"}, shows which directions the colours represent)
#' @param ...
#' additional plotting arguments for \code{image}
#'
#' @return NULL
#'
#' @importFrom graphics image
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' # plot model domain
#' fnms <- system.file(c("rflow_mf_demo.bas",
#'                       "rflow_mf_demo.dis"),
#'                     package = "Rflow")
#'
#' dis <- read.DIS(fnms[2L])
#' bas <- read.BAS(fnms[1L], dis)
#'
#' # grid divider co-ordinates (the example has a regular grid)
#' gccs <- with(dis, seq(0, by = DELR, length.out = extent["NCOL"] + 1L))
#' grcs <- with(dis, seq(0, by = DELC, length.out = extent["NROW"] + 1L))
#'
#' MFimage(bas$IBOUND[,, 1L], gccs, grcs, c(-1, 1),
#'         c("blue", "grey", "transparent"))
#'
MFimage <- function(mtx, gccs, grcs, zlim = "auto",
                    col = colorRampPalette(c("blue", "white", "red"))(101L),
                    show.range = identical(zlim, "angle"), ...){
  # case in which matrix is only NA - make blank plot
  if(all(is.na(mtx))){
    cat("all NAs\n")
    image(gccs, grcs, mtx, zlim = 0:1)
    return(invisible())
  }

  # automatic range
  if(identical(zlim, "auto")) zlim <- range(pretty(range(mtx, na.rm = T)))

  # automatic symmetrical range
  if(identical(zlim, "sym")){
    rg1 <- max(abs(pretty(range(mtx, na.rm = T))))
    zlim <- c(-rg1, rg1)
  }

  # variable is an angle
  if(identical(zlim, "angle")){
    col <- colorRampPalette(c("white", "blue", "yellow",
                              "red", "white"))(101L)

    if(show.range){
      cat("white \u2190; blue \u2193; yellow \u2192; red \u2191")
    }

    zlim <- pi*c(-1, 1)
  }else if(show.range){cat(zlim, sep = " to "); cat("; ")}

  # plot image
  image(flip(mtx, 2L), x = gccs, y = grcs, col = col, zlim = zlim, ...)
}

#' Plot MODFLOW variable matrix from NetCDF
#'
#' @param nc
#' NetCDF object;
#' an open connection to a MODFLOW/ MT3DMS data set
#' @param variable
#' character string;
#' name of the variable to be plotted
#' @param start
#' integer \code{[4]};
#' index in each dimension at which to start reading (see
#'  \code{\link{var.get.nc, pkg = RNetCDF}}); \code{NA} implies start
#' @param count
#' integer \code{[4]};
#' numbers of indices to read in each dimension (see
#'  \code{\link{var.get.nc, pkg = RNetCDF}}); \code{NA} implies to end
#' @inheritParams MFimage
#' @param zlim
#' @param col
#' @param show.range
#' @inheritDotParams MFimage
#' @param ...
#'
#' @return NULL
#'
#' @import RNetCDF
#' @export
#'
#' @examples
#'
#' library("RNetCDF")
#' fnm <- system.file("rflow_mf_demo.nc", package = "Rflow")
#' mfdata <- open.nc(fnm)
#'
#' # plot the flux from river (blue is negative)
#' MFncimage(mfdata, "RiverLeakage",
#'           c(NA, NA, 1L, 10L),
#'           c(NA, NA, 1L, 1L),
#'           "sym", show.range = TRUE)
#'
MFncimage <- function(nc, variable,
                      start = c(NA, NA, 1L, 1L),
                      count = c(NA, NA, 1L, 1L), zlim = "auto",
                      col = colorRampPalette(c("blue", "white", "red"))(101L),
                      show.range = identical(zlim, "angle"), ...){
  mtx <- var.get.nc(nc, variable, start, count)
  xy0 <- c(att.get.nc(nc, "NC_GLOBAL", "origin-x"),
           att.get.nc(nc, "NC_GLOBAL", "origin-y"))

  MFimage(mtx,
          var.get.nc(nc, "gccs", start[1L], count[1L]) + xy0[1L],
          var.get.nc(nc, "grcs", start[2L], count[2L]) + xy0[2L],
          zlim, col, show.range, ...)
}
