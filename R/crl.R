# Rflow package - cellref.loc

#' Find grid reference of point in grid
#'
#' Find column, row or other discrete divider reference in a grid, from its
#'  co-ordinate.
#'
#' @param x
#' numeric \code{[]]};
#' locations in relevant dimension
#' @param gcs
#' numeric \code{[]};
#' grid line co-ordinates in relevant dimension.  Must have length one more
#'  than the number of discretised units in the relevant dimension.  For
#'  example, when finding the column reference, the length of gcs must be
#'  one more than the number of columns, including the left and right edges
#'  of the model grid.
#' @param rev
#' logical \code{[1]};
#' whether the grid lines are numbered in reverse order of co-ordinate
#'
#' @return
#' An integer vector of the same length as x, showing which grid cell
#' reference the values of x are found in.  NA shows that values of x are
#' outside the grid boundaries.
#'
#' Values of x that are exactly on the grid lines are assigned to the
#' smaller reference (or greater if rev = TRUE).
#'
#' @export
#'
#' @examples
#' # some of these values may be outside the grid boundary
#' x <- sample(seq(0, 11, .1), 20)
#' # irregular grids are supported, but gcs must be in ascending order
#' gcs <- sort(sample(seq(0, 10, .5), 11))
#'
#' cellref.loc(x, gcs)
#' cellref.loc(x, gcs, TRUE)
#'
#' # for MODFLOW models:
#' library("RNetCDF")
#' mfdata <- open.nc(system.file("rflow_mf_demo.nc", package = "Rflow"))
#' #
#' # - example column reference
#' gccs <- c(var.get.nc(mfdata, "gccs") +
#'             att.get.nc(mfdata, "NC_GLOBAL", "origin-x"))
#' (C <- cellref.loc(625, gccs))
#' #
#' # - example row reference
#' #  -- rows are numbered in order of decreasing y, so put rev = TRUE
#' grcs <- c(var.get.nc(mfdata, "grcs") +
#'             att.get.nc(mfdata, "NC_GLOBAL", "origin-y"))
#' (R <- cellref.loc(825, grcs, TRUE))
#' #
#' # - example layer reference
#' #  -- layer divides are in descending order, so reverse the layer
#' #      divides and also use rev = TRUE
#' ldivs.cr <- c(var.get.nc(mfdata, "elev", c(C, R, NA), c(1L, 1L, NA)))
#' (L <- cellref.loc(25, rev(ldivs.cr), TRUE))
#' #
#' # - example time step reference
#' #  -- remember to prefix the start time
#' mftime <- c(0, var.get.nc(mfdata, "time")) +
#'             att.get.nc(mfdata, "NC_GLOBAL", "start_time")
#' (ts <- cellref.loc(1500, mftime))
#'
cellref.loc <- function(x, gcs, rev = FALSE){
  refs <- as.integer(if(rev){
    length(gcs) - approx(gcs, 1:length(gcs), x, "constant", f = 0)$y
  }else{
    approx(gcs, 1:length(gcs), x, "constant", f = 0)$y
  })
  refs[x >= max(gcs)] <- NA_integer_
  refs
}
