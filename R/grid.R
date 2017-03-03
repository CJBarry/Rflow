# Rflow package - cellref.loc

#' Find grid reference of point in grid
#'
#' Find column, row or other discrete divider reference in a grid, from its
#'  co-ordinate.  May also be used in the time dimension.  Vectorised in
#'  \code{x}.
#'
#' @param x
#' numeric \code{[]};
#' locations in relevant dimension
#' @param gcs
#' numeric \code{[]};
#' grid line co-ordinates in relevant dimension.  Must have length one more
#'  than the number of discretised units in the relevant dimension.  For
#'  example, when finding the column reference, the length of \code{gcs} must be
#'  one more than the number of columns, including the left and right edges
#'  of the model grid.
#' @param rev
#' logical \code{[1]};
#' whether the grid lines are numbered in reverse order of co-ordinate.
#'  Generally \code{FALSE} for column references and \code{TRUE} for row
#'  references.  See the examples for more guidance on usage with MODFLOW
#'  grids.
#'
#' @return
#' An integer vector of the same length as \code{x}, showing which grid cell
#'  reference the values of \code{x} are found in.  \code{NA} shows that values
#'  of \code{x} are outside the grid boundaries.
#'
#' Values of \code{x} that are exactly on the grid lines are assigned to the
#'  smaller reference (or greater if \code{rev = TRUE}, as a side effect of
#'  \code{rev}).
#'
#' @importFrom stats approx
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

#' @name gcs
#' @rdname gcs
#'
#' @title Grid line co-ordinates
#'
#' @description
#' Find the co-ordinates of column (\code{gccs}) or row (\code{grcs})
#'  dividers from the DIS package or the NetCDF-format output data.
#'
#' @param data
#' DIS.MFpackage or NetCDF object;
#' MODFLOW data relating to the model
#' @param absolute
#' logical \code{[1]};
#' \code{TRUE} if the absolute, or global co-ordinates should be returned.
#'  The default is \code{FALSE}, so that the first co-ordinate is always 0
#' @param x0,y0
#' numeric \code{[1]};
#' If \code{data} is a DIS package and \code{absolute = TRUE}, then the
#'  value of \eqn{x} (\code{gccs}) or \eqn{y} (\code{grcs}) at the model
#'  origin (bottom left corner) must be given.  Not needed with NetCDFs
#'  because this information is stored in the NetCDF output generated from
#'  \code{Rflow}.
#'
#' @return
#' \code{gccs}: numeric \code{[NCOL + 1]}\cr
#' \code{grcs}: numeric \code{[NROW + 1]}
#'
#' @examples
#' # using DIS.MFpackage
#' dis <- read.DIS(system.file("rflow_mf_demo.dis",
#'                             package = "Rflow"))
#' gccs(dis)
#' grcs(dis)
#'
#' # arbitrary origin
#' gccs(dis, TRUE, 1000)
#'
#' # using NetCDF
#' mfdata <- RNetCDF::open.nc({
#'   system.file("rflow_mf_demo.nc", package = "Rflow")
#' })
#'
#' gccs(mfdata)
#' grcs(mfdata)
#'
#' # origin is in fact 0,0 for the demo, so these will show no difference
#' gccs(mfdata, TRUE)
#' grcs(mfdata, TRUE)
#'
NULL

#' @rdname gcs
#' @importFrom RNetCDF var.get.nc
#' @importFrom RNetCDF att.get.nc
#' @export
gccs <- function(data, absolute = FALSE, x0){
  switch(class(data),
         NetCDF = {
           if(absolute && !missing(x0)) warning({
             "Rflow::gccs: with NetCDF origin x is read from the data set, so x0 is ignored"
           })

           as.vector(var.get.nc(data, "gccs") + if(absolute){
             att.get.nc(data, "NC_GLOBAL", "origin-x")
           }else 0)
         },
         DIS.MFpackage = {
           if(absolute && missing(x0)) stop({
             "Rflow::gccs: with a DIS.MFpackage, x0 must be given for absolute co-ordinates"
           })

           if(length(data$DELR) == 1L){
             unname(seq(0, by = data$DELR,
                        length.out = data$extent["NCOL"] + 1L)) +
               if(absolute) x0 else 0
           }else{
             c(0, cumsum(data$DELR)) + if(absolute) x0 else 0
           }
         },
         stop("data must be a DIS.MFpackage or a NetCDF"))
}

#' @rdname gcs
#' @importFrom RNetCDF var.get.nc
#' @importFrom RNetCDF att.get.nc
#' @export
grcs <- function(data, absolute = FALSE, y0){
  switch(class(data),
         NetCDF = {
           if(absolute && !missing(y0)) warning({
             "Rflow::grcs: with NetCDF origin y is read from the data set, so y0 is ignored"
           })

           as.vector(var.get.nc(data, "grcs") + if(absolute){
             att.get.nc(data, "NC_GLOBAL", "origin-y")
           }else 0)
         },
         DIS.MFpackage = {
           if(absolute && missing(x0)) stop({
             "Rflow::grcs: with a DIS.MFpackage, y0 must be given for absolute co-ordinates"
           })

           if(length(data$DELC) == 1L){
             unname(seq(0, by = data$DELC,
                        length.out = data$extent["NROW"] + 1L)) +
               if(absolute) y0 else 0
           }else{
             c(0, cumsum(data$DELC)) + if(absolute) y0 else 0
           }
         },
         stop("data must be a DIS.MFpackage or a NetCDF"))
}
