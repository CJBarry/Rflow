% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid.R
\name{cellref.loc}
\alias{cellref.loc}
\title{Find grid reference of point in grid}
\usage{
cellref.loc(x, gcs, rev = FALSE)
}
\arguments{
\item{x}{numeric \code{[]};
locations in relevant dimension}

\item{gcs}{numeric \code{[]};
grid line co-ordinates in relevant dimension.  Must have length one more
 than the number of discretised units in the relevant dimension.  For
 example, when finding the column reference, the length of \code{gcs} must be
 one more than the number of columns, including the left and right edges
 of the model grid.}

\item{rev}{logical \code{[1]};
whether the grid lines are numbered in reverse order of co-ordinate.
 Generally \code{FALSE} for column references and \code{TRUE} for row
 references.  See the examples for more guidance on usage with MODFLOW
 grids.}
}
\value{
An integer vector of the same length as \code{x}, showing which grid cell
 reference the values of \code{x} are found in.  \code{NA} shows that values
 of \code{x} are outside the grid boundaries.

Values of \code{x} that are exactly on the grid lines are assigned to the
 smaller reference (or greater if \code{rev = TRUE}, as a side effect of
 \code{rev}).
}
\description{
Find column, row or other discrete divider reference in a grid, from its
 co-ordinate.  May also be used in the time dimension.  Vectorised in
 \code{x}.
}
\examples{
# some of these values may be outside the grid boundary
x <- sample(seq(0, 11, .1), 20)
# irregular grids are supported, but gcs must be in ascending order
gcs <- sort(sample(seq(0, 10, .5), 11))

cellref.loc(x, gcs)
cellref.loc(x, gcs, TRUE)

# for MODFLOW models:
library("RNetCDF")
mfdata <- open.nc(system.file("rflow_mf_demo.nc", package = "Rflow"))
#
# - example column reference
gccs <- c(var.get.nc(mfdata, "gccs") +
            att.get.nc(mfdata, "NC_GLOBAL", "origin-x"))
(C <- cellref.loc(625, gccs))
#
# - example row reference
#  -- rows are numbered in order of decreasing y, so put rev = TRUE
grcs <- c(var.get.nc(mfdata, "grcs") +
            att.get.nc(mfdata, "NC_GLOBAL", "origin-y"))
(R <- cellref.loc(825, grcs, TRUE))
#
# - example layer reference
#  -- layer divides are in descending order, so reverse the layer
#      divides and also use rev = TRUE
ldivs.cr <- c(var.get.nc(mfdata, "elev", c(C, R, NA), c(1L, 1L, NA)))
(L <- cellref.loc(25, rev(ldivs.cr), TRUE))
#
# - example time step reference
(ts <- cellref.loc(1500, mftstime(mfdata)))

}
