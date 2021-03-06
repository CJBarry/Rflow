% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DIS.R
\name{read.DIS}
\alias{read.DIS}
\title{Read a discretisation (DIS) package file}
\usage{
read.DIS(file)
}
\arguments{
\item{file}{character string; file name}
}
\value{
a list of class "DIS.MFpackage" with elements:\cr
\code{$extent}: a vector whose named elements are \code{NLAY},
 \code{NROW}, \code{NCOL}, \code{NPER}, \code{t_unit}, \code{l_unit}\cr
$LAYCBD: integer vector indicating the equation type used for each layer
  (confined or unconfined, see MODFLOW-2000 documentation)\cr
\code{$DELR}: column spacing (along rows), may be a single constant value\cr
\code{$DELC}: row spacing (along columns), may be a single constant value\cr
\code{$elev}: a vector of layer divide elevations or a 3D array of distributed
  layer divide elevations by cell\cr
\code{$sps}: a data.frame with stress period descriptions:\cr
\code{..$PERLEN} (num): stress period length\cr
\code{..$NSTP} (int): number of time steps in the stress period\cr
\code{..$TSMULT} (num): time step multiplier\cr
\code{..$TR} (log): does the stress period use transient equations?\cr
}
\description{
Reads information from a MODFLOW-2000 DIS package file.  The DIS
 package contains information about the finite difference grid,
 including layer elevations, and stress period set up.  It also tells
 the layer types.
}
\examples{
fnm <- system.file("rflow_mf_demo.dis", package = "Rflow")

dis <- read.DIS(fnm)
class(dis)
str(dis)

}
