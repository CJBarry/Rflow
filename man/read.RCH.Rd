% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCH.R
\name{read.RCH}
\alias{read.RCH}
\title{Read recharge (RCH) package file}
\usage{
read.RCH(filename, dis, show.prog = FALSE)
}
\arguments{
\item{filename}{character string;
name of file to read}

\item{dis}{object of class DIS.MFpackage;
the DIS package for the model, which gives the array extent}

\item{show.prog}{logical \code{[1]};
whether to print a \code{.} for each stress period that is read, because
 the function can take some time with large models}
}
\value{
a list with class RCH.MFpackage, with elements:\cr
  \code{$header} 1-row data.frame:\cr
  \code{..$NRCHOP} int: recharge option (1, 2 or 3)\cr
  \code{..$IRCHCB} int: unit number to which output data is saved; 0
   signals that it is not saved\cr
  \code{$data} num \code{[NCOL, NROW, NPER]}: the recharge input, in
   units of volume per area per time
}
\description{
Reads information from a MODFLOW RCH package file.  The RCH package
 specifies a distributed fixed flux per unit area.
}
\examples{

}

