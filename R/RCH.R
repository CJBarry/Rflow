# Rflow package - functions for RCH MODFLOW package

#' Read recharge (RCH) package file
#'
#' Reads information from a MODFLOW RCH package file.  The RCH package
#'  specifies a distributed fixed flux per unit area.
#'
#' @param filename
#' character string;
#' name of file to read
#' @param dis
#' object of class DIS.MFpackage;
#' the DIS package for the model, which gives the array extent
#' @param show.prog
#' logical \code{[1]};
#' whether to print a \code{.} for each stress period that is read, because
#'  the function can take some time with large models
#'
#' @return
#' a list with class RCH.MFpackage, with elements:\cr
#'   \code{$header} 1-row data.frame:\cr
#'   \code{..$NRCHOP} int: recharge option (1, 2 or 3)\cr
#'   \code{..$IRCHCB} int: unit number to which output data is saved; 0
#'    signals that it is not saved\cr
#'   \code{$data} num \code{[NCOL, NROW, NPER]}: the recharge input, in
#'    units of volume per area per time
#'
#' @import data.table
#' @export
#'
#' @examples
#'
read.RCH <- function(filename, dis, show.prog = FALSE){
  # get number of stress periods
  nSP <- nrow(dis$sps)

  con <- file(filename, "rt")
  on.exit(close(con))

  # find how many lines to skip
  # note that this function has not been programmed to read recharge
  #  parameters
  skip <- 0L
  while({
    ln1 <- readLines(con, 1L)
    is.comment(ln1, "#") || grepl("PARAMETER", ln1, ignore.case = TRUE)
  }) skip <- skip + 1L

  # start over and skip
  close(con)
  on.exit(NULL)
  con <- file(filename, "rt")
  on.exit(close(con))
  if(skip) readLines(con, skip)

  # read global header
  hd <- read.fwf(con, n = 1L, widths = c(10L, 10L),
                 col.names = c("NRCHOP", "IRCHCB"))

  # assign array
  RCHar <- array(dim = dis$extent[c("NCOL", "NROW", "NPER")])

  for(sp in 1:nSP){
    sphd <- read.fwf(con, n = 1L, widths = c(10L, 10L),
                     col.names = c("INRECH", "INIRCH"))

    if(sphd$INRECH >= 0L){
      arhd <- readLines(con, 1L)
      arlns <- expected.lines.RIARRAY(arhd, dis$extent["NCOL"],
                                      dis$extent["NROW"], 1L)

      artxt <- c(arhd, if(arlns > 0L) readLines(con, arlns))
      RCHar[,, sp] <-
        unname(interpret.RIARRAY(artxt, TRUE, dis$extent["NCOL"]))
    }else{
      # use data from previous stress period
      RCHar[,, sp] <- RCHar[,, sp - 1L]
    }

    if(show.prog) cat(".")
  }; cat("\n")

  structure(list(header = hd, data = RCHar), class = "RCH.MFpackage")
}
