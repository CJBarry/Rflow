# Rflow package - functions for the HFB MODFLOW package

#' Read a horizontal flow barrier (HFB) package file
#'
#' Reads information from a MODFLOW-2000 HFB6 package file.
#'  The HFB package places hydraulic barriers of set conductance between
#'  specified horizontally adjacent cells and is often used to represent
#'  artificial confining walls or geological faults.
#'
#' @param filename
#' character string
#'
#' @return
#' list with class HFB.MFpackage, with elements:\cr
#' \code{$header}, 1-row data.frame:\cr
#' \code{..$NPHFB} int: number of HFBs specified as parameters (note that
#'  this feature is not yet supported by \code{read.HFB})\cr
#' \code{..$MXFB} int: the maximum number of flow barriers that will be
#'  defined by parameters (again, not yet supported)\cr
#' \code{..$NHFBNP} int: the number of barriers not defined by paramters\cr
#' \code{$data}, data.table:\cr
#' \code{..$L} int: layer\cr
#' \code{..$IROW1} int: row on side 1 of the barrier\cr
#' \code{..$ICOL1} int: column on side 1 of the barrier\cr
#' \code{..$IROW2} int: row on side 2 of the barrier\cr
#' \code{..$ICOL2} int: column on side 2 of the barrier\cr
#' \code{..$Hydchr} num: Hydraulic characteristic of the barrier, which is
#'  the hydraulic conductivity divided by the width
#'
#' @import data.table
#' @export
#'
#' @examples
read.HFB <- function(filename){
  con <- file(filename, "rt")
  on.exit(close(con))

  # is the first line a comment?
  com <- is.comment(readLines(con, 1L), "#")

  close(con)
  on.exit(NULL)

  con <- file(filename, "rt")
  on.exit(close(con))
  if(com) readLines(con, 1L)

  # free format read
  hd <- read.table(con, nrows = 1L,
                   col.names = c("NPHFB", "MXFB", "NHFBNP"))

  if(hd$NPHFB > 0L) stop({
    "some barriers are defined by parameters - this feature needs to be added to read.HFB"
  })

  # free format read
  # assuming steady state data
  data <- data.table(read.table(con, nrows = hd$NHFBNP,
                                col.names = c("L", "IROW1", "ICOL1",
                                              "IROW2", "ICOL2", "Hydchr")))

  structure(list(header = hd, data = data), class = "HFB.MFpackage")
}

#' Write a hydraulic flow barrier (HFB) package file
#'
#' Writes information from a HFB.MFpackage list object to a MODFLOW-readable
#'  HFB package file.
#'
#' @param HFB
#' HFB.MFpackage object, as would be returned by \code{\link{read.HFB}}
#' @param filename
#' character string
#' @param title
#' character string;
#' optional title to put into package file as a comment
#'
#' @return NULL
#'
#' @import data.table
#' @export
#'
#' @examples
write.HFB <- function(HFB, filename, title = ""){
  con <- file(filename, "wt")
  on.exit(close(con))

  writeLines({
    c("# MODFLOW 2000 HFB package file created by write.HFB function, in R",
      paste("#", title))
  }, con)

  writeLines(as.character(c(HFB$header, recursive = TRUE)),
             con, sep = "  ")
  writeLines("", con)

  HFB$data[, {
    ff.digit <- c(rep(0L, 5L), 3L)
    ff.format <- c(rep("d", 5L), "e")

    ffs <- vapply(1:6,
                  function(col){
                    formatC(.SD[[col]], width = 10L,
                            digits = ff.digit[col],
                            format = ff.format[col])
                  }, character(.N))

    writeLines(apply(ffs, 1L, paste, collapse = ""), con)
  }, .SDcols = c("L", "IROW1", "ICOL1", "IROW2", "ICOL2", "Hydchr")]

  writeLines("  0", con)

  invisible()
}

