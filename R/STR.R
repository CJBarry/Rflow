# Rflow package - functions for the STR MODFLOW package

#' Read stream (STR) package file
#'
#' Reads information from a MODFLOW STR package file.  The stream package
#'  routes stream flow through connected cells and simultaneously
#'  calculates water exchange between streams and the aquifer.
#'
#' @inheritParams read.GHB
#'
#' @return
#' list with class STR.MFpackage:\cr
#' \code{$header} 1-row data.frame:\cr
#'   \code{..$MXSTRM} (int): maximum number of stream reaches\cr
#'   \code{..$NSS} (int): maximum number of segments\cr
#'   \code{..$NTRIB} (int): maximum number of tributary segments that join
#'    one downstream segment\cr
#'   \code{..$NDIV} (int): flag, which when positive, specifies that
#'    diversions from segments are to be simulated\cr
#'   \code{..$ICALC} (int): flag, which when positive, specifies that
#'    stream stages in reaches are to be calculated\cr
#'   \code{..$CONST} (num): Constant value used in calculating stream stage in
#'    reaches. It is specified whenever ICALC is greater than zero.  A unit
#'    conversion between stream flow units and the units used for the MODFLOW
#'    simulation.\cr
#'   \code{..$ISTCB1} (int): unit number to which to save stream-aquifer
#'    exchange\cr
#'   \code{..$ISTCB2} (int): unit number to which to save stream flow \cr
#' \code{$spheaders} data.frame with 1 row for each stress period:\cr
#'   \code{..$ITMP} (int): if < 0, read stream data from the previous
#'    stress period, otherwise it is the number of reaches to read for the
#'    current stress period\cr
#'   \code{..$IRDFLG} (int): if > 0, don't print input data for this sp\cr
#'   \code{..$IPTFLG} (int): if > 0, don't print results for this sp\cr
#' \code{$data} data.table:\cr
#'   \code{..$sp} (int); stress period (the key)\cr
#'   \code{..$L} (int): layer\cr
#'   \code{..$R} (int): row\cr
#'   \code{..$C} (int): column\cr
#'   \code{..$Seg} (int): number assigned to a group of reaches\cr
#'   \code{..$Reach} (int): sequential number in a segment that begins with one
#'    for the farthest upstream reach and continues in downstream order to the
#'    last reach in the segment\cr
#'   \code{..$Flow} (num): streamflow, in length cubed per time, entering
#'    a segment\cr
#'   \code{..$Stage} (num): stream stage, in units of length\cr
#'   \code{..$Cond} (num): streambed hydraulic conductance\cr
#'   \code{..$Sbot} (num): elevation of the bottom of the streambed\cr
#'   \code{..$Stop} (num): elevation of the top of the streambed\cr
#'   if \code{ICALC > 0}\{\cr
#'   \code{..$Width} (num): width of the stream channel\cr
#'   \code{..$Slope} (num): slope of the stream channel in each reach\cr
#'   \code{..$Rough} (num): Manning's roughness coefficient for each stream
#'    reach\cr
#'   \}\cr
#'   if \code{NTRIB > 0}\{\cr
#'   \code{..$ITrib1} (int): number of the first tributary segment\cr
#'   ...\cr
#'   \code{..$ITrib<NTrib>} (int):\cr
#'   \}\cr
#'   if \code{NDIV > 0}\{\cr
#'   \code{..$Iupseg} (int): number of the upstream segment from which
#'    water is diverted\cr
#'   \}
#'
#' @import data.table
#' @export
#'
#' @references
#' Prudic, D.E., 1989. Documentation of a Computer Program to Simulate Stream-Aquifer Relations Using a Modular, Finite-Difference Ground-water Flow Model, Carson City, Nevada. Available at: https://pubs.er.usgs.gov/publication/ofr88729.
read.STR <- function(filename, nSP){
  # get number of stress periods if not given already
  nSP <- switch(class(nSP),
                integer = nSP,
                numeric = as.integer(nSP),
                DIS.MFpackage = nrow(nSP$sps),
                stop("Rflow::read.STR: invalid nSP"))

  con <- file(filename, "rt")
  on.exit(close(con))

  # find how many lines to skip
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
  hd <- read.fwf(con, n = 1L, widths = rep(10L, 8L),
                 col.names = c("MXSTRM", "NSS", "NTRIB", "NDIV",
                               "ICALC", "CONST", "ISTCB1", "ISTCB2"))

  spheaders <- data.table(ITMP = integer(nSP),
                          IRDFLG = integer(nSP),
                          IPTFLG = integer(nSP))

  # read data sets for each stress period
  lst <- vector("list", nSP)
  for(sp in 1:nSP){
    sphd <- read.fwf(con, n = 1L, widths = rep(10L, 3L))

    spheaders[sp] <- sphd

    rd <- sphd[, 1L] > 0L
    nrw <- sphd[, 1L]

    lst[[sp]] <- if(rd){
      dt <- data.table(sp = rep(sp, nrw), key = "sp")

      dt[, {
        c("L", "R", "C", "Seg",
          "Reach", "Flow", "Stage",
          "Cond", "Sbot", "Stop")
      } := read.fwf(con, n = nrw,
                    widths = c(5L, 5L, 5L, 5L, 5L,
                               15L, 10L, 10L, 10L, 10L))]


      if(hd$ICALC){
        dt[, c("Width", "Slope", "Rough") := {
          read.fwf(con, n = nrw, widths = c(10L, 10L, 10L))
        }]
      }

      if(hd$NTRIB){
        dt[, paste0("ITrib", 1:hd$NTRIB) := {
          read.fwf(con, n = nrw, widths = rep(5L, hd$NTRIB))
        }]
      }

      if(hd$NDIV){
        dt[, Iupseg := {
          read.fwf(con, n = nrw, widths = 10L)
        }]
      }

      dt
    }else lst[[sp - 1L]]
  }

  structure(list(header = hd, spheaders = spheaders,
                 data = rbindlist(lst)), class = "STR.MFpackage")
}

#' Write a stream (STR) package file
#'
#' Writes information from a STR.MFpackage list object to a MODFLOW-readable
#'  STR package file.
#'
#' @param STR
#' STR.MFpackage object
#' @param filename
#' character string
#'
#' @return NULL
#'
#' @import data.table
#' @export
#'
#' @examples
write.STR <- function(STR, filename){
  stopifnot(identical(class(STR), "STR.MFpackage"))

  con <- file(filename, "wt")
  on.exit(close(con))

  # title and PARAMETER specifications
  writeLines(c({
    "# MODFLOW 2000 stream package file created by write.STR function, in R"
  }, {
    "PARAMETER  0  0"
  }), con)

  hd <- mapply(formatC, STR$header,
               digits = c(rep(0L, 5L), 4L, rep(0L, 2L)),
               width = rep(10L, 8L),
               format = c(rep("d", 5L), "f", rep("d", 2L)))

  writeLines(hd, con, sep = "")
  writeLines("", con)

  STR$data[, {
    sphd <- mapply(formatC, STR$spheaders[sp,],
                   digits = rep(0L, 3L),
                   width = rep(10L, 3L),
                   format = rep("d", 3L))

    writeLines(sphd, con, sep = "")
    writeLines("", con)

    if(STR$spheaders[sp, ITMP > 0L]){
      # segments
      sdt <- .SD[, list(L, R, C, Seg, Reach, Flow,
                        Stage, Cond, Sbot, Stop)]
      ff.width <- c(5L, 5L, 5L, 5L, 5L,
                    15L, 10L, 10L, 10L, 10L)
      ff.digits <- c(0L, 0L, 0L, 0L, 0L,
                     4L, 3L, 3L, 3L, 3L)
      ff.format <- c("d", "d", "d", "d", "d",
                     "f", "f", "f", "f", "f")

      ffs <- vapply(1:10,
                    function(col){
                      formatC(sdt[[col]],
                              width = ff.width[col],
                              digits = ff.digits[col],
                              format = ff.format[col])
                    }, character(.N))

      writeLines(apply(ffs, 1L, paste, collapse = ""), con)

      # stagepars
      if(STR$header$ICALC){
        sdt <- .SD[, list(Width, Slope, Rough)]
        ff.width <- c(10L, 10L, 10L)
        ff.digits <- c(3L, 3L, 3L)
        ff.format <- c("f", "f", "f")

        ffs <- vapply(1:3,
                      function(col){
                        formatC(sdt[[col]],
                                width = ff.width[col],
                                digits = ff.digits[col],
                                format = ff.format[col])
                      }, character(.N))

        writeLines(apply(ffs, 1L, paste, collapse = ""), con)
      }

      # tributaries
      if(ntrib <- STR$header$NTRIB){
        sdt <- .SD[, .SD, .SDcols = paste0("ITrib", 1:ntrib)]
        ff.width <- rep(5L, ntrib)
        ff.digits <- rep(0L, ntrib)
        ff.format <- rep("d", ntrib)

        ffs <- vapply(1:ntrib,
                      function(col){
                        formatC(sdt[[col]],
                                width = ff.width[col],
                                digits = ff.digits[col],
                                format = ff.format[col])
                      }, character(.N))

        writeLines(apply(ffs, 1L, paste, collapse = ""), con)
      }

      # diversions
      if(STR$header$NDIV){
        ff.width <- 10L
        ff.digits <- 0L
        ff.format <- "d"

        ffs <- formatC(Iupseg,
                       width = ff.width,
                       digits = ff.digits,
                       format = ff.format)

        writeLines(ffs, con)
      }
    }
  }, by = sp]

  invisible()
}

