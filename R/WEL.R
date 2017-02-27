# Rflow package - functions for the WEL MODFLOW package

#' Read well (WEL) package file
#'
#' Reads information from a MODFLOW WEL package file.  The WEL package
#'  assigns fixed sink/ source fluxes to individual cells and as such is
#'  particularly appropriate for representing wells of known abstraction or
#'  injection rate.
#'
#' @param filename
#' character string;
#' WEL package file
#' @param nSP
#' integer \code{[1]} or DIS.MFpackage object;
#' number of stress periods or corresponding DIS package from which this
#'  can be read
#' @param MF2k
#' logical \code{[1]};
#' is this a MODFLOW 2000+ WEL package?
#'
#' @return
#' list \code{[2]} of class WEL.MFpackage, with elements:\cr
#' \code{$header}:\cr
#'   \code{..$MXACTW} the maximum number of active wells in a stress period\cr
#'   \code{..$IWELCB} unit number to which well flux array is saved (if
#'    \eqn{> 0})\cr
#' \code{$spheaders}:\cr
#' headers for each stress period, a data.frame with columns:\cr
#'   \code{..$ITMP} (int): number of active wells in stress period, or (if
#'    \eqn{< 0})\cr
#'   indicates to read from previous stress period\cr
#'   \code{..$NP} (int): number of parameters to read for this stress period
#'    (always 0 if \code{MF2k = FALSE})\cr
#' \code{$read}: logical vector of length nSP indicating for which stress periods
#'  well information is to be read\cr
#' \code{$data}:\cr
#' a data.table with the following columns\cr
#'   \code{..$sp} (int): stress period (key)\cr
#'   \code{..$C} (int): column\cr
#'   \code{..$R} (int): row\cr
#'   \code{..$L} (int): layer\cr
#'   \code{..$Q} (num): flux (negative indicates abstraction, positive indicates
#'    injection)\cr
#'
#' @import data.table
#' @export
#'
#' @examples
#' fnms <- system.file(c("rflow_mf_demo.dis",
#'                       "rflow_mf_demo.wel"), package = "Rflow")
#'
#' # get model information from DIS package file
#' dis <- read.DIS(fnms[1L])
#'
#' wel <- read.WEL(fnms[2L], dis)
#' # or, if you already know how many stress periods there are
#' wel <- read.WEL(fnms[2L], 15L)
#'
#' class(wel)
#' str(wel)
#'
read.WEL <- function(filename, nSP, MF2k = TRUE){
  # get number of stress periods
  nSP <- switch(class(nSP),
                integer = nSP,
                numeric = as.integer(nSP),
                DIS.MFpackage = nrow(nSP$sps),
                stop("number of stress periods not given"))

  if(MF2k){
    con <- file(filename, "rt")
    on.exit(close(con))

    # find how many lines to skip
    # note that this function has not been programmed to read parameter
    #  wells
    skip <- 0L
    while({
      ln1 <- readLines(con, 1L)
      Rflow:::is.comment(ln1, "#") || grepl("PARAMETER", ln1, ignore.case = TRUE)
    }) skip <- skip + 1L

    # start over and skip
    close(con)
    on.exit(NULL)
  }else skip <- 0L

  con <- file(filename, "rt")
  on.exit(close(con))
  if(skip) readLines(con, skip)

  # read global header
  hd <- read.fwf(con, n = 1L, widths = c(10L, 10L),
                 col.names = c("MXACTW", "IWELCB"))

  read <- logical(nSP)
  spheaders <- data.frame(ITMP = integer(nSP),
                          NP = integer(nSP))
  lst <- vector("list", nSP)

  # for each stress period, read ...
  for(sp in 1:nSP){
    sphd <- read.fwf(con, n = 1L, widths = c(10L, if(MF2k) 10L),
                     col.names = c("ITMP", if(MF2k) "NP"))

    spheaders$ITMP[sp] <- sphd$ITMP
    if(MF2k && sphd$NP){
      stop("read.WEL is not yet programmed to read parameter wells\n")
    }

    read[sp] <- as.logical(sphd$ITMP)

    if(read[sp]){
      lst[[sp]] <- read.fwf(con, n = sphd$ITMP, widths = rep(10L, 4L),
                            col.names = c("L", "R", "C", "Q"))
      setDT(lst[[sp]])
      lst[[sp]][, sp := sp]
    }
  }

  data <- rbindlist(lst)
  setkey(data, sp)
  setcolorder(data, c("sp", "C", "R", "L", "Q"))

  structure(list(header = hd, spheaders = spheaders, read = read,
                 data = data), class = "WEL.MFpackage")
}

#' Write a well (WEL) package file
#'
#' Writes information from a WEL.MFpackage list object to a MODFLOW-readable
#'  WEL package file.
#'
#' @param WEL
#' object of class WEL.MFpackage
#' @param filename
#' character string;
#' file to which package information is to be written
#' @inheritParams write.BAS
#' @param title
#' @param MF2k
#'
#' @return
#' NULL
#'
#' @import data.table
#' @export
#'
#' @examples
#' # read WEL package
#' fnm <- system.file("rflow_mf_demo.wel", package = "Rflow")
#' wel <- read.WEL(fnm, 15L)
#'
#' # meaningful modification e.g. add a stress period with no abstraction
#' wel$spheader <- rbind(wel$spheader, list(0L, 0L))
#' wel$read <- c(wel$read, TRUE)
#'
#' # write modified wel package to file
#' write.WEL(wel, "RFLOW_EXAMPLE.wel", "example: added stress period")
#'
write.WEL <- function(WEL, filename, title = "", MF2k = TRUE){
  stopifnot(identical(class(WEL), "WEL.MFpackage"))

  con <- file(filename, "wt")
  on.exit(close(con))

  # title and PARAMETER specifications
  if(MF2k) writeLines(c({
    "# MODFLOW 2000 stream package file created by write.WEL function, in R"
  }, paste("#", title), {
    "PARAMETER  0  0"
  }), con)

  # write header
  writeLines(c(formatC(max(WEL$data[, .N, by = sp]$N), digits = 0L,
                       width = 10L, format = "d"),
               formatC(WEL$header$IWELCB[1L], digits = 0L,
                       width = 10L, format = "d")),
             con, sep = "")
  writeLines("", con)

  # number of stress periods
  nSP <- length(WEL$read)

  for(spn in 1:nSP){
    if(WEL$read[spn]){
      WEL$data[sp == spn & Q != 0, if(!.N){
        writeLines("         0         0", con)
      }else{
        writeLines(c(formatC(.N, digits = 0L, width = 10L, format = "d"),
                     "         0"),
                   con, sep = "")
        writeLines("", con)

        # segments
        sdt <- .SD[, list(L, R, C, Q)]
        ff.width <- rep(10L, 4L)
        ff.digits <- c(0L, 0L, 0L, 3L)
        ff.format <- c("d", "d", "d", "f")

        ffs <- vapply(1:4,
                      function(col){
                        formatC(sdt[[col]],
                                width = ff.width[col],
                                digits = ff.digits[col],
                                format = ff.format[col])
                      }, character(.N))

        if(is.matrix(ffs)){
          writeLines(apply(ffs, 1L, paste, collapse = ""), con)
        }else{
          writeLines(paste(ffs, collapse = ""), con)
        }
      }]
    }else{
      writeLines("        -1         0", con)
    }
  }

  invisible()
}
