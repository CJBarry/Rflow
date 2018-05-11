# Rflow package - functions for the BAS MODFLOW package

# read and write BAS package ----------------------------------------------

#' Read basic (BAS) package file
#'
#' Reads information from a MODFLOW-2000 BAS package file.  The BAS
#'  package contains information about the active region, the head assigned
#'  to no-flow cells and starting heads.
#'
#' @param filename
#' character string;
#' file to be read
#' @param dis
#' character string or object of class DIS.MFpackage;
#' a DIS package file (as file path or object), from which model dimensions
#'  will be read
#' @param MF96
#' logical [1];
#' whether to expect an MF96 formatted BAS file
#'
#' @return
#' object of class BAS.MFpackage:\cr
#'   \code{$Options}: chr \code{[]}; options such as SHOWPROGRESS, FREE\cr
#'   \code{$IBOUND} : num \code{[NCOL, NROW, NLAY]} or \code{[NLAY]}; flags for
#'    inactive (0), active (1) and fixed head cells (-1)\cr
#'   \code{$HNOFLO} : num \code{[1]}; head assigned to no flow cells\cr
#'   \code{$STRT}   : num \code{[NCOL, NROW, NLAY]} or \code{[NLAY]}; initial
#'    heads\cr
#'
#' @export
#'
#' @examples
#'
read.BAS <- function(filename, dis, MF96 = FALSE){
  dis <- switch(class(dis),
                character = read.DIS(dis),
                DIS.MFpackage = dis,
                stop("invalid dis to read.BAS"))

  txt <- readLines(filename)

  # remove comment/ title lines
  if(MF96) txt <- txt[-(1:2)] else txt <- txt[!is.comment(txt, "#")]

  # options
  Options <- if(!MF96) scan(text = txt[1L], what = "character", quiet = TRUE)

  # extent (only for MF96 - otherwise this information is in the DIS package)
  if(MF96){
    extent <- as.integer(read.fws(txt[1L]))
    names(extent) <- c("NLAY", "NROW", "NCOL", "NPER", "ITMUNI", rep("", length(extent) - 5L))
  }else extent <- NULL

  # unit numbers (only for MF96 - otherwise this information is in the NAM file)
  unitNos <- if(MF96) as.integer(read.fws(txt[2L], 3L))

  # IBOUND array
  ln <- 2L
  IBOUND <- replicate(dis$extent["NLAY"], {
    nIB <- expected.lines.RIARRAY(txt[ln], dis$extent["NCOL"], dis$extent["NROW"], 1L)
    lns <- seq(ln, ln + nIB, 1L)
    ar <- interpret.RIARRAY(txt[lns], T, dis$extent["NCOL"], 1L)
    ln <<- ln + nIB + 1L
    ar
  }, simplify = FALSE)

  lgs <- lengths(IBOUND)

  if(all(lgs == 1L)){
    IBOUND <- unname(c(IBOUND, recursive = TRUE))
  }else{
    ar <- array(dim = dis$extent[c("NCOL", "NROW", "NLAY")])
    for(l in 1:length(IBOUND)) ar[,, l] <- IBOUND[[l]]
    IBOUND <- ar; rm(ar)
  }

  #no flow cell head label
  HNOFLO <- as.numeric(txt[ln])

  #starting head
  ln <- ln + 1L
  STRT <- replicate(dis$extent["NLAY"], {
    nST <- expected.lines.RIARRAY(txt[ln], dis$extent["NCOL"], dis$extent["NROW"], 1L)
    lns <- ln:(ln + nST)
    ar <- interpret.RIARRAY(txt[lns], T, dis$extent["NCOL"], 1L)
    ln <<- ln + nST + 1L
    ar
  }, simplify = FALSE)

  lgs <- lengths(STRT)

  if(all(lgs == 1L)){
    STRT <- unname(c(STRT, recursive = TRUE))
  }else{
    ar <- array(dim = dis$extent[c("NCOL", "NROW", "NLAY")])
    for(l in 1:length(STRT)) ar[,, l] <- STRT[[l]]
    STRT <- ar; rm(ar)
  }

  structure(mget(c("Options", "IBOUND", "HNOFLO", "STRT")),
            class = "BAS.MFpackage")
}

#' Write MODFLOW BAS package file
#'
#' Writes information from a BAS.MFpackage list object to a MODFLOW-readable
#'  BAS package file.
#'
#' @param BAS
#' object of class BAS.MFpackage, as would be returned by
#'  \code{\link{read.BAS}}
#' @param filename
#' character string;
#' file to write to, ideally ending in \code{".bas"}
#' @param unit
#' integer \code{[1]};
#' unit number that the BAS file will be assigned to
#' @param title
#' character string;
#' optional title to prefix to the file, for your reference (not used if MF2k =
#'  FALSE)
#' @param MF2k
#' logical;
#' should the file be for MODFLOW 2000, or MODFLOW 96 (in which case stress
#'  period information is written from sps)
#' @param sps
#' data.frame or object of class DIS.MFpackage;
#' only required if \code{MF2k = FALSE}; a data frame of stress period
#'  information with the same columns as described in \code{\link{read.DIS}},
#'  except that the \code{TR} column is not required (ignored if provided)
#' @param NCOL,NROW,NLAY
#' integer \code{[1]}; only required if \code{MF2k = FALSE} and
#'   \code{class(sps) != "DIS.MFpackage"}
#'
#' @return \code{NULL}
#'
#' @import data.table
#' @export
#'
#' @examples
#'
write.BAS <- function(BAS, filename, unit = 1L, title = "", MF2k = TRUE, sps,
                      NCOL, NROW, NLAY){
  stopifnot(is(BAS, "BAS.MFpackage"))

  if(!MF2k){
    if(class(sps) == "DIS.MFpackage"){
      NCOL <- sps$extent[c("NCOL")]
      NROW <- sps$extent[c("NROW")]
      NLAY <- sps$extent[c("NLAY")]
    }

    sps <- switch(class(sps)[1L],
                  data.frame = sps,
                  data.table = sps,
                  DIS.MFpackage = sps$sps,
                  stop("write.BAS: invalid sps"))

    setDF(sps)
    sps <- sps[, c("PERLEN", "NSTP", "TSMULT")]
  }

  con <- file(filename, "wt")
  on.exit(close(con))

  writeLines(paste("#", {
    substr(gsub("\n", "", title), 1L, if(MF2k) 1000L else 78L)
  }), con)

  writeLines(c(if(MF2k){
    "# MODFLOW 2000 basic package file created by write.BAS function, in R"
  }else{
    "# MODFLOW 1996 basic package file created by write.BAS function, in R"
  }), con)

  # MF96 item 3
  # - NLAY, NROW, NCOL, NPER, ITMUNI
  if(!MF2k){
    writeLines(paste0(formatC(NLAY, 0L, 10L, "d"),
                      formatC(NROW, 0L, 10L, "d"),
                      formatC(NCOL, 0L, 10L, "d"),
                      formatC(nrow(sps), 0L, 10L, "d"),
                      formatC(0L, 0L, 10L, "d")), con)
  }

  # MF2k item 1
  # MF96 item 4
  # - Options
  if(!all(BAS$Options %in% c("XSECTION", "CHTOCH",
                             "FREE", "SHOWPROGRESS"))){
    warning("write.BAS: invalid BAS options are being written to BAS file")
  }
  writeLines(paste(BAS$Options, collapse = " "), con)

  # MF96 item 5
  # - IAPART, ISTRT
  if(!MF2k){
    writeLines("         1         1", con)
  }

  # MF2k item 2
  # MF96 item 6
  # - IBOUND array
  writeLines(if(is.vector(BAS$IBOUND)){
    RIARRAY.splitlayers(length(BAS$IBOUND), CNSTNT = BAS$IBOUND,
                        compress = TRUE, FMTIN_type = "I", FMTIN_w = 3L,
                        flag.no = unit)
  }else{
    RIARRAY.splitlayers(dim(BAS$IBOUND)[3L], arr = BAS$IBOUND,
                        compress = TRUE, FMTIN_type = "I", FMTIN_w = 3L,
                        flag.no = unit)
  }, con)

  # MF2k item 3
  # MF96 item 7
  # - HNOFLO
  writeLines(formatC(BAS$HNOFLO, 2L, 10L, "e"), con)

  # MF2k item 4
  # MF96 item 8
  # - STRT/ Shead
  writeLines(if(is.vector(BAS$STRT)){
    RIARRAY.splitlayers(length(BAS$STRT), CNSTNT = BAS$STRT,
                        compress = TRUE, FMTIN_type = "e", FMTIN_w = 12L,
                        FMTIN_d = 4L, flag.no = unit)
  }else{
    RIARRAY.splitlayers(dim(BAS$STRT)[3L], arr = BAS$STRT,
                        compress = TRUE, FMTIN_type = "e", FMTIN_w = 12L,
                        FMTIN_d = 4L, flag.no = unit)
  }, con)

  # MF96 item 9
  # - stress period information
  # - PERLEN, NSTP, TSMULT
  if(!MF2k){
    setDT(sps)
    on.exit(setDF(sps), add = TRUE)
    ff.widths <- c(10L, 10L, 10L)
    ff.digit <- c(2L, 0L, 2L)
    ff.format <- c("e", "d", "e")
    sps[, {

      ffs <- matrix("", .N, 4L)
      ffs[, 1:3] <- vapply(1:3,
                           function(col){
                             formatC(.SD[[col]], width = 10L,
                                     digits = ff.digit[col],
                                     format = ff.format[col])
                           }, character(.N))

      writeLines(apply(ffs, 1L, paste, collapse = ""), con)
    }]
  }

  invisible()
}
