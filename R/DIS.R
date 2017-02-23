# Rflow package - functions for the DIS MODFLOW package

#' Read a discretisation (DIS) package file
#'
#' Reads information from a MODFLOW-2000 DIS package file.  The DIS
#'  package contains information about the finite difference grid,
#'  including layer elevations, and stress period set up.  It also tells
#'  the layer types.
#'
#' @param file
#' character string; file name
#'
#' @return
#' a list of class "DIS.MFpackage" with elements:\cr
#' \code{$extent}: a vector whose named elements are \code{NLAY},
#'  \code{NROW}, \code{NCOL}, \code{NPER}, \code{t_unit}, \code{l_unit}\cr
#' $LAYCBD: integer vector indicating the equation type used for each layer
#'   (confined or unconfined, see MODFLOW-2000 documentation)\cr
#' \code{$DELR}: column spacing (along rows), may be a single constant value\cr
#' \code{$DELC}: row spacing (along columns), may be a single constant value\cr
#' \code{$elev}: a vector of layer divide elevations or a 3D array of distributed
#'   layer divide elevations by cell\cr
#' \code{$sps}: a data.frame with stress period descriptions:\cr
#' \code{..$PERLEN} (num): stress period length\cr
#' \code{..$NSTP} (int): number of time steps in the stress period\cr
#' \code{..$TSMULT} (num): time step multiplier\cr
#' \code{..$TR} (log): does the stress period use transient equations?\cr
#'
#' @export
#'
#' @examples
#'
read.DIS <- function(file){
  txt <- readLines(file)
  txt <- txt[!is.comment(txt, "#")]

  #gets the numbers from the first (non-comment) entry
  extent <- as.integer(str_extract_all(txt[1], "\\d+", simplify = T))
  names(extent) <- c("NLAY", "NROW", "NCOL", "NPER", "t_unit", "l_unit")

  #gets the layer types
  LAYCBD <- c(sapply(str_extract_all(txt[2], "\\d+"), as.integer))

  #row spacing
  nDRlns <- expected.lines.RIARRAY(txt[3], extent["NCOL"], 1L, 1L)
  DELR <- interpret.RIARRAY(txt[3 + 0:nDRlns], T)

  #column spacing
  nDClns <- expected.lines.RIARRAY(txt[3L + 1L + nDRlns],
                                   extent["NROW"], 1L, 1L)
  DELC <- interpret.RIARRAY(txt[4 + nDRlns + 0:nDClns], T)

  ln <- 4L + nDRlns + nDClns

  #layer divide elevations
  ellns <- integer(extent["NLAY"] + 1L)
  elev <- rep(list(NULL), extent["NLAY"] + 1L)
  for(ld in 0:extent["NLAY"]){
    ellns[ld + 1] <- expected.lines.RIARRAY(txt[ln + 1], extent["NCOL"], extent["NROW"], 1L)
    elev[[ld + 1]] <- interpret.RIARRAY(txt[ln + 1 + 0:ellns[ld + 1]], T, extent["NCOL"])
    ln <- ln + 1L + ellns[ld + 1]
  }

  if(all(sapply(elev, length) == 1)) elev <- unname(do.call(c, elev)) else{
    varis <- vapply(elev, is.matrix, logical(1))
    elev[!varis] <- lapply(elev[!varis], matrix, extent["NCOL"], extent["NROW"])
    # note that rows and columns are being reversed in this R interface
    elev <- do.call(abind, c(elev, list(along = 3L)))
  }

  sps <- as.data.frame(scan(text = trimws(txt[ln + 1:extent["NPER"]]),
                            what = list(double(), integer(), double(), character()), quiet = T))
  names(sps) <- c("PERLEN", "NSTP", "TSMULT", "TR")
  sps$TR <- ifelse(str_to_upper(sps$TR) == "TR", T, F) #convert to logical - true if transient

  structure(mget(c("extent", "LAYCBD", "DELR", "DELC", "elev", "sps")),
            class = "DIS.MFpackage")
}

#' Write a discretisation (DIS) package file
#'
#' Writes information from a DIS.MFpackage list object to a MODFLOW-readable
#'  DIS package file.
#'
#' @param DIS object of class DIS.MFpackage, as would be read by
#'   \code{\link{read.DIS}}
#' @param filename
#' character string;
#' file name to write to, ideally ending in \code{".dis"}
#' @param title character string optional title to put at the start of the
#'   package file
#'
#' @return \code{NULL}
#'
#' @export
#'
#' @examples
#'
write.DIS <- function(DIS, filename, title){
  con <- file(filename, "wt")
  on.exit(close(con))

  if(!missing(title)){
    writeLines(paste("#", title), con)
  }

  writeLines({
    "# MODFLOW 2000 DIS package file created by write.DIS function, in R"
  }, con)

  # model dimensions
  writeLines(paste0(" ", paste(DIS$extent, collapse = "  ")), con)

  # layer types
  writeLines(paste0(" ", paste(DIS$LAYCBD, collapse = " ")), con)

  # column and row spacings
  writeLines(if(isTRUE(all.equal(names(DIS$DELR), "CNSTNT"))){
    RIARRAY(CNSTNT = DIS$DELR, FMTIN_type = "e")
  }else{
    RIARRAY(DIS$DELR, FMTIN_type = "e", flag.no = 29L)
  }, con)
  writeLines(if(isTRUE(all.equal(names(DIS$DELC), "CNSTNT"))){
    RIARRAY(CNSTNT = DIS$DELC, FMTIN_type = "e")
  }else{
    RIARRAY(DIS$DELC, FMTIN_type = "e")
  }, con)

  # elevations
  writeLines(RIARRAY.splitlayers(DIS$extent["NLAY"] + 1L, arr = DIS$elev, flag.no = 29L),
             con)

  # stress periods
  setDT(DIS$sps)
  on.exit(setDF(DIS$sps), add = TRUE)
  ff.widths <- c(10L, 4L, 10L)
  ff.digit <- c(3L, 0L, 3L)
  ff.format <- c("e", "d", "e")
  DIS$sps[, {

    ffs <- matrix("", .N, 4L)
    ffs[, 1:3] <- vapply(1:3,
                         function(col){
                           formatC(.SD[[col]], width = 10L,
                                   digits = ff.digit[col],
                                   format = ff.format[col])
                         }, character(.N))

    ffs[, 4L] <- ifelse(TR, "TR", "SS")

    writeLines(apply(ffs, 1L, paste, collapse = " "), con)
  }]

  invisible()
}
