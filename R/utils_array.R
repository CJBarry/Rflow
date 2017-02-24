# Rflow package - array reading and writing utilities

#' Number of lines in a MODFLOW text array
#'
#' Uses the header file to determine how many lines to read for a text
#'  array, excluding the header itself.
#'
#' @param header.txt
#' @param ncol
#' @param nrow
#' @param nlay
#'
#' @return
#'
#' @importFrom stringr str_extract_all
#'
#' @examples
expected.lines.RIARRAY <- function(header.txt, ncol, nrow, nlay = 1L){
  # constant value
  if(identical(as.integer(substr(header.txt[1], 1, 10)), 0L)) return(0L)

  fmtin <- substr(header.txt, 21, 40)
  fmtin <- as.numeric(str_extract_all(trimws(fmtin), "\\d+", T))

  valperline <- fmtin[1]

  # last term 0 if exact multiple, 1 otherwise
  lineperrow <- ncol %/% valperline + as.logical(ncol %% valperline)

  return(as.integer(lineperrow*nrow*nlay))
}

#' Convert lines of text to a numeric array
#'
#' @param txt.lines
#' @param header
#' @param ncol
#' @param nlay
#'
#' @return
#'
#' @importFrom stringr str_c
#'
#' @examples
interpret.RIARRAY <- function(txt.lines, header = T, ncol, nlay = 1L){
  # only a simple array or constant value coded thus far, "FREE" in the
  #  fmtin will give an error currently
  if(header){
    # if constant value
    if(identical(as.integer(substr(txt.lines[1], 1, 10)), 0L)){
      return(c(CNSTNT = as.numeric(substr(txt.lines[1], 11, 20))))
    }

    # multiplier value read
    # - note that a value of 0 is valid, but is converted to 1 by MODFLOW
    multiplier <- as.numeric(substr(txt.lines[1L], 11L, 20L))
    if(multiplier == 0) multiplier <- 1
    if(is.na(multiplier)){
      multiplier <- 1
      warning("value read for multiplier; corrected to 1")
    }

    # finds how many characters are occupied by each number (fmtin[2])
    fmtin <- substr(txt.lines[1], 21, 40)
    fmtin <- as.numeric(str_extract_all(trimws(fmtin), "\\d+", T))

    txt.merge <- str_c(txt.lines[-1], collapse = "")
    nc <- nchar(txt.merge)

    w <- fmtin[2]
    vals <- vapply(seq(1L, nc, w), function(start){
      substr(txt.merge, start, start + w - 1L)
    }, character(1))
    vals <- multiplier*
      (if(length(fmtin) == 3L) as.numeric(vals) else as.integer(vals))
  }else{
    # relies on space between each value if no header to give fixed width
    txt.merge <- str_c(txt.lines, collapse = "")

    vals <- numeric.extract(txt.merge, "")
  }

  if(!missing(ncol)){
    if((lv <- length(vals)) %% nlay*ncol == 0){
      mtx <- drop(array(vals, c(ncol, lv/(nlay*ncol), nlay)))
      return(mtx)
    }else warning("number of values is not a multiple of nlay*ncol\n")
  }

  return(vals)
}

#' Skip over unneed binary arrays
#'
#' @param con
#' file connection
#' @param nvec
#' integer \code{[]};
#' numbers of entries
#' @param bnvec
#' integer \code{[]};
#' numbers of bytes per entry
#'
#' @return
#' \code{NULL}
#'
skiparr <- function(con, nvec, bnvec){
  Map(readBin, n = nvec, size = bnvec, MoreArgs = list(con = con, what = "integer"))
  invisible()
}

#' Flip an array
#'
#' reflects an array in a given dimension
#'
#' @param arr
#' array
#' @param by
#' integer \code{[1]}
#'
#' @return
#' array
#'
#' @examples
#' mtx <- matrix(1:4, 2, 2)
#' mtx
#' flip(mtx, 1L)
#' flip(mtx, 2L)
#'
flip <- function(arr, by){
  if(is.vector(arr)) return(rev(arr))
  dims <- dim(arr)
  args <- c(list(arr), rep(list(bquote()), by - 1L),
            list(dims[by]:1), rep(list(bquote()), length(dims) - by),
            drop = FALSE)
  do.call(`[`, args)
}
