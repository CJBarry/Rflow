# Rflow package - array reading and writing utilities


# reading arrays ----------------------------------------------------------

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


# writing arrays ----------------------------------------------------------

#' Formatted array
#'
#' @param arr
#' @param df
#' @param CNSTNT
#' @param FMTIN_type
#' @param FMTIN_w
#' @param FMTIN_d
#' @param FMTIN_N
#' @param IPRN
#' @param header
#' @param nrowmod
#' @param ncolmod
#' @param bgC
#' @param flag.no
#'
#' @return
#' character string
#'
#' @examples
#' Rflow:::RIARRAY(matrix(1:4, 2, 2), FMTIN_type = "I", FMTIN_w = 3L)
#' Rflow:::RIARRAY(matrix(1:4, 2, 2), FMTIN_type = "e", FMTIN_w = 12L)
#' Rflow:::RIARRAY(CNSTNT = 1L, FMTIN_type = "I", FMTIN_w = 3L)
#'
RIARRAY <- function(arr = NULL, df = NULL, CNSTNT = 0,
                    FMTIN_type = "e", FMTIN_w = 12, FMTIN_d = 4, FMTIN_N = ncol(arr),
                    IPRN = -1, header = T, nrowmod, ncolmod, bgC = 0, flag.no = NA){

  #see if arr was intended as CNSTNT
  if(identical(names(arr), "CNSTNT")){CNSTNT <- arr; arr <- NULL}

  #infer reading in style (array, blocks, zoned array, constant respectively)
  IREAD <- if(!is.null(arr) && is.null(df)) 100L else
    if(is.null(arr) && !is.null(df)) 101L else
      if(!is.null(arr) && !is.null(df)) 102L else 0L

  if(identical(IREAD, 0L) & !header){warning("constant value so header reset to TRUE\n"); header <- T}
  if(IREAD == 0L) flag.no <- 0L

  if(is.vector(arr)) arr <- as.array(arr) #ensure array

  #switch to row-major order and collapse to matrix
  #needs dummy array if constant
  if(!is.null(arr)) arr <- matrix(arr, prod(dim(arr)[-1]), dim(arr)[1], byrow = T) else arr <- t(NA)

  #allocate character vector
  txt <- if(identical(IREAD, 0L)) character(1) else if(identical(IREAD, 100L)){
    character(prod(dim(arr)[-2]) + header)
  }else if(identical(IREAD, 101L)){
    character(nrow(df) + 2L + header)
  }else character(prod(dim(arr)[-1]) + 2L + header)

  #header
  if(header){
    if(identical(IREAD, 102L)) FMTIN_type <- "I"
    incd <- ifelse(identical(str_to_upper(FMTIN_type), "I"), F, T)
    FMTIN <- if(IREAD %in% c(0L, 100L, 102L)){
      paste0("(", FMTIN_N, str_to_upper(FMTIN_type),
             FMTIN_w, ifelse(incd, paste0(".", FMTIN_d), ""), ")")
    }else ""
    txt[1] <- paste0(FFI(ifelse(is.na(flag.no), IREAD, flag.no), 10), FFgen(CNSTNT, FMTIN_type, 10, 2), FMTIN, str_dup(" ", 20L - nchar(FMTIN)), FFI(IPRN, 10))
  }

  if(identical(IREAD, 0L)) return(txt)

  if(identical(IREAD, 100L)){ #distributed array
    #array
    for(ln in 1:nrow(arr)){
      txt[ln + as.integer(header)] <-
        cc(paste0(vapply(arr[ln,], FFgen, character(1L), FMTIN_type, FMTIN_w, FMTIN_d),
                  c(rep("", min(ncol(arr), FMTIN_N) - 1L), "\n")))
    }
  }

  if(identical(IREAD, 101L)){ #blocks
    if(missing(nrowmod) || missing(ncolmod)) stop("RIARRAY: block format (IREAD = 101) requires nrowmod and ncolmod to be specified for background block.\n")
    if(identical(ncol(df), 3L)) df <- df[, c(1, 1, 2, 2, 3)] #if blocks are specified as single cells
    txt[1 + header] <- as.character(nrow(df) + 1L) #NBLOCK
    txt[2 + header] <- str_c(1, nrowmod, 1, ncolmod, bgC, sep = " ") #background block
    for(ln in 1:nrow(df)){
      txt[ln + 2L + as.integer(header)] <- str_c(df[ln,], collapse = " ")
    }
  }

  if(identical(IREAD, 102L)){ #zone array
    #target is vector of zone values in order - would normally be given as a vector in the first place
    if(is.data.frame(df)){
      df <- if(identical(ncol(df), 1L)) df[, 1] else{
        #assumes that first column is zone number (if data type is integer) and last column is zone value, then reorders
        if(is.integer(df[, 1])) df[with(df, order(df[, 1])), ncol(df)] else df[, ncol(df)]
      }
    }
    if(!is.vector(df)) stop("RIARRAY: zone array format (IREAD = 102) requires vector or data frame input to df.\n")
    txt[1 + header] <- length(df)
    txt[2 + header] <- str_c(df, collapse = " ")
    for(ln in 1:nrow(arr)){
      txt[ln + 2L + as.integer(header)] <-
        cc(paste0(vapply(arr[ln,], FFI, character(1L), FMTIN_w),
                  c(rep("", min(ncol(arr), FMTIN_N)- 1L), "\n")))
    }
  }

  # ensure no double line-breaks
  str_c(rmnewline(txt), collapse = "\n")
}

# usually it is required that distributed parameters are given with one array per layer
# here is a wrapper for enabling this
# compress = TRUE means that layers with constant values are represented as CNSTNT
# if collapse = TRUE, paste elements together into single string

#' Formatted array, split by layer
#'
#' @param NLAY
#' @param CNSTNT
#' @param arr
#' @param compress
#' @param collapse
#' @param ...
#'
#' @importFrom abind abind
#'
#' @return
#' character string
#'
RIARRAY.splitlayers <- function(NLAY, CNSTNT = rep(0, NLAY), arr = NULL,
                                compress = TRUE, collapse = TRUE, ...){
  # CNSTNT can be given as a single value to be recycled
  CNSTNT <- expand.vec(CNSTNT, NLAY)

  # allows an array to be given without naming the argument
  if(is.array(CNSTNT)) arr <- CNSTNT

  # promote arr to three dimensions if necessary
  if(length(dim(arr)) == 2L) arr <- abind(arr, along = 3L)

  txt <- if(is.null(arr)){
    vapply(1:NLAY, function(l) RIARRAY(CNSTNT = CNSTNT[l], ...), character(1L))
  }else{
    vapply(1:NLAY, function(l) if(compress && diff(range(arr[,, l])) == 0){
      # all elements of this layer were equal, so write as constant value
      RIARRAY(CNSTNT = arr[1L, 1L, l], ...)
    }else RIARRAY(arr = arr[,, l], ...), character(1L))
  }

  if(collapse) paste(rmnewline(txt), collapse = "\n") else rmnewline(txt)
}

# easy concatenate
# will paste vectors as well as individual arguments into a single string

#' Easy concatenate
#'
#' @param ...
#'
#' @return
#' character string
#'
#' @importFrom plyr splat
#'
cc <- function(...){
  splat(paste0)(c(...))
}

#' Easy concatenate
#'
#' @param ...
#'
#' @return
#' character string
#'
#' @importFrom plyr splat
#'
cc_ <- function(...){
  splat(paste)(c(...))
}

#' Remove a tailing line end
#'
#' @param s
#' vector of character strings
#'
#' @return
#' vector of character strings with any \code{"\\n"}s removed from the ends
#'
rmnewline <- function(s){
  ifelse(substring(s, nchar(s)) == "\n",
         substring(s, 1L, nchar(s) - 1L), s)
}


# binary arrays -----------------------------------------------------------

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


# manipulating arrays -----------------------------------------------------

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
#' Rflow:::flip(mtx, 1L)
#' Rflow:::flip(mtx, 2L)
#'
flip <- function(arr, by){
  if(is.vector(arr)) return(rev(arr))
  dims <- dim(arr)
  args <- c(list(arr), rep(list(bquote()), by - 1L),
            list(dims[by]:1), rep(list(bquote()), length(dims) - by),
            drop = FALSE)
  do.call(`[`, args)
}
