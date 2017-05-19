# which.nc

#' "which" for NetCDF arrays
#'
#' @description
#' A version of "which" for NetCDF datasets that does not require prior
#' loading of arrays and avoids loading the whole of large arrays at once.
#'
#' @param ncfile NetCDF connection object
#' @param variable variable name within NetCDF data set
#' @param FUN function to operate on array, which should produce a logical
#' vector or array
#' @param arr.ind,useNames passed to which
#' @param size.threshold how large can the array be before the function
#' performs on chunks?
#' @param ... additional arguments to FUN
#'
#' @return
#' integer vector (arr.ind = FALSE) or matrix (arr.ind = TRUE), as with
#' `which`
#'
#' @import RNetCDF
#' @export
#'
#' @examples
#' library(RNetCDF)
#'
#' mfdata <- open.nc(system.file("rflow_mf_demo.nc", package = "Rflow"))
#'
#' # find which cells have a well in, by C, R, L reference
#' unique(which.nc(mfdata, "Wells", `!=`, 0, arr.ind = TRUE)[, 1:3],
#'        MARGIN = 1L)
#'
which.nc <- function(ncfile, variable, FUN = Negate(is.na), ...,
                     arr.ind = FALSE, useNames = TRUE,
                     size.threshold = 1e7){
  # get dimensions of array
  ards <- sapply(lapply(var.inq.nc(ncfile, variable)$dimids,
                        dim.inq.nc, ncfile = ncfile),
                 `[[`, "length")
  nards <- length(ards)
  arsize <- prod(ards)
  FUN <- match.fun(FUN)

  if(arsize > size.threshold){
    # large array, perform in chunks
    # - number of chunks
    nchunk <- arsize%/%size.threshold + 1L
    #
    # - higest dimension
    bigd <- ards[length(ards)]
    #
    # - max. index range for highest dimension in a chunk
    chunkd <- ceiling(bigd/nchunk)

    # define chunk indices - roughly equal sections along highest index
    chunks <- lapply(1:nchunk, function(i){
      (1:bigd)[((1:bigd)%/%chunkd + 1L) == i]
    })
    chunks <- chunks[lengths(chunks) != 0L]

    sections <- lapply(chunks, function(ch){
      x <- var.get.nc(ncfile, variable, c(rep(NA, nards - 1L), ch[1L]),
                      c(rep(NA, nards - 1L), ch[length(ch)] - ch[1L] + 1L),
                      collapse = FALSE)
      wh <- which(FUN(x, ...), arr.ind, useNames)
      if(arr.ind){
        wh[, nards] <- wh[, nards] + ch[1L] - 1L
      }else wh <- wh + prod(ards[-nards])*(ch[1L] - 1L)
      wh
    })

    if(arr.ind){
      do.call(rbind, sections)
    }else c(sections, recursive = TRUE)
  }else{
    # small array, perform all together
    x <- var.get.nc(ncfile, variable, collapse = FALSE)
    which(FUN(x, ...), arr.ind, useNames)
  }
}
