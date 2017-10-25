## functions that read MODFLOW-SURFACT data


#' Read Fracture Well 5 output
#' 
#' A wrapper for readCBB.arr that correctly formats the result for the
#'  fracture well 5 output from MODFLOW-SURFACT.
#'
#' @param file
#' character string;
#' file name to read
#' @param well.name
#' character string vector;
#' names of wells in correct order, if you know them
#' @param flux.bn
#' integer \code{[1]};
#' number of bytes occupied by a flow array element (generally 4)
#' @param show.help
#' logical \code{[1]};
#' display progress markers and help message at the end?
#' @param sp_ts
#' character strings \code{[]} or integer matrix \code{[, 2]};
#' either \code{"all"} or a set of <stress period _ time step>s to read
#' @param artys
#' character strings \code{[]};
#' \code{"all"} or output array types to read (e.g. "FlowRightFace",
#'  "RiverLeakage")
#'
#' @return
#' 4D numeric array \code{[NLAY, NWells, NTS, NARTY]}.  The third and
#'  fourth dimensions are given appropriately formatted dimnames, such that
#'  specific output types may be subsetted by name.
#'
#' @import abind
#' @export
#'
#' @examples
read.CW5 <- function(file, well.names = NULL, flux.bn = 4L,
                     sp_ts = "all", artys = "all", show.help = TRUE){
  bud <- readCBB.arr(file, flux.bn, show.help, sp_ts = sp_ts, artys = artys)
  
  # different dimensions mean different things in CW5 output files, so
  #  adjust accordingly
  bud <- adrop(bud, c(FALSE, FALSE, TRUE, FALSE, FALSE))
  dnm <- dimnames(bud)
  names(dim(bud)) <- c("L", "well", "sp_ts", "")
  dimnames(bud)[[1L]] <- paste0("L", seq_len(dim(bud)[1L]))
  dimnames(bud)[2L] <- list(well.names)
  dimnames(bud)[3:4] <- dnm[3:4]
  
  bud
}
