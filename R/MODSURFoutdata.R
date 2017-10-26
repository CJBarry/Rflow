## functions that read MODFLOW-SURFACT data


#' Read Fracture Well 5 Flow output
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


#' Read Fracture Well 5 Head and Concentration output
#'
#' @param file
#' character string;
#' file name to read
#' @param well.name
#' character string vector;
#' names of wells in correct order, if you know them
#' @param hd.bn 
#' integer, 4 or 8; 
#' how many bytes to the head records occupy (generally 4)? 
#' @param time.bn 
#' integer, 4 or 8; 
#' how many bytes do the time records occupy (generally 4)? 
#' @param show.help
#' logical;
#' show progress markers? 
#'
#' @return
#' a list with one or two elements:\cr 
#'   \code{$hw5}: num \code{[NWell, NTS, NType]}; head values
#'   \code{$time}: num \code{[NTS]}; time at end of each time step, 
#'    relative to model start 
#'    
#' Note that \code{NTS} is the number of output time steps, rather than
#'  necessarily the number of modelled time steps.
#'
#' @export
#'
read.HW5 <- function(file, well.names = NULL, hd.bn = 4L, time.bn = 4L,
                     show.help = TRUE){
  to.read <- file(file, "rb")
  on.exit(close(to.read))
  tssp1 <- readBin(to.read, "integer", 2L, 4L)
  # only keep total elapsed time
  times1 <- readBin(to.read, "double", 2L, time.bn)[2L]
  type1 <- nicearname(readChar(to.read, 16L))
  # number of wells and two 1s
  W11 <- readBin(to.read, "integer", 3L, 4L)
  
  nwell <- W11[1L]
  hw5.1 <- readBin(to.read, "double", nwell, hd.bn)
  
  # pre-allocating - vast speed-ups from this step (no Option Explicit in R) 
  fs <- file.info(file)$size 
  bpa <- prod(W11)*hd.bn + 36L + 2L*time.bn # bytes per array (including metadata) 
  nar.est <- as.integer(fs/bpa) + 10L # small overestimate for safety 
  
  hw5 <- vector("list", 100L)
  hw5[[1L]] <- array(NA_real_, c(nwell, nar.est))
  names(hw5)[1L] <- type1
  hw5[[type1]][, 1L] <- hw5.1
  
  tssp <- matrix(NA_integer_, nar.est, 2L)
  tssp[1L,] <- tssp1
  colnames(tssp) <- c("ts", "sp")
  
  times <- c(times1, double(nar.est - 1L)*NA_real_)
  names(times)[1L] <- paste(rev(tssp1), collapse = "_")
  
  an <- 1L
  repeat{
    # time step, stress period
    tssp.new <- readBin(to.read, "integer", 2L, 4L)
    if(identical(tssp.new, integer(0))){if(show.help) cat("\nend of file\n"); break}
    times.new <- readBin(to.read, "double", 2L, time.bn)[2L]
    
    tssp.name <- paste(rev(tssp.new), collapse = "_")
    if(!tssp.name %in% names(times)){
      an <- an + 1L
      
      tssp[an,] <- tssp.new
      
      times[an] <- times.new
      names(times)[an] <- tssp.name
    }
    
    type <- nicearname(readChar(to.read, 16L))
    readBin(to.read, "integer", 3L, 4L)
    if(type %in% names(hw5)){
      hw5[[type]][, an] <- readBin(to.read, "double", nwell, hd.bn)
    }else{
      NtypesExisting <- sum(!is.na(names(hw5)) & names(hw5) != "")
      PutIn <- NtypesExisting + 1L
      
      names(hw5)[PutIn] <- type
      hw5[[type]] <- array(NA_real_, c(nwell, nar.est))
      hw5[[type]][, an] <- readBin(to.read, "double", nwell, hd.bn)
    }
    
    if(show.help) cat(".")
  }
  
  # remove over-allocated bits
  hw5 <- hw5[!vapply(hw5, is.null, logical(1L))]
  tssp <- tssp[seq_len(an),, drop = FALSE]
  times <- times[seq_len(an)]
  for(i in seq_along(hw5)) hw5[[i]] <- hw5[[i]][, seq_len(an), drop = FALSE]
  
  # collapse hw5 into single array
  hw5 <- structure(c(hw5, recursive = TRUE),
                   dim = c(Nwell = nwell, Nts = an, Ntype = length(hw5)),
                   dimnames = list(well = well.names,
                                   sp_ts = names(times),
                                   type = names(hw5)))
  
  # return list of array and times
  list(hw5 = hw5, time = times)
}
