## functions that read MODFLOW data


readHDS.arr.old <- function(file, conc = FALSE, time.only = FALSE,
                        time.bn = 4L, hd.bn = 4L, show.help = TRUE,
                        nf.to.NA = FALSE, nf.val = 999,
                        CRs = list("all", "all"), lays = "all",
                        sp_ts = "all"){
  conc <- conc[1] #just in case
  if(!identical(sp_ts, "all")){
    #if sp_ts is given as a character vector (apart from "all"), then it is expected that each item will be of the form SP_TS[_TTS] where SP and TS may be converted to integers
    if(is.character(sp_ts)) sp_ts <- matrix(as.integer(do.call(rbind, str_split(sp_ts, "_"))), ncol = ifelse(conc, 3L, 2L))
    if(is.vector(sp_ts)) sp_ts <- t(sp_ts) #ensure matrix if selecting only some time steps
  }

  #get all spatial indices?
  names(CRs) <- c("C", "R")
  allspis <- vapply(CRs, identical, logical(1L), "all")
  CRsfi <- llply(CRs, function(is) if(identical(is, "all")) bquote() else is) #used for convenient indexing

  to.read <- file(file, "rb")
  i <- "integer"; d <- "double"
  tssp <- readBin(to.read, i, ifelse(conc, 3L, 2L), 4L) #[transport timestep], timestep and stress period number
  if(!conc) readBin(to.read, d, 1L, time.bn) #time in current stress period (ignored)
  times <- readBin(to.read, d, 1L, time.bn) #total elapsed time
  type <- readChar(to.read, 16L)
  CR <- structure(readBin(to.read, i, 2L, 4L), names = c("C", "R")) #number of columns and rows
  lay <- readBin(to.read, i, 1L, 4L) #layer number

  if(!identical(lays, "all") && !1L %in% lays) lay <- NA_integer_
  if(!identical(sp_ts, "all") && !any(apply(sp_ts, 1L, function(st) all(st == c(1L, 1L, if(conc) 1L))))) tssp[] <- NA_integer_

  #pre-allocating - vast speed-ups from this step (no Option Explicit in R)
  fs <- file.info(file)$size
  bpa <- prod(CR)*hd.bn + ifelse(conc, 40L, 36L) + ifelse(conc, 1L, 2L)*time.bn #bytes per array (including metadata)
  nar.est <- as.integer(fs/bpa) + 10L #small overestimate for safety

  tssp <- rbind(tssp, matrix(nrow = nar.est - 1L, ncol = ifelse(conc, 3L, 2L)))
  times <- c(times, rep(NA, nar.est - 1L))
  lay <- c(lay, rep(NA, nar.est - 1L))
  if(!time.only) hds <- array(dim = c(ifelse(allspis, CR, lengths(CRs)), nar.est))
  if(show.help) cat("model has ", CR[2], " rows, ", CR[1], " columns\n", sep = "")
  if(!time.only && (identical(lays, "all") || lay[1] %in% lays)){
    if(all(allspis)) hds[,, 1] <- readBin(to.read, d, prod(CR), hd.bn) else{
      hds[,, 1] <- do.call(`[`, c(list(array(readBin(to.read, d, prod(CR), hd.bn), CR)), CRsfi))
    }
  } else readBin(to.read, d, prod(CR), hd.bn)

  if(show.help) cat("reading arrays:\n")
  an <- 2L
  repeat{
    #time step, stress period
    tssp.new <- readBin(to.read, i, ifelse(conc, 3L, 2L), 4L)
    if(identical(tssp.new, integer(0))){if(show.help) cat("\nend of file\n"); break}
    if(!identical(sp_ts, "all") && !any(apply(sp_ts, 1L, function(st) all(st == rev(tssp.new))))){
      #if all the requested stress periods are now past, might as well skip to end (assumes time steps are in order)
      if(all(sp_ts[, 1] < tssp.new[2])){cat("\nread all requested time steps\n"); break}

      #this timestep not requested, so skip
      skiparr(to.read, c(ifelse(conc, 1L, 2L), 4L, 3L, prod(CR)), c(time.bn, 4L, 4L, hd.bn))
      next
    }
    tssp[an,] <- tssp.new

    #time and array type (latter ignored)
    if(!conc) readBin(to.read, d, 1L, time.bn) #time in current stress period (ignored)
    times[an] <- readBin(to.read, d, 1L, time.bn) #same time value for each layer: resolved later
    readChar(to.read, 16L) #assumed all the same in an hds file

    #layer
    lay.new <- readBin(to.read, i, 3L, 4L)[3L]
    if(!identical(lays, "all") && !lay.new %in% lays){
      #this layer not requested, so skip
      skiparr(to.read, prod(CR), hd.bn)
      next
    }
    lay[an] <- lay.new

    if(time.only) readBin(to.read, d, prod(CR), hd.bn) else{
      if(all(allspis)) hds[,, an] <- readBin(to.read, d, prod(CR), hd.bn) else{
        hds[,, an] <- do.call(`[`, c(list(array(readBin(to.read, d, prod(CR), hd.bn), CR)), CRsfi))
      }
    }
    an <- an + 1L
    if(show.help) cat(".")
  }

  #remove unneeded indices resulting from overestimate
  keep <- !is.na(rowSums(tssp)) & !is.na(lay)
  tssp <- unname(tssp[keep,, drop = F])
  times <- times[keep]; lay <- lay[keep]
  times <- times[seq(length(unique(lay)), length(lay), length(unique(lay)))] #only keep one time record per layer
  if(!time.only) hds <- hds[,, keep, drop = F]

  #find unique tssp and array types
  untssp <- unique.matrix(tssp, MARGIN = 1L) #this is a matrix still; the function does not "drop"
  unts <- 1:nrow(untssp) #unique timestep numbers
  ts <- vapply(1:nrow(tssp), function(r){
    unts[Reduce(`&`, lapply(1:ifelse(conc, 3, 2), function(c) untssp[, c] == tssp[r, c]))]
  }, integer(1)) #assign unique timestep numbers

  names(times) <- apply(flip(untssp, 2), 1, paste, collapse = "_")
  if(time.only){close(to.read); return(times)}

  hds <- lapply(unts, function(tsn) hds[,, ts == tsn, drop = F])
  hds <- structure(do.call(c, hds), dim = c(ifelse(allspis, CR, lengths(CRs)), L = length(unique(lay)), ts = length(unts)),
                   dimnames = c(ifelse(allspis, llply(CR, seq), CRs),
                                list(unique(lay)), list(names(times))))
  names(dimnames(hds)) <- c("C", "R", "L", ifelse(conc, "sp_ts_tts", "sp_ts"))

  #convert no-flows to NA if requested
  if(nf.to.NA) hds[hds == nf.val] <- NA

  lst <- list(hds, times)
  names(lst) <- c(str_to_title(str_replace_all(type, " ", "")), "time")
  if(show.help){
    cat(CR[1], "rows,", CR[2], "columns,", length(unique(lay)), "layers\n")
    cat("The array values are stored in a four-dimensional array in item $", names(lst)[1], ": [column, row, layer, timestep] and the time values for each timestep are stored in a vector in item $time.\n", sep = "")
    if(conc) cat("Flow Timesteps:\n") else cat("Stress periods:\n")
    cat(vapply(unique(untssp[, ifelse(conc, 3, 2)]), function(speq) paste0(speq, ": ", max(untssp[untssp[, ifelse(conc, 3, 2)] == speq, 1]), ifelse(conc, " t", " "), "ts"), character(1)), sep = c(rep("; ", 4), "\n"))
  }
  close(to.read)
  return(lst)
}
#' Read Head Save array
#'
#' @param file
#' character;
#' file name of head save array
#' @param conc
#' logical;
#' is the file in fact an unformatted concentration file (MT3D output)?
#' @param time.only
#' logical;
#' only read the model time information and discard the head data
#' @param time.bn
#' integer, 4 or 8;
#' how many bytes do the time records occupy (generally 4)?
#' @param hd.bn
#' integer, 4 or 8;
#' how many bytes to the head records occupy (generally 4)?
#' @param show.help
#' logical;
#' show what is being read and print message describing output?
#' @param nf.to.NA
#' logical;
#' convert no flow cells to NA
#' @param nf.val
#' logical;
#' if nf.to.NA, what head value indicates a no flow cell?
#' @param CRs
#' length-2 list;
#' elements are either "all" (default) or vectors of column/ row indices
#'  which are requested (not currently implemented)
#' @param lays
#' either "all" (default) or integer vector of requested layers
#' @param sp_ts
#' either "all", a length-2 integer vector or a 2-column matrix;
#' which stress period-time step references to return
#'
#'
#' @return
#' a list with one or two elements:\cr
#'   \code{$data}: num \code{[NCOL, NROW, NLAY, NTS]}; head values, only if
#'    \code{time.only = FALSE}
#'   \code{$time}: num \code{[NTS]}; time at end of each time step,
#'    relative to model start
#'
#' @import plyr
#' @import abind
#' @import stringr
#' @export
#'
#' @examples
readHDS.arr <- function(file, conc = FALSE, time.only = FALSE,
                        time.bn = 4L, hd.bn = 4L, show.help = TRUE,
                        nf.to.NA = FALSE, nf.val = 999,
                        CRs = list("all", "all"), lays = "all",
                        sp_ts = "all"){
  to.read <- file(file, "rb")
  on.exit(close(to.read))
  tssp1 <- readBin(to.read, "integer", 2L, 4L)
  # only keep total elapsed time
  times1 <- readBin(to.read, "double", 2L, time.bn)[2L]
  type1 <- nicearname(readChar(to.read, 16L))
  # number of columns and rows
  CsRs <- readBin(to.read, "integer", 2L, 4L)
  # first layer number
  lay1 <- readBin(to.read, "integer", 1L, 4L)

  hds1 <- readBin(to.read, "double", nval <- prod(CsRs), hd.bn)

  # pre-allocating - vast speed-ups from this step (no Option Explicit in R)
  fs <- file.info(file)$size
  bpa <- nval*hd.bn + 36L + 2L*time.bn # bytes per array (including metadata)
  nar.est <- as.integer(fs/bpa) + 10L # small overestimate for safety

  if(!identical(sp_ts, "all") && !is.matrix(sp_ts)) sp_ts <- t(sp_ts)

  # determine whether to include in result
  retain <- `&&`(identical(lays, "all") || lay1 %in% lays,
                 identical(sp_ts, "all") || any(sp_ts[, 1L] == tssp1[2L] &
                                                  sp_ts[, 2L] == tssp1[1L]))


  hds <- vector("list", 100L)
  hds[[1L]] <- array(NA_real_, c(nval, nar.est))
  if(retain){
    names(hds)[1L] <- type1
    hds[[type1]][, 1L] <- hds1
  }

  tssp <- matrix(NA_integer_, nar.est, 2L)
  if(retain){
    tssp[1L,] <- tssp1
    colnames(tssp) <- c("ts", "sp")
  }

  times <- double(nar.est)*NA_real_
  if(retain){
    times[1L] <- times1
    names(times)[1L] <- paste(rev(tssp1), collapse = "_")
  }

  AtsL <- matrix(NA_integer_, nar.est, 2L)
  if(retain) AtsL[1L,] <- c(1L, lay1)

  tsn <- 1L
  an <- if(retain) 1L else 0L
  lay.old <- 0L
  repeat{
    # time step, stress period
    tssp.new <- readBin(to.read, "integer", 2L, 4L)
    if(identical(tssp.new, integer(0))){if(show.help) cat("\nend of file\n"); break}
    times.new <- readBin(to.read, "double", 2L, time.bn)[2L]

    tssp.name <- paste(rev(tssp.new), collapse = "_")

    type <- nicearname(readChar(to.read, 16L))
    lay <- readBin(to.read, "integer", 3L, 4L)[3L]

    # determine whether to include in result
    retain <- `&&`(identical(lays, "all") || lay %in% lays,
                   identical(sp_ts, "all") || any(sp_ts[, 1L] == tssp.new[2L] &
                                                    sp_ts[, 2L] == tssp.new[1L]))

    if(retain && !tssp.name %in% names(times)){
      tsn <- tsn + 1L
      an <- an + 1L

      tssp[tsn,] <- tssp.new

      times[tsn] <- times.new
      names(times)[tsn] <- tssp.name
    }else if(retain && lay != lay.old) an <- an + 1L

    if(retain) AtsL[an,] <- c(tsn, lay)

    ar <- readBin(to.read, "double", nval, hd.bn)
    if(retain){
      if(nf.to.NA) ar[ar == nf.val] <- NA_real_
      if(type %in% names(hds)){
        hds[[type]][, an] <- ar
      }else{
        NtypesExisting <- sum(!is.na(names(hds)) & names(hds) != "")
        PutIn <- NtypesExisting + 1L

        names(hds)[PutIn] <- type
        hds[[type]] <- array(NA_real_, c(nval, nar.est))
        hds[[type]][, an] <- ar

      }

      lay.old <- lay
    }

    if(show.help) cat(".")
  }

  # if nothing has been retained
  if(an == 0L){
    warning("Rflow::readHDS.arr: no existing layers or timesteps requested, returning dummy output")
    return(list(data = NULL, time = double(0L)))
  }

  # remove over-allocated bits
  hds <- hds[!vapply(hds, is.null, logical(1L))]
  tssp <- tssp[seq_len(an),, drop = FALSE]
  times <- times[seq_len(an)]
  AtsL <- AtsL[seq_len(an),, drop = FALSE]
  for(i in seq_along(hds)) hds[[i]] <- hds[[i]][, seq_len(an), drop = FALSE]

  # make with dimensions [NCOL,NROW,NLAY,Nts,Ntype]
  arr <- array(NA_real_, c(CsRs,
                           length(unique(AtsL[, 2L])),
                           length(times), length(hds)),
               list(NULL, NULL,
                    paste0("L", sort(unique(AtsL[, 2L]))),
                    names(times), names(hds)))

  # fill array
  for(i in seq_along(hds)) for(j in which(!is.na(times))){
    arr[,, paste0("L", AtsL[j, 2L]), AtsL[j, 1L], i] <- hds[[i]][, j]
  }
  arr <- arr[,,, !is.na(times),, drop = FALSE]
  times <- times[!is.na(times)]

  # return list of array and times
  list(data = arr, time = times)
}



#returns a 5D array [col, row, lay, ts, arr type]
#Be more cautious with converting no flows to NA, because nfs will have value 0, within the reasonable range of active cells.  Consider instead inferring nfs from the hds array.  If a particular array type at1 does not feature within a certain timestep ts1, then [,,, ts1, at1] will be returned filled with NA.


#' Read Cell-By-Cell Budget Array
#'
#' @param file
#' character string;
#' file name to read
#' @param flux.bn
#' integer \code{[1]};
#' number of bytes occupied by a flow array element (generally 4)
#' @param show.help
#' logical \code{[1]};
#' display progress markers and help message at the end?
#' @param nf.to.NA
#' logical \code{[1]};
#' make cells for no-flow cells into NA?  (requires a head array to be
#'  given to \code{hds}; can be slow)
#' @param hds
#' list, as would be output from \code{\link{readHDS.arr}}
#' @param CRLs
#' list \code{[3]};
#' each element should either by \code{"all"} or a subset of columns (first
#'  element), rows (second element) or layers (third element) to return a
#'  sub-array
#' @param sp_ts
#' character strings \code{[]} or integer matrix \code{[, 2]};
#' either \code{"all"} or a set of <stress period _ time step>s to read
#' @param artys
#' character strings \code{[]};
#' \code{"all"} or output array types to read (e.g. "FlowRightFace",
#'  "RiverLeakage")
#'
#' @return
#' 5D numeric array \code{[NCOL, NROW, NLAY, NTS, NARTY]}.  The fourth and
#'  fifth dimensions are given appropriately formatted dimnames, such that
#'  specific output types may be subsetted by name.
#'
#' @import plyr
#' @import abind
#' @import stringr
#' @export
#'
#' @examples
readCBB.arr <- function(file, flux.bn = 4L, show.help = TRUE,
                        nf.to.NA = F, hds,
                        CRLs = list("all", "all", "all"), sp_ts = "all",
                        artys = "all"){
  if(!identical(sp_ts, "all")){
    #if sp_ts is given as a character vector (apart from "all"), then it is expected that each item will be of the form SP_TS where SP and TS may be converted to integers
    if(is.character(sp_ts))sp_ts <- matrix(as.integer(do.call(rbind, str_split(sp_ts, "_"))), ncol = 2L)
    if(is.vector(sp_ts)) sp_ts <- t(sp_ts) #ensure matrix if selecting only some time steps
  }

  #get all spatial indices?
  names(CRLs) <- c("C", "R", "L")
  allspis <- vapply(CRLs, identical, logical(1L), "all")
  CRLsfi <- llply(CRLs, function(is) if(identical(is, "all")) bquote() else is) #used for convenient indexing

  i <- "integer"; d <- "double"
  to.read <- file(file, "rb")

  #first array
  tssp <- readBin(to.read, i, 2L, 4L) #stress period; timestep
  arty <- nicearname(readChar(to.read, 16L)) #array type
  spds <- structure(abs(readBin(to.read, i, 3L, 4L)), names = c("C", "R", "L")) #spatial dimensions (c, r, l)

  #pre-allocating - vast speed-ups from this step (no Option Explicit in R)
  fs <- file.info(file)$size
  if(.Platform$OS.type == "windows" && fs/(2^20) > (memory.limit() - memory.size())/4) memory.limit(memory.size() + fs*4/(2^20)) #ensure array won't occupy more than a quarter of remaining space
  bpa <- prod(spds)*flux.bn + 36L #bytes per array
  nar.est <- as.integer(fs/bpa) + 10L #small overestimate for safety
  cbb <- array(dim = c(ifelse(allspis, spds, lengths(CRLs)), nar.est))

  #should we bother reading the first array?
  if((!identical(sp_ts, "all") && !any(apply(sp_ts, 1L, function(st) all(st == c(1L, 1L))))) ||
     (!identical(artys, "all") && !arty %in% artys)){
    tssp[] <- NA_integer_
    arty <- NA_character_
    readBin(to.read, d, prod(spds), flux.bn) #skip
  }else if(all(allspis)) cbb[,,, 1] <- readBin(to.read, d, prod(spds), flux.bn) else{
    cbb[,,, 1] <- do.call(`[`, c(list(array(readBin(to.read, d, prod(spds), flux.bn), spds)), CRLsfi))
  }

  tssp <- rbind(tssp, matrix(NA, nar.est - 1L, 2L))
  arty <- c(arty, rep(NA, nar.est - 1L))

  #read in data
  if(show.help) cat("reading arrays: ")
  an <- 2L
  repeat{
    #time step, stress period
    tssp.new <- readBin(to.read, i, 2L, 4L)
    if(identical(tssp.new, integer(0))){if(show.help) cat("\nend of file\n"); break}
    if(!identical(sp_ts, "all") && !any(apply(sp_ts, 1L, function(st) all(st == rev(tssp.new))))){
      #if all the requested stress periods are now past, might as well skip to end (assumes time steps are in order)
      if(all(sp_ts[, 1] < tssp.new[2])){cat("\nread all requested time steps\n"); break}

      #this timestep not requested, so skip
      skiparr(to.read, c(4L, 3L, prod(spds)), c(4L, 4L, flux.bn))
      next
    }
    tssp[an,] <- tssp.new

    #array type (which variable is being read?)
    arty.new <- nicearname(readChar(to.read, 16L))
    if(!identical(artys, "all") && !arty.new %in% artys){
      #this array type not requested, so skip
      skiparr(to.read, c(3L, prod(spds)), c(4L, flux.bn))
      tssp[an,] <- NA
      next
    }
    arty[an] <- arty.new

    readBin(to.read, i, 3L, 4L) #skip
    if(all(allspis)) cbb[,,, an] <- readBin(to.read, d, prod(spds), flux.bn) else{
      cbb[,,, an] <- do.call(`[`, c(list(array(readBin(to.read, d, prod(spds), flux.bn), spds)), CRLsfi))
    }
    an <- an + 1L
    if(show.help) cat(".")
  }

  if(show.help) cat("arrays successfully read and bound, now sorting\n")

  #remove unneeded indices
  tssp <- unname(tssp[keep <- !is.na(arty),, drop = F])
  cbb <- cbb[,,, keep, drop = F]
  arty <- arty[keep]

  #find unique tssp and array types
  untssp <- unique.matrix(tssp, MARGIN = 1) #this is a matrix still
  unts <- 1:nrow(untssp) #unique timestep numbers
  ts <- vapply(1:nrow(tssp), function(r) unts[untssp[, 1] == tssp[r, 1] & untssp[, 2] == tssp[r, 2]], integer(1)) #assign unique timestep numbers
  unarty <- unique(arty)

  #make table (matrix) of T/F indicating whether each array type is represented by each timestep
  #columns for array types, rows for timesteps
  repd <- vapply(seq_along(unts), function(tsp){
    vapply(seq_along(unarty), function(at){
      unts[tsp] %in% ts[arty == unarty[at]]
    }, logical(1))
  }, logical(length(unarty)))
  if(!is.vector(repd)) repd <- t(repd) else repd <- as.matrix(repd)
  dimnames(repd) <- list(unts, unarty)

  #rearrange by array type
  #1. add in NA arrays where necessary, whilst updating the array type and timestep vectors/ matrices
  if(!all(repd)){
    pos <- function(tsp) length(ts[ts <= tsp]) + 1
    lapply(unarty, function(at){
      toins <- as.integer(dimnames(repd[!repd[, at], at, drop = F])[[1]])
      lapply(toins, function(tsp){
        arty <<- c(arty[1:pos(tsp)], at, arty[(pos(tsp) + 1):length(arty)])
        cbb <<- abind(cbb[,,, 1:pos(tsp), drop = F], array(NA, c(ifelse(allspis, spds, lengths(CRLs)), 1)), cbb[,,, (pos(tsp) + 1):dim(cbb)[4], drop = F], along = 4)
        tssp <<- rbind(tssp[1:pos(tsp),], untssp[unts == tsp,], tssp[(pos(tsp) + 1):length(ts),])
        ts <<- c(ts[1:pos(tsp)], tsp, ts[(pos(tsp) + 1):length(ts)])
      })
    })
  }

  #2. rearrange array by array type and promote to 5D
  #array type sequence
  atseq <- cbind.data.frame(pos = 1:length(arty), at = vapply(arty, function(at) which(unarty == at), integer(1)))
  atseq <- atseq[with(atseq, order(at)),]

  cbb <- array(cbb[,,, atseq$pos], c(ifelse(allspis, spds, lengths(CRLs)), length(unts), length(unarty)))

  #3. label time and array type dimensions, as well as giving fixed numbers to CRL dimensions
  dimnames(cbb) <- c(ifelse(allspis, llply(spds, seq), CRLs), list(apply(flip(untssp, 2), 1, paste, collapse = "_")), list(unarty))
  names(dimnames(cbb)) <- c("C", "R", "L", "sp_ts", "par")

  #convert no flows to NA if requested, using a hds array for which this has already been done
  if(nf.to.NA){
    if(!missing(hds)){
      if(is.list(hds)) hds <- do.call(`[`, c(hds["Head"], CRLsfi, list(if(sp_ts == "all") bquote() else sp_ts)))
      if(!any(is.na(hds))) cat("no no-flows detected in hds array\n") else{
        for(at in seq_len(dim(cbb)[5])) cbb[,,,, at][is.na(hds)] <- NA
      }
    }else cat("no hds array given, so nf conversion not performed\n")
  }

  #return values
  close(to.read)
  if(show.help) cat({
    "The returned object is a 5-dimensional array, where the dimensions represent [column, row, layer, timestep, array type].  The fifth dimension is named according to the saved value type: "
  }, str_c(unarty, collapse = ", "), ".\n", sep = "")
  if(show.help){
    cat("Stress periods:\n")
    cat(vapply(unique(untssp[, 2]), function(sp) paste0(sp, ": ", max(untssp[untssp[, 2] == sp, 1]), " ts"), character(1)), sep = c(rep(";", 4), "\n"))
  }
  return(cbb)
}

#' Read Unformatted Concentration File
#'
#' @param file
#' character string
#' @param time.only
#' logical \code{[1]};
#' if \code{TRUE}, then only return the time at the end of each time step
#' @param time.bn
#' integer \code{[1]};
#' number of bytes in a time record (generally 4)
#' @param conc.bn
#' integer \code{[1]};
#' number of bytes in an array record (generally 4)
#' @param show.help
#' logical \code{[1]};
#' show progress markers and explanatory notes?
#' @param nf.to.NA
#' logical \code{[1]};
#' convert inactive cells to NA?
#' @param nf.val
#' numeric \code{[1]};
#' value denoting inactive cells (commonly -1)
#' @param lays
#' \code{"all"} or integer \code{[]};
#' return all layers, or else a subset
#' @param sp_ts_tts
#' \code{"all"}, character strings \code{[]} or integer matrix \code{[, 3]};
#' return all transport time steps, or a subset of
#'  <stress period_time step_transport time step>s
#'
#' @return
#' either a list with two elements (with \code{time.only = FALSE}):\cr
#'   \code{$Concentration}: num \code{[NCOL, NROW, NLAY, NTTS]}\cr
#'   \code{$time}: num \code{[NTTS]}; time at end of transport time steps,
#'    relative to model start\cr
#' or just the time values
#'
#' @import plyr
#' @import abind
#' @import stringr
#' @export
#'
#' @examples
readUCN.arr <- function(file, time.only = F, time.bn = 4L, conc.bn = 4L, show.help = T, nf.to.NA = T, nf.val = -1, lays = "all", sp_ts_tts = "all"){
  fs <- file.info(file)$size

  if(.Platform$OS.type == "windows" && fs/(2^20) > (memory.limit() - memory.size())/4) memory.limit(memory.size() + fs*4/(2^20)) #ensure array won't occupy more than a quarter of remaining space

  ucn <- readHDS.arr(file, conc = T, time.only = time.only, time.bn = time.bn, hd.bn = conc.bn, show.help = show.help, nf.to.NA = nf.to.NA, nf.val = nf.val, lays = lays, sp_ts = sp_ts_tts)

  return(ucn)
}
