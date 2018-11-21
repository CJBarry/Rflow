## functions that read MODFLOW data


#' Read Head Save array.
#'
#' Compatible with MODFLOW-SURFACT output (.HDS and .CON) which may have
#'  more than one output type.  Also compatible with MODHMS outputs which
#'  may have outputs with different dimensions (channel segment heads).
#'
#' @param file
#' character;
#' file name of head save array
#' @param conc
#' logical;
#' is the file in fact an unformatted concentration file (MT3D .UCN output)?
#'  (Use FALSE for MODFLOW-SURFACT .CON output, for which the format is the
#'  same as .HDS)
#' @param flux
#' logical;
#' is the file in fact a budget file (.CBB or similar output)?  This can
#'  fail if there are a different number of array types for each time step
#'  (e.g. if the first time step is transient, then there will be no
#'  storage array for that time step).  In this case, use readCBB.arr
#'  instead.
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
#' either "all" (default) or integer vector of requested layers (not
#'  currently implemented)
#' @param sp_ts
#' either "all", a 2-column integer matrix (columns for stress period, time
#'  step), a character vector (elements in the form "sp_ts") or an integer
#'  vector giving global time step numbers
#' @param artypes
#' character [];
#' names of array types to be read (e.g. "Head", "FlowRightFace"); either
#'  in capitals (" FLOW RIGHT FACE") or 'nice' format ("FlowRightFace");
#'  "all" (default) returns all array types
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
readHDS.arr <- function(file, conc = FALSE, flux = FALSE, time.only = FALSE,
                        time.bn = 4L, hd.bn = 4L, show.help = TRUE,
                        nf.to.NA = FALSE, nf.val = 999,
                        dry.to.NA = nf.to.NA, dry.val = nf.val,
                        CRs = list("all", "all"), lays = "all",
                        sp_ts = "all", artypes = "all"){
  if(flux && time.only) stop("readHDS.arr: flux and time.only set to TRUE, but budget files do not store time.")

  if(file.info(file)$size == 0) stop(sprintf("readHDS.arr: %s is empty", file))

  to.read <- file(file, "rb")
  on.exit(close(to.read))

  # estimate number of arrays
  #
  # - time step, stress period of first array
  tssp <- tssp1 <- readBin(to.read, "integer", 2L + conc, 4L)
  #
  if(!flux) readBin(to.read, "numeric", 2L - conc, time.bn)
  readChar(to.read, 16L)
  dimset1 <- readBin(to.read, "integer", 2L + flux, 4L)
  l1 <- if(!flux) readBin(to.read, "integer", 1L, 4L) else 1L
  readBin(to.read, "numeric", prod(dimset1), hd.bn)
  #
  # - find dimension sets of first group of arrays (until next stress
  #    period is found or end of file is reached)
  dimsets <- list(dimset1)
  ls <- l1
  repeat{
    tssp <- readBin(to.read, "integer", 2L + conc, 4L)
    if(!length(tssp) || !all(tssp == tssp1)) break

    if(!flux) readBin(to.read, "numeric", 2L - conc, time.bn)
    readChar(to.read, 16L)

    dimset <- readBin(to.read, "integer", 2L + flux, 4L)
    dimsets <- c(dimsets, list(dimset))
    if(!flux){
      l <- readBin(to.read, "integer", 1L, 4L)
      ls <- c(ls, l)
    }

    readBin(to.read, "numeric", prod(dimset), hd.bn)
  }

  # dimension set organise
  # - unique dimension sets
  UDS <- unique(dimsets)
  NUDS <- length(UDS)
  #
  # - dimension set integer labels
  iDS <- sapply(dimsets, function(ds) which(sapply(UDS, identical, ds)))
  UiDS <- unique(iDS)
  #
  UDSCounts <- sapply(UiDS, function(i) sum(i == iDS))

  # initialise
  # - estimated maximum number of nodes
  #  -- for HDS format, this is nodes within a layer
  Nnodes_max <- max(sapply(dimsets, prod))
  Nnodes_min <- min(sapply(dimsets, prod))
  #
  #  -- reasonable safe number of nodes, given that some time steps may have
  #      different sets of arrays (e.g. storage)
  Nnodes_sensible <- (sum(sapply(dimsets, prod)) + Nnodes_min*2L)/(length(dimsets) + 2L)
  #
  # - number of layers (assume all layers were registered in the first
  #    time step)
  Nlay <- if(flux) 1L else max(ls)
  #
  # - number of array types - allow for four extra which may not have
  #    featured in the first time step (e.g. Storage)
  Narty.est <- if(identical(artypes, "all")){
    if(flux) length(dimsets) + 4L else{
      # assumes that all array types feature in layer 1
      # - OLF head (e.g.) only features in layer 1
      length(dimsets[ls == 1L]) + 4L
    }
  }else length(artypes)
  if(artypes != "all") artypes <- nicearname(artypes)
  #
  # - number of time steps - estimated based on file size and number of
  #    nodes within each array in the first time step
  nbpts_sensible <- length(dimsets)*((2L + conc)*4L +
                                       (if(flux) 0L else (2L - conc)*time.bn) +
                                       16L + 3L*4L) + sum(sapply(dimsets, prod))*hd.bn
  nts.est <- if(identical(sp_ts, "all")){
    ceiling(file.info(file)$size/nbpts_sensible*1.2)
  }else switch(class(sp_ts)[1L],
               matrix = nrow(sp_ts),
               character = length(sp_ts),
               integer = length(sp_ts),
               numeric = length(sp_ts),
               stop("readHDS.arr: invalid `sp_ts`"))
  if(is.character(sp_ts) && sp_ts != "all"){
    sp_ts <- do.call(rbind, strsplit(sp_ts, "_"))
    mode(sp_ts) <- "integer"
  }
  #
  # - stress period, timestep matrix
  #  -- include space for a global time step counter
  spts_mtx <- matrix(NA_integer_, nts.est, 3L + conc, dimnames = list(NULL, c("sp", "ts", if(conc) "tts", "ts_global")))
  #
  # - model time
  if(!flux) time <- numeric(nts.est)
  #
  # - dimensions log
  #  -- rows for dimension set identifier, then dim1 (usually col, perhaps
  #      NSeg), dim2 (usually row), dim3 (layer, in the case of budget
  #      file)
  dims <- array(NA_integer_, c(3L + flux, Narty.est))
  #
  #  -- unique dimension sets
  Udims <- list()
  #
  # - layer number log
  if(!flux) lays <- array(NA_integer_, c(Nlay, nts.est, Narty.est))
  #
  # - array type log
  artys <- array(NA_character_, c(if(!flux) Nlay else 1L, nts.est, Narty.est))
  Uartys <- rep(NA_character_, Narty.est)
  #
  # - values
  if(!time.only) vals <- array(NA_real_, c(Nnodes_max, if(!flux) Nlay else 1L, nts.est, Narty.est))

  # reset reading position
  close(to.read); to.read <- file(file, "rb")

  # get values
  tssp_prev <- tssp_prevRead <- c(0L, 0L, if(conc) 0L)
  ts_global <- 0L
  ts_read <- 0L
  force(nf.to.NA); force(nf.val)
  repeat{
    # time step, stress period
    tssp <- readBin(to.read, "integer", 2L + conc, 4L)
    if(!length(tssp)) break # end of file reached
    if(!all(tssp == tssp_prev)) ts_global <- ts_global + 1L

    readAr <- if(identical(sp_ts, "all")) TRUE else
      switch(class(sp_ts)[1L],
             matrix = if(conc){
               any(tssp[3L] == sp_ts[, 1L] & tssp[2L] == sp_ts[, 2L] & tssp[1L] == sp_ts[, 3L])
             }else{
               any(tssp[2L] == sp_ts[, 1L] & tssp[1L] == sp_ts[, 2L])
             },
             integer = ts_global %in% sp_ts,
             numeric = ts_global %in% sp_ts)

    # model time
    if(!flux) t <- readBin(to.read, "numeric", 2L - conc, time.bn)[2L - conc]

    # array type
    at <- nicearname(readChar(to.read, 16L))

    readAr <- readAr && (artypes[1L] == "all" ||
                           str_to_upper(nicearname(at)) %in% str_to_upper(nicearname(artypes)))

    if(readAr){
      atn <- which(Uartys == at)
      if(!any(atn)){
        Uartys[which(is.na(Uartys))[1L]] <- at
        atn <- which(Uartys == at)
      }
    }

    if(readAr && !all(tssp == tssp_prevRead)){
        ts_read <- ts_read + 1L
        spts_mtx[ts_read,] <- c(rev(tssp), ts_global)
      }

    if(!all(tssp == tssp_prev)){
      # sometimes MODHMS datasets have extraneous extra time steps
      if(ts_read > nts.est) break

      if(readAr && !flux) time[ts_read] <- t
      an <- 1L
    }else an <- an + 1L

    # dimensions
    spds <- readBin(to.read, "integer", 2L + flux, 4L)
    l <- if(!flux) readBin(to.read, "integer", 1L, 4L) else 1L

    if(prod(spds) != 0){
      if(readAr && !time.only){
        v <- readBin(to.read, "numeric", prod(spds), hd.bn)
        if(nf.to.NA) v[v == nf.val] <- NA_real_
        if(dry.to.NA && dry.val != nf.val) v[v == dry.val] <- NA_real_
      }else readBin(to.read, "numeric", prod(spds), hd.bn)


      if(readAr){
        dims[-1L, atn] <- spds
        if(!length(Udims) || !any(Idim <- which(sapply(Udims, identical, spds)))){
          Idim <- length(Udims) + 1L
          Udims[[Idim]] <- spds
        }
        dims[1L, atn] <- Idim

        if(!flux) lays[l, ts_read, atn] <- l

        artys[l, ts_read, atn] <- at

        if(!time.only) vals[seq_along(v), l, ts_read, atn] <- v
      }
    }

    # store ts, sp for this array
    tssp_prev <- tssp
    if(readAr) tssp_prevRead <- tssp

    if(show.help && readAr) cat(".")
  }; if(show.help) cat("\n")

  spts_mtx <- spts_mtx[seq_len(ts_read),, drop = FALSE]
  spts_labels <- apply(spts_mtx[, c("sp", "ts", if(conc) "tts"), drop = FALSE], 1L, paste, collapse = "_")
  if(!flux){
    time <- time[seq_len(ts_read)]
    names(time) <- spts_labels
    if(time.only) return(time)
  }
  if(ts_read == 0L) return(if(flux) list(data = NULL) else list(data = NULL, time = time))

  # sorting into arrays
  OutList <- vector("list", length(Udims))
  names(OutList) <- if(length(Udims) > 1L) paste0("data", seq_along(Udims)) else "data"
  #
  # - by dimension set
  for(ds in seq_along(Udims)){
    OutList[[ds]] <- vals[seq_len(prod(Udims[[ds]])),, seq_len(ts_read),
                          ds == dims[1L,] & !is.na(dims[1L,]), drop = FALSE]
    dim(OutList[[ds]]) <- if(flux){
      c(Udims[[ds]], dim(OutList[[ds]])[3:4])
    }else{
      c(Udims[[ds]], dim(OutList[[ds]])[2:4])
    }

    dimnames(OutList[[ds]])[4:5] <- list(spts_labels, Uartys[ds == dims[1L,] & !is.na(dims[1L,])])

    # remove extraneous layers and array types for this dimension set
    OutList[[ds]] <- OutList[[ds]][,,
                                   if(flux) bquote() else !apply(is.na(OutList[[ds]]), 3L, all),,
                                   !apply(is.na(OutList[[ds]]), 5L, all), drop = FALSE]
  }
  if(!flux) OutList$time <- time

  # #
  # #  -- number of layers applicable to each dimension set
  # if(!flux) UDSnlay <- sapply(UiDS, function(i) max(lays[[1L]][i == iDS]))
  # #
  # #  -- array types applicable to each dimension set
  # UDSarty <- lapply(UiDS, function(i) unique(artys[[1L]][i == iDS]))
  # UDSnarty <- lengths(UDSarty)
  # #
  # # - actual number of time steps
  # nts <- sum(!is.na(spts_mtx[, "ts_global"]))
  # keep <- seq_len(nts)
  # #
  # # - trim outputs to actual time steps
  # spts_mtx <- spts_mtx[keep,, drop = FALSE]
  # if(!flux) time <- time[keep]
  # artys <- artys[keep]
  # if(!flux) lays <- lays[keep]
  # vals <- vals[keep]
  # #
  # # - stress period, timestep labels
  # spts_labels <- apply(spts_mtx[, c("sp", "ts", if(conc) "tts"), drop = FALSE], 1L, paste, collapse = "_")
  # if(!flux) names(time) <- spts_labels
  # #
  # # - results list initialise
  # OutList <- c(Map(function(dimset_, idimset_) if(!time.only){
  #   array(c(lapply(vals, `[`, iDS == idimset_), recursive = TRUE),
  #         dim = c(dimset_, if(!flux) UDSnlay[idimset_], nts, UDSnarty[idimset_]),
  #         dimnames = list(NULL, NULL, NULL, spts_labels, UDSarty[[idimset_]]))
  # }, UDS, UiDS),
  # if(!flux) list(time = time))
  # if(!time.only){
  #   names(OutList)[1:NUDS] <- if(NUDS == 1L) "data" else paste0("data", 1:NUDS)
  # }else OutList <- OutList$time

  OutList
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
#' @param CRLs
#' list \code{[3]};
#' each element should either by \code{"all"} or a subset of columns (first
#'  element), rows (second element) or layers (third element) to return a
#'  sub-array (not currently implemented)
#' @param sp_ts
#' character strings \code{[]} or integer matrix \code{[, 2]};
#' either \code{"all"} or a set of <stress period _ time step>s to read
#'  (not currently implemented)
#' @param artys
#' character strings \code{[]};
#' \code{"all"} or output array types to read (e.g. "FlowRightFace",
#'  "RiverLeakage") (not currently implemented)
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
                        nf.to.NA = FALSE, hds,
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
