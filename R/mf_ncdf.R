# Rflow package - converting binary MODFLOW results to NetCDF

#' MODFLOW results to NetCDF data set
#'
#' @param dir
#' character string;
#' directory in which reside the MODFLOW output and input files
#' @param mfrt
#' character string;
#' non-extended filename to use as root for finding default MODFLOW input and
#'  output file names
#' @param ncrt
#' character string;
#'  filename (minus extension) for resulting NetCDF data set (default is \code{mfrt})
#' @param DIS,BAS,LPF,BCF
#' character string;
#' File names for DIS, BAS, LPF and BCF package files.  The latter three
#'  are only used for finding \code{HNOFLO} and \code{HDRY} if these are
#'  not given.
#' @param HDS
#' character string;
#' File name for head array.
#' @param CBB,CBW,CRC,CRV,CBD,CBG,CS1
#' character string;
#' File names for budget files, cell-by-cell flows and those relating to
#'  specific boundary conditions.  If a file doesn't exist it is simply
#'  ignored.
#' @param HNOFLO
#' numeric \code{[1]};
#' value to interpret as inactive cells in head array, automatically gleaned
#'  from BAS package file if not given
#' @param HDRY
#' numeric \code{[1]};
#' value to interpret as dry cell in head array, automatically gleaned from LPF
#'  or BCF package file if not given
#' @param title
#' character string;
#' title of the model, optional but recommended
#' @param author
#' character string;
#' who constructed the MODFLOW model?  (Not who performed the conversion to
#'  NetCDF.)
#' @param datum
#' character string;
#' description of the reference datum for height/ head values, for example
#'  "Ordnance Datum"
#' @param length.unit,time.unit
#' character string;
#' implied length and time units of the model, recommended to use full name
#'  (e.g. \code{"metre"} not \code{"m"})
#' @param xyt0
#' numeric \code{[3]};
#' global x and y values at the bottom left corner of the model and the absolute
#'  time value that the start of the model refers to
#' @param start_date
#' character string;
#' purely explanatory, enables the author to describe what start_time refers to
#' @param updating
#' logical \code{[1]};
#' Get meta data (global attributes) from existing NetCDF data set?  Will
#'  only attempt if the NetCDF data set already exists, so the default is
#'  TRUE.  This can save rewriting meta data such as title and origin.
#' @param split.tss
#' integer \code{[]};
#' time steps at which to split into multiple NetCDFs; the end time step of the
#'  split group is to be given; an integer vector if this feature is to be used
#' @param ch
#' logical \code{[1]};
#' whether to save the constant head flux (in a simulation with no CHs, there
#'  will still be arrays for constant head in the budget files, which waste space
#'  as all the elements are 0; it is not possible to know \emph{a priori} whether
#'  the CH array is needed, so the user should say if s/he knows)
#'
#' @return
#' \code{NULL}
#'
#' Prints a summary of the final NetCDF data set using print.nc.
#'
#' @import RNetCDF
#' @export
#'
#' @examples
#' # find the directory
#' wkdir <- paste0(system.file("inst", package = "Rflow"))
#'
#' # make fresh
#' GW.nc(wkdir, "rflow_mf_demo", "RFLOW_EXAMPLE", title = "example",
#'       author = "model builder's name", datum = "OD or something",
#'       xyt0 = c(0, 0, 0), start_date = "the start time descriptively")
#'
#' # update results, but keep metadata, with updating = TRUE (default)
#' GW.nc(wkdir, "rflow_mf_demo", "RFLOW_EXAMPLE")
#'
#' library(RNetCDF)
#' mfdata <- open.nc(paste0(wkdir, "/rflow_mf_demo.nc"))
#' print.nc(mfdata)
#'
GW.nc <- function(dir, mfrt, ncrt = mfrt,
                  DIS = paste0(mfrt, ".dis"), HDS = paste0(mfrt, ".hds"),
                  CBB = paste0(mfrt, ".cbb"), CBW = paste0(mfrt, ".cbw"),
                  CRC = paste0(mfrt, ".crc"), CRV = paste0(mfrt, ".crv"),
                  CBD = paste0(mfrt, ".cbd"), CBG = paste0(mfrt, ".cbg"),
                  CS1 = paste0(mfrt, ".cs1"), BAS = paste0(mfrt, ".bas"),
                  LPF = paste0(mfrt, ".lpf"), BCF = paste0(mfrt, ".bcf"),
                  HNOFLO, HDRY, title = "untitled MODFLOW data set",
                  author = "", datum = "", length.unit = "metre",
                  time.unit = "day", xyt0 = c(0, 0, 0), start_date = "",
                  updating = TRUE, .FillValue = 1e30, split.tss = NULL,
                  ch = TRUE){
  od <- getwd()
  setwd(dir)
  on.exit(setwd(od), add = TRUE)
  spl <- !is.null(split.tss) && class(split.tss) != "spl.instructions"
  spl.mode <- class(split.tss)[1L] == "spl.instructions"

  if(!spl.mode){

    ## if updating, get global attributes from existing data set, unless
    ##  they have been explicitly provided
    suff <- if(updating && spl && file.exists(paste0(ncrt, "_1.nc")))
      "_1" else ""
    if(updating && file.exists(paste0(ncrt, suff, ".nc"))){
      nc <- open.nc(paste0(ncrt, suff, ".nc"))
      if(missing(title)) title <- att.get.nc(nc, "NC_GLOBAL", "title")
      if(missing(author)) author <- att.get.nc(nc, "NC_GLOBAL", "author")
      if(missing(datum)) datum <- att.get.nc(nc, "NC_GLOBAL", "datum")
      if(missing(xyt0)){
        xyt0[1L] <- att.get.nc(nc, "NC_GLOBAL", "origin-x")
        xyt0[2L] <- att.get.nc(nc, "NC_GLOBAL", "origin-y")
        xyt0[3L] <- att.get.nc(nc, "NC_GLOBAL", "start_time")
      }
      close.nc(nc)
    }

    ## check that files exist before proceeding
    fnm.vec <- c(DIS = if(is.character(DIS)) DIS else "",
                 head = HDS, budget = CBB, wells = CBW,
                 recharge = CRC, river = CRV, drain = CBD, ghb = CBG,
                 stream = CS1)
    exst <- file.exists(fnm.vec)
    fi <- file.info(fnm.vec[exst])[, "ctime", drop = F]
    if(any(!exst[(if(is.character(DIS)) 1L else 2L):3]))
      stop(paste0("Not found files for: ",
                  str_c(names(fnm.vec[1:3])[!exst[1:3]]), collapse = ", "),
           "\n")

    # don't attempt to read files that aren't filled etc. (recognised by
    #  creation time - assumed that all binary save files finished being
    #  written within one minute; they are written simultaneously, a
    #  timestep at a time)
    cbbtime <- fi[CBB, "ctime"]
    oftimes <- logical(6)
    oftimes[(1:6)[exst[4:9]]] <- fi[(4:9)[exst[4:9]], "ctime"]
    exst[4:9] <- exst[4:9] &
      (oftimes > cbbtime - 60 & oftimes < cbbtime + 60)
    exst[is.na(exst)] <- FALSE

    # boundary condition flows to read
    BCflow.files <<- (fnm.vec[-(1:2)])[exst[-(1:2)]]
    on.exit(rm(BCflow.files, pos = .GlobalEnv), add = TRUE)

    ## read discretisation information and assign dimensions
    # also read values assigned to no flow and dry cells
    DIS <- switch(class(DIS),
                  character = read.DIS(DIS),
                  DIS.MFpackage = DIS,
                  stop({
                    "DIS must be character string giving path to the DIS package file, or a preloaded DIS package of class DIS.MFpackage, as returned by read.DIS"
                  }))

    if(missing(HNOFLO)){
      HNOFLO <- switch(class(BAS),
                       character = if(file.exists(BAS)){
                         read.BAS(BAS, DIS)$HNOFLO
                       }else 999,
                       BAS.MFpackage = BAS$HNOFLO,
                       stop({
                         "BAS must be character string giving path to the DIS package file, or a preloaded DIS package of class DIS.MFpackage, as returned by read.BAS"
                       }))
    }

    # find relative time values at the end of time steps
    mftime <- mftstime(DIS)

    if(missing(HDRY)){
      HDRY <- if(file.exists(LPF)){
        scan(LPF, list(integer(), double(), integer()), 1L,
             comment.char = "#")[[2L]]
      }else if(file.exists(BCF)){
        scan(BCF, list(integer(), double(), integer()), 1L,
             comment.char = "#")[[2L]]
      }else HNOFLO
    }

  }

  # recursive routine for when the data should be split into multiple
  #  NetCDFs
  # - creates a set of instructions to be passed back into the start of
  #    GW.nc
  if(spl){
    nspl <- length(split.tss)
    instruct <- Map(function(i, start, end){
      structure(list(suff = paste0("_", i),
                     startts = start,
                     endts = end),
                class = "spl.instructions")
    }, 1:nspl, c(0L, split.tss[-nspl]) + 1L, split.tss)

    sp.start_date <- if(i == 1L) start_date else if({
      "td" %in% rownames(installed.packages()) &&
        length.unit %in% c("d", "day")
    }) td::invtd(xyt0[3L] + mftime[instruct[[i - 1L]]$endts]) else ""

    for(i in 1:nspl){
      # should update to pass all arguments
      GW.nc(".", mfrt, paste(mfrt, i, sep = "_"), updating = FALSE,
            split.tss = instruct[[i]], title = title, datum = datum,
            xyt0 = c(xyt0[1:2], xyt0[3L] + if(i == 1L) 0 else{
              mftime[instruct[[i - 1L]]$endts]
            }),
            author = author, time.unit = time.unit,
            length.unit = length.unit,
            start_date = sp.start_date,
            HNOFLO = HNOFLO, HDRY = HDRY,
            DIS = DIS, BAS = BAS, BCF = BCF, LPF = LPF,
            HDS = HDS, CBB = CBB, CBW = CBW, CRC = CRC, CBD = CBD,
            CRV = CRV, CS1 = CS1, CBG = CBG)
    }

    return(invisible())
  }

  # create NetCDF
  nc <- create.nc(paste0(ncrt, ".nc"), large = TRUE)
  on.exit(close.nc(nc), add = TRUE)

  ## set global attributes
  att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", title)
  att.put.nc(nc, "NC_GLOBAL", "author", "NC_CHAR", author)
  att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR",
             paste("Created on", date(),
                   "by GW.nc function in R, coded by Christopher Barry"))
  att.put.nc(nc, "NC_GLOBAL", "origin-x", "NC_DOUBLE", as.double(xyt0[1L]))
  att.put.nc(nc, "NC_GLOBAL", "origin-y", "NC_DOUBLE", as.double(xyt0[2L]))
  att.put.nc(nc, "NC_GLOBAL", "datum", "NC_CHAR", datum)
  att.put.nc(nc, "NC_GLOBAL", "start_time", "NC_DOUBLE",
             as.double(xyt0[3L]))
  # note that the start date is purely explanatory, whereas the start time
  #  gives the numerical time value that should be taken to be the model
  #  start
  att.put.nc(nc, "NC_GLOBAL", "start_date", "NC_CHAR", start_date)
  att.put.nc(nc, "NC_GLOBAL", "note", "NC_CHAR",
             paste("co-ordinate variables (gccs, grcs, elev and time),",
                   "as well as the Head values, are relative to the",
                   "origin, datum and start time as appropriate"))
  if(spl.mode){
    att.put.nc(nc, "NC_GLOBAL", "subset", "NC_CHAR",{
      paste("subset of data from time step",
            split.tss$startts, "to", split.tss$endts)
    })
    # oddly, att.put.nc only accepts double values for the value argument,
    #  even when type = "NC_INT"
    att.put.nc(nc, "NC_GLOBAL", "subset_start_ts", "NC_INT",
               as.double(split.tss$startts))
    att.put.nc(nc, "NC_GLOBAL", "subset_end_ts", "NC_INT",
               as.double(split.tss$endts))
  }

  # define dimensions
  dim.def.nc(nc, "NCOL", DIS$extent["NCOL"])
  dim.def.nc(nc, "NROW", DIS$extent["NROW"])
  dim.def.nc(nc, "NLAY", DIS$extent["NLAY"])
  dim.def.nc(nc, "NTS", sum(DIS$sps$NSTP))
  if(spl.mode){
    dim.def.nc(nc, "sNTS", split.tss$endts - split.tss$startts + 1L)
  }
  dim.def.nc(nc, "NCOL+1", DIS$extent["NCOL"] + 1L)
  dim.def.nc(nc, "NROW+1", DIS$extent["NROW"] + 1L)
  dim.def.nc(nc, "NLAY+1", DIS$extent["NLAY"] + 1L)
  # max string length
  dim.def.nc(nc, "msl", 16L)
  # number of variables
  dim.def.nc(nc, "vars", unlim = TRUE)

  # grid dimensions and layer elevations
  var.def.nc(nc, "gccs", "NC_DOUBLE", "NCOL+1")
  var.def.nc(nc, "grcs", "NC_DOUBLE", "NROW+1")
  var.def.nc(nc, "elev", "NC_FLOAT", c("NCOL", "NROW", "NLAY+1"))

  att.put.nc(nc, "gccs", "description", "NC_CHAR",
             "column divider co-ordinates")
  att.put.nc(nc, "gccs", "units", "NC_CHAR", length.unit)
  att.put.nc(nc, "grcs", "description", "NC_CHAR",
             "row divider co-ordinates")
  att.put.nc(nc, "grcs", "units", "NC_CHAR", length.unit)
  att.put.nc(nc, "elev", "description", "NC_CHAR",
             "layer elevations above datum")
  att.put.nc(nc, "elev", "units", "NC_CHAR", length.unit)
  nts <- if(spl.mode) "sNTS" else "NTS"

  # expand if necessary
  rsp <- rev(DIS$DELC)
  if(identical(names(rsp), "CNSTNT")) rsp <- rep(rsp, DIS$extent["NROW"])
  grcs <- cumsum(c(0, rsp))
  csp <- DIS$DELR
  if(identical(names(csp), "CNSTNT")) csp <- rep(csp, DIS$extent["NCOL"])
  gccs <- cumsum(c(0, csp))
  if(is.vector(elev <- DIS$elev)){
    elev <- sapply(elev, matrix,
                   DIS$extent["NCOL"], DIS$extent["NROW"],
                   simplify = "array")
  }

  var.put.nc(nc, "gccs", gccs)
  var.put.nc(nc, "grcs", grcs)
  var.put.nc(nc, "elev", elev)

  # how many time steps have there been in the previous stresss periods?
  prev.ts <- cumsum(c(0, DIS$sps$NSTP[-DIS$extent["NPER"]]))

  ## assign time and head arrays
  var.def.nc(nc, "sp_ts", "NC_CHAR", c("msl", nts))
  var.def.nc(nc, "time", "NC_DOUBLE", nts)
  var.def.nc(nc, "Head", "NC_FLOAT", c("NCOL", "NROW", "NLAY", nts))

  att.put.nc(nc, "sp_ts", "description", "NC_CHAR",
             "stress period_time step")
  att.put.nc(nc, "time", "description", "NC_CHAR", {
    "time since start of MODFLOW simulation at end of time step"
  })
  att.put.nc(nc, "time", "units", "NC_CHAR", time.unit)
  att.put.nc(nc, "Head", "description", "NC_CHAR",
             "piezometric head, above datum")
  att.put.nc(nc, "Head", "units", "NC_CHAR", length.unit)
  att.put.nc(nc, "Head", "missing_value", "NC_FLOAT", as.double(HNOFLO))
  att.put.nc(nc, "Head", "HDRY", "NC_FLOAT", as.double(HDRY))
  att.put.nc(nc, "Head", "note", "NC_CHAR", {
    "HNOFLO, the no-flow value for head is used for the missing_value"
  })

  ## read in heads and time
  to.read <- file(HDS, "rb")
  cat("Head save file:\n")

  dcounts <- c(DIS$extent[c("NCOL", "NROW")], 1L, 1L)
  oldts <- tsi <- 0L
  repeat{
    # meta data, including elapsed time and layer number
    tssp <- readBin(to.read, "integer", 2L, 4L)

    # has end of file been reached?
    if(!length(tssp)) break

    # cumulative time step number
    curts <- prev.ts[tssp[2L]] + tssp[1L]

    if(!spl.mode || (curts >= split.tss$startts && curts <= split.tss$endts)){
      # is this a new time step?
      if(curts != oldts) tsi <- tsi + 1L

      time <- readBin(to.read, "double", 2L, 4L)[2L]
      readChar(to.read, 16L)
      lay <- readBin(to.read, "integer", 3L, 4L)[3L]

      # data array for this layer and time step
      ar <- readBin(to.read, "double",
                    prod(DIS$extent[c("NCOL", "NROW")]), 4L)

      dstarts <- c(1L, 1L, lay, tsi)

      # write to the NetCDF dataset
      var.put.nc(nc, "sp_ts", paste(rev(tssp), collapse = "_"),
                 c(1L, dstarts[4L]), c(16L, 1L))
      var.put.nc(nc, "time", time, dstarts[4L], 1L)
      var.put.nc(nc, "Head", ar, dstarts, dcounts)
      cat(".")
      oldts <- curts
    }else{
      # skip this array, or stop reading if got the necessary time steps
      if(curts > split.tss$endts) break
      readBin(to.read, "double", 2L, 4L)
      readChar(to.read, 16L)
      readBin(to.read, "integer", 3L, 4L)
      readBin(to.read, "double", prod(DIS$extent[c("NCOL", "NROW")]), 4L)
    }
  }; cat("\n")
  rm(ar); close(to.read)

  ## read flux data
  artys <- character(0L)
  cat("Flux files (the first arrays may be slow):\n")

  dcounts <- c(DIS$extent[c("NCOL", "NROW", "NLAY")], 1L)
  oldts <- tsi <- 0L
  for(fnm in BCflow.files){
    to.read <- file(fnm, "rb")

    repeat{
      tssp <- readBin(to.read, "integer", 2L, 4L)

      # end of file?
      if(!length(tssp)) break

      # cumulative time step number
      curts <- prev.ts[tssp[2L]] + tssp[1L]

      if(!spl.mode || (curts >= split.tss$startts && curts <= split.tss$endts)){
        # is this a new time step?
        if(curts != oldts) tsi <- tsi + 1L

        # data type (arty means array type)
        arty <- nicearname(readChar(to.read, 16L))
        if(!arty %in% artys && !(arty == "ConstantHead" && !ch)){
          # this a growing object, but only small and only a bit
          artys <- c(artys, arty)

          var.def.nc(nc, arty, "NC_FLOAT", c("NCOL", "NROW", "NLAY", nts))
          att.put.nc(nc, arty, "_FillValue", "NC_FLOAT",
                     as.double(.FillValue))
          att.put.nc(nc, arty, "units", "NC_CHAR",
                     paste0(length.unit, "^3 per ", time.unit))
        }

        # array dimensions (know already)
        readBin(to.read, "integer", 3L, 4L)

        ar <- readBin(to.read, "double",
                      prod(DIS$extent[c("NCOL", "NROW", "NLAY")]), 4L)

        dstarts <- c(1L, 1L, 1L, tsi)
        if(!(arty == "ConstantHead"&& !ch)){
          var.put.nc(nc, arty, ar, dstarts, dcounts)
          cat(".")
        }
        oldts <- curts
      }else{
        # skip this array, or stop reading if got the necessary time steps
        if(curts > split.tss$endts) break
        readChar(to.read, 16L)
        readBin(to.read, "integer", 3L, 4L)
        readBin(to.read, "double",
                prod(DIS$extent[c("NCOL", "NROW", "NLAY")]), 4L)
      }
    }
    close(to.read)
  }; cat("\n")
  rm(ar)

  # which data arrays are present?
  var.def.nc(nc, "outvars", "NC_CHAR", c("msl", "vars"))

  att.put.nc(nc, "outvars", "description", "NC_CHAR",
             "output variable arrays in this data set")

  var.put.nc(nc, "outvars", c("Head", artys),
             count = c(NA, length(artys) + 1L))

  print.nc(nc)
}

#' MT3DMS results to NetCDF
#'
#' @param dir
#' character string;
#' directory in which input and output files are stored
#' @param mtrt
#' character string;
#' default file name, minus extensions, for the data sets
#' @param ntts
#' integer \code{[]};
#' @param gw.nc
#' @param file
#' @param Nspecies
#' @param species.names
#' @param files
#' @param CNOFLO
#' @param title
#' @param author
#' @param conc.unit
#' @param length.unit
#' @param time.unit
#' @param updating
#'
#' @return
#' \code{NULL}
#'
#' prints summary of NetCDF file created
#'
#' @import RNetCDF
#' @export
#'
#' @examples
#'
MT3DMS.nc <- function(dir, mtrt, ntts = Inf, gw.nc, file = paste0(mtrt, ".nc"),
                      Nspecies = 10L, species.names = as.character(1:Nspecies),
                      files = cbind(mob = paste0(mtrt, 1:Nspecies, ".ucn"),
                                    immob = paste0(mtrt, "_sorbed", 1:Nspecies, ".ucn")),
                      CNOFLO = -1, title = "untitled MT3DMS data set",
                      author = "", conc.unit = "kg/m^3",
                      length.unit = "metre", time.unit = "day",
                      updating = TRUE){
  od <- getwd()
  setwd(dir)
  on.exit(setwd(od), add = TRUE)

  mobile <- if(NCOL(files) > 1L){
    cbind(!logical(nrow(files)), logical(nrow(files)))
  }else !logical(length(files))

  any.immob <- NCOL(files) > 1L
  files <- as.matrix(files)
  if(ncol(files) == 1L) files <- cbind(files, "")

  if(updating && file.exists(file)){
    nc <- open.nc(file)
    if(missing(title)) title <- att.get.nc(nc, "NC_GLOBAL", "title")
    if(missing(author)) author <- att.get.nc(nc, "NC_GLOBAL", "author")
    if(missing(length.unit)) length.unit <- att.get.nc(nc, "gccs", "units")
    if(missing(time.unit)) time.unit <- att.get.nc(nc, "time", "units")
    close.nc(nc)
  }

  nc <- create.nc(file)
  on.exit(print.nc(nc), add = TRUE)
  on.exit(close.nc(nc), add = TRUE)

  # get co-ordinate data from the MODFLOW data set
  gwdata <- switch(class(gw.nc),
                   character = {
                     on.exit(close.nc(gwdata), add = TRUE)
                     open.nc(gw.nc)
                   },
                   NetCDF = gw.nc,
                   stop({
                     "gw.nc should be a NetCDF dataset or a character string path to one"
                   }))

  # copy spatial dimensions
  ds.to.copy <- c("NCOL", "NROW", "NLAY", "NCOL+1", "NROW+1", "NLAY+1")
  ds.info <- lapply(ds.to.copy, dim.inq.nc, ncfile = gwdata)
  l_ply(ds.info, function(lst) dim.def.nc(nc, lst$name, lst$length))

  nC <- dim.inq.nc(gwdata, "NCOL")$length
  nR <- dim.inq.nc(gwdata, "NROW")$length
  nL <- dim.inq.nc(gwdata, "NLAY")$length

  # new attributes
  att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", title)
  att.put.nc(nc, "NC_GLOBAL", "author", "NC_CHAR", author)
  att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR",
             paste("Created on", date(),
                   "by MT3DMS.nc function in R, coded by Christopher Barry"))
  att.put.nc(nc, "NC_GLOBAL", "note", "NC_CHAR", {
    "co-ordinate variables (gccs, grcs, elev and time) are relative to the origin, datum and start time as appropriate"
  })

  # copy certain attributes
  att.copy.nc(gwdata, "NC_GLOBAL", "origin-x", nc, "NC_GLOBAL")
  att.copy.nc(gwdata, "NC_GLOBAL", "origin-y", nc, "NC_GLOBAL")
  att.copy.nc(gwdata, "NC_GLOBAL", "datum", nc, "NC_GLOBAL")
  att.copy.nc(gwdata, "NC_GLOBAL", "start_time", nc, "NC_GLOBAL")
  try(att.copy.nc(gwdata, "NC_GLOBAL", "start_date", nc, "NC_GLOBAL"), TRUE)

  # define new dimensions
  dim.def.nc(nc, "NTTS", unlim = TRUE)
  dim.def.nc(nc, "msl", 16L)

  # copy co-ordinate information
  var.def.nc(nc, "gccs", "NC_DOUBLE", "NCOL+1")
  var.def.nc(nc, "grcs", "NC_DOUBLE", "NROW+1")
  var.def.nc(nc, "elev", "NC_DOUBLE", c("NCOL", "NROW", "NLAY+1"))

  var.put.nc(nc, "gccs", var.get.nc(gwdata, "gccs"))
  var.put.nc(nc, "grcs", var.get.nc(gwdata, "grcs"))
  var.put.nc(nc, "elev", var.get.nc(gwdata, "elev"))

  # time and time point labels
  var.def.nc(nc, "time", "NC_DOUBLE", "NTTS")
  var.def.nc(nc, "sp_ts_tts", "NC_CHAR", c("msl", "NTTS"))

  att.put.nc(nc, "sp_ts_tts", "description", "NC_CHAR",
             "stress period_flow time step_transport time step")

  # copy units and descriptions
  for(vn in c("gccs", "grcs", "elev", "time")){
    try(att.copy.nc(gwdata, vn, "description", nc, vn), TRUE)
    try(att.copy.nc(gwdata, vn, "units", nc, vn), TRUE)
  }

  times <- double(1000L) + NA
  ttstssps <- character(1000L); ttstssps[] <- NA_character_

  first <- TRUE
  for(i in 1:Nspecies){
    for(mim in c("mobile", "immobile")){
      vn <- paste0("C_", species.names[i], "_", substr(mim, 1L, 1L))
      if(file.exists(files[i, switch(mim, mobile = 1L, immobile = 2L)])){
        cat("found results for", mim, "species:", species.names[i], "\n")

        var.def.nc(nc, vn, "NC_FLOAT",
                   c("NCOL", "NROW", "NLAY", "NTTS"))
        att.put.nc(nc, vn, "longname", "NC_CHAR",
                   paste("concentration of", species.names[i]))
        # att.put.nc(nc, vn, "mobile", "NC_INT",
        #            switch(mim, mobile = 1L, immobile = 0L))
        att.put.nc(nc, vn, "mobileTXT", "NC_CHAR", mim)
        att.put.nc(nc, vn, "units", "NC_CHAR", conc.unit)
        att.put.nc(nc, vn, "missing_value", "NC_FLOAT", CNOFLO)

        to.read <- file(files[i, switch(mim, mobile = 1L, immobile = 2L)],
                        "rb")

        an <- 0L
        repeat{
          ttstssp <- readBin(to.read, "integer", 3L, 4L)
          if(!length(ttstssp) || an > ntts - 1L) break
          time <- readBin(to.read, "double", 1L, 4L)

          readChar(to.read, 16L) # "   CONCENTRATION"
          lay <- readBin(to.read, "integer", 3L, 4L)[3L] # C, R, L

          # update time details
          if(lay == 1L) an <- an + 1L
          if(first){
            ttstssps[an] <- paste(rev(ttstssp), collapse = "_")
            times[an] <- time
          }

          ar <- readBin(to.read, "double", nC*nR, 4L)

          var.put.nc(nc, vn, ar,
                     c(1L, 1L, lay, an), c(nC, nR, 1L, 1L))
          cat(".")

        }; close(to.read); cat("\n")

        first <- FALSE
      }else{
        cat("not found results for", mim, "species:", species.names[i], "\n")
      }
    }
  }

  var.put.nc(nc, "time", times[!is.na(times)])
  var.put.nc(nc, "sp_ts_tts", ttstssps[!is.na(ttstssps)])
  invisible()
}
