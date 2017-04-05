# Rflow package - finding the phreatic groundwater top

#' Phreatic groundwater top
#'
#' @param mfdata
#' NetCDF object (see \code{\link{GW.nc}});
#' MODFLOW results
#' @param fnm
#' character string;
#' file name to save to
#' @param nts.dtit
#' \code{"NTS"} (almost always) or \code{"sNTS"};
#' name of the dimension giving the number of time steps (\code{"sNTS"}) is
#'  used when the NetCDF has been split (see \code{\link{GW.nc}})
#'
#' @return
#' NetCDF object with "wtop" variable representing the top of saturated
#'  groundwater in each cell, either the water table or the cell top,
#'  whichever is lower
#'
#' @import RNetCDF
#' @export
#'
#' @examples
#' mfdata <- RNetCDF::open.nc(system.file("rflow_mf_demo.nc", package = "Rflow"))
#'
#' wtop <- get.wtop.nc(mfdata, "RFLOW_EXAMPLE_wtop.nc")
#' print.nc(wtop)
get.wtop.nc <- function(mfdata, fnm, nts.dtit = "NTS"){
  if(file.exists(fnm)){
    wtop <- open.nc(fnm)
    on.exit(close.nc(wtop))
  }else{
    wtop <- create.nc(fnm, large = TRUE)
    on.exit(close.nc(wtop))
    att.put.nc(wtop, "NC_GLOBAL", "title", "NC_CHAR",
               "height of saturated groundwater above datum")
    att.put.nc(wtop, "NC_GLOBAL", "history", "NC_CHAR",
               paste("Created on", date(), "by MassTrack"))

    ncol <- dim.inq.nc(mfdata, "NCOL")$length
    nrow <- dim.inq.nc(mfdata, "NROW")$length
    nlay <- dim.inq.nc(mfdata, "NLAY")$length
    dim.def.nc(wtop, "NCOL", ncol)
    dim.def.nc(wtop, "NROW", nrow)
    dim.def.nc(wtop, "NLAY", nlay)
    dim.def.nc(wtop, nts.dtit, nmfts <- dim.inq.nc(mfdata, nts.dtit)$length)

    var.def.nc(wtop, "wtop", "NC_FLOAT",
               c("NCOL", "NROW", "NLAY", nts.dtit))
    att.copy.nc(mfdata, "Head", "missing_value", wtop, "wtop")

    # lt is layer top
    lt <- c(var.get.nc(mfdata, "elev", count = c(NA, NA, nlay)))
    for(i in 1:nmfts){
      var.put.nc(wtop, "wtop", {
        # wt is water head
        wt <- c(var.get.nc(mfdata, "Head",
                           c(NA, NA, 1L, i), c(NA, NA, nlay, 1L)))

        # layer top or water head, whichever is less
        ifelse(is.na(wt), NA, ifelse(lt > wt, wt, lt))
      }, c(1L, 1L, 1L, i), c(ncol, nrow, nlay, 1L))
    }

    att.copy.nc(mfdata, "Head", "units", wtop, "wtop")
    att.copy.nc(mfdata, "NC_GLOBAL", "datum", wtop, "NC_GLOBAL")
  }

  on.exit(NULL)
  wtop
}
