# Rflow package - functions for the GHB MODFLOW package

#' Read general head boundary (GHB) package file
#'
#' Reads information from a MODFLOW GHB package file.  The GHB package
#'  specifies the simplest kind of head dependent flux boundary condition
#'  that includes a conductance term.  The flux is calculated by:\deqn{
#'    Q_{GHB} = \frac{h_{model} - h_{GHB}}{Cond}
#'  }
#'
#' @param filename
#' character string;
#' the package file to read
#' @param nSP
#' integer \code{[1]} or DIS.MFpackage object;
#' number of stress periods or corresponding DIS package from which this
#'  can be read
#' @param FREE
#' logical \code{[1]};
#' if TRUE, read assuming free format; normally FREE = TRUE is okay even if
#'  the file is in fixed format, providing entries have a space between
#'  them
#'
#' @return
#' object of class GHB.MFpackage:\cr
#'   \code{$header} (1-row data.frame):\cr
#'     \code{..$MXACTB} (int): the maximum number of GHB boundary cells in any
#'       stress period\cr
#'     \code{..$IGHBCB} (int): unit number to which to save flow to GHB cells\cr
#'   \code{$spheaders} (data.frame with \code{nSP} rows):\cr
#'     \code{..$ITMP} (int): number of active GHB cells in each stress period
#'      or, if \eqn{<1}, signals to reuse information from previous stress
#'      period\cr
#'     \code{..$NP} (int): number of parameters active in each stress period (not
#'       that parameter-defined boundaries are not currently supported in
#'       \code{read.GHB})\cr
#'   \code{$data} (data.table):\cr
#'     \code{..$L} (int): layer\cr
#'     \code{..$R} (int): row\cr
#'     \code{..$C} (int): column\cr
#'     \code{..$BHead} (double): head at boundary above model datum
#'      (\eqn{h_{GHB}})\cr
#'     \code{..$Cond} (double): conductance of boundary (units
#'      \eqn{length^2/time})\cr
#'     \code{..$...}: any auxiliary variables that are included in the file
#'      but which don't play a role in the MODFLOW simulation
#'
#' @import data.table
#' @export
#'
#' @examples
#'
read.GHB <- function(filename, nSP, FREE = TRUE){
  # get number of stress periods
  nSP <- switch(class(nSP),
                integer = nSP,
                numeric = as.integer(nSP),
                DIS.MFpackage = nrow(nSP$sps),
                stop("number of stress periods not given"))

  con <- file(filename, "rt")
  on.exit(close(con))

  # find how many lines to skip
  skip <- 0L
  while({
    ln1 <- readLines(con, 1L)
    is.comment(ln1, "#") || grepl("PARAMETER", ln1, ignore.case = TRUE)
  }) skip <- skip + 1L

  # start over and skip
  close(con)
  on.exit(NULL)
  con <- file(filename, "rt")
  on.exit(close(con))
  if(skip) readLines(con, skip)

  header <- read.table(con, nrows = 1L)
  names(header)[1:2] <- c("MXACTB", "IGHBCB")
  aux <- if(length(header) > 2L) sapply(header[-(1:2)], as.character)

  auxs <- if(!is.null(aux)){
    unname({
      aux[auxread <-
            which(aux %in% c("AUX", "AUXILIARY", "aux", "auxiliary")) + 1L]
    })
  }else character(0L)

  if(!is.null(aux)){
    names(header)[2L + auxread] <- paste0("AUX", 1:length(auxread))
    header <- header[, c(1:2, 2L + auxread)]
  }

  spheaders <- data.frame(ITMP = integer(nSP), NP = integer(nSP))
  ghb <- vector("list", nSP)

  for(sp in 1:nSP){
    sphd <- if(FREE){
      read.table(con, nrows = 1L, col.names = c("ITMP", "NP"))
    }else{
      read.fwf(con, c(10L, 10L), n = 1L, col.names = c("ITMP", "NP"))
    }
    if(sphd$NP > 0L) stop({
      "read.GHB does not yet support boundaries specified as parameters"
    })

    spheaders[sp,] <- sphd

    ghb[[sp]] <- if(FREE){
      read.table(con, nrows = sphd$ITMP,
                 col.names = c("L", "R", "C", "BHead", "Cond", auxs))
    }else{
      read.fwf(con, rep(10L, 5L + length(auxs)), n = sphd$ITMP,
               col.names = c("L", "R", "C", "BHead", "Cond", auxs))
    }
    setDT(ghb[[sp]])
    ghb[[sp]][, "sp" := sp]
  }

  data <- rbindlist(ghb)
  setcolorder(data, c("sp", "L", "R", "C", "BHead", "Cond", auxs))
  setkey(data, sp)

  structure(list(header = header, spheaders = spheaders,
                 data = data), class = "GHB.MFpackage")
}

#' Write a general head boundary (GHB) package file
#'
#' Writes information from a GHB.MFpackage list object to a MODFLOW-readable
#'  GHB package file.
#'
#' @param GHB
#' GHB.MFpackage object, as would be returned by \code{\link{read.GHB}}
#' @param filename
#' character string;
#' file to write to
#' @param title
#' character string;
#' optional title to prefix to the file, for your reference (not used if
#'  \code{MF2k = FALSE})
#' @param MF2k
#' logical \code{[1]};
#' should the file be readable by MODFLOW 2000 and higher, or
#'  (\code{FALSE}) keep back-compatibility with MF96; in the case of the
#'  GHB package, this means no title lines
#'
#' @return
#' NULL
#'
#' @import data.table
#' @export
#'
#' @examples
#'
write.GHB <- function(GHB, filename, title, MF2k = TRUE,
                      write.auxiliaries = TRUE){
  con <- file(filename, "wt")
  on.exit(close(con))

  if(!missing(title) && MF2k){
    writeLines(paste("#", title), con)
  }

  if(MF2k) writeLines(c({
    "# MODFLOW 2000 general head boundary package file created by write.GHB function, in R"
  }, {
    "PARAMETER  0  0"
  }), con)

  # auxiliary parameters
  if(write.auxiliaries){
    auxs <- sapply(GHB$header[, -(1:2)], as.character)
  }

  with(GHB$header, writeLines({
    formatC(c(MXACTB, IGHBCB), width = 10L, digits = 0L, format = "d")
  }, sep = "", con))
  if(write.auxiliaries && ncol(GHB$header) > 2L){
    writeLines(paste0(" AUXILIARY ",
                      sapply(GHB$header[, -(1:2)], as.character)),
               sep = "", con)
  }
  writeLines("", con)

  nSP <- nrow(GHB$spheaders)

  setcolorder(GHB$data, c("sp", "L", "R", "C", "BHead", "Cond", auxs,
                          recursive = TRUE))
  if(write.auxiliaries && ncol(GHB$data) > 6L){
    ff.width <- rep(10L, ncol(GHB$data) - 1L)
    ff.digits <- c(0L, 0L, 0L, 3L, 3L,
                   sapply(sapply(GHB$data[, -(1:6), with = FALSE], class),
                          switch, numeric = 3L, 0L))
    ff.format <- c("d", "d", "d", "f", "f",
                   sapply(sapply(GHB$data[, -(1:6), with = FALSE], class),
                          switch, numeric = "f", integer = "d", "s"))

    cols <- 2:ncol(GHB$data)
  }else{
    ff.width <- c(10L, 10L, 10L, 10L, 10L)
    ff.digits <- c(0L, 0L, 0L, 3L, 3L)
    ff.format <- c("d", "d", "d", "f", "f")

    cols <- 2:6
  }

  for(spn in 1:nSP){
    with(GHB$spheaders[spn,], writeLines({
      formatC(c(ITMP, NP), width = 10L, digits = 0L, format = "d")
    }, sep = "", con))
    writeLines("", con)

    # can't use data.table's efficient `by` subsetting here because some
    #  stress periods may have no entries, which would result in them being
    #  omitted
    GHB$data[sp == spn, {
      ffs <- vapply(1:length(cols), function(coln){
        formatC(.SD[[cols[coln]]],
                width = ff.width[coln],
                digits = ff.digits[coln],
                format = ff.format[coln])
      }, character(.N))

      writeLines(apply(ffs, 1L, paste, collapse = ""), con)
    }]
  }

  invisible()
}
