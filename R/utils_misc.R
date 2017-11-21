# Rflow package - miscellaneous utility functions


#' Check if line of text is a comment
#'
#' @param line
#' character string
#' @param comm.char
#' character string;
#' character (or string is permissible) to identify comment line
#' @param allow.ws
#' logical \code{[1]};
#' if FALSE, the comment character must be the first character, including
#'  spaces, for the line to be identified as a comment
#'
#' @return
#' logical \code{[1]}
#'
#' @examples
#' Rflow:::is.comment(" #  ")
#' Rflow:::is.comment(" #  ", allow.ws = TRUE)
#' Rflow:::is.comment("! test", "!")
#' Rflow:::is.comment("code", "")
#'
is.comment <- function(line, comm.char = "#", allow.ws = FALSE){
  if(length(line) > 1L){
    return(vapply(line, is.comment, logical(1L), comm.char, allow.ws))
  }
  if(length(comm.char) > 1L){
    return(any(vapply(comm.char, is.comment, logical(1L), line = line, allow.ws)))
  }
  if(identical(comm.char, "")) return(FALSE)
  if(allow.ws) line = trimws(line)
  substring(line, 1, nchar(comm.char)) %in% comm.char
}

#' Extract numerical values from a string
#'
#' @param x
#' character string;
#' @param sep
#' character;
#' 0 or 1 characters
#'
#' @return
#' numeric \code{[]}
#'
#' @examples
#' # sep = " " insists on exactly one space between numbers
#' Rflow:::numeric.extract("1 2  3", " ")
#'
#' # sep = "" allows any amount of white space between numbers
#' Rflow:::numeric.extract("1 2  3", "")
#'
numeric.extract <- function(x, sep = " "){
  scan(text = trimws(x), sep = sep, what = numeric(), quiet = T)
}

#' Format array headers from MODFLOW
#'
#' @param arty
#' character string
#'
#' @return
#' character string
#'
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_to_title
#'
#' @examples
#' Rflow:::nicearname(" FLOW RIGHT FACE") # "FlowRightFace"
#'
nicearname <- function(arty){
  str_replace_all(str_to_title(arty), " ", "")
}

#' Enforce recycling to a target length
#'
#' @param vec
#' @param target.length
#'
#' @return
#' vector \code{[target.length]}
#'
expand.vec <- function(vec, target.length){
  if(identical(length(vec), 1L)){
    return(rep(vec, target.length))
  }else return(vec)
}
