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
#' is.comment(" #  ")
#' is.comment(" #  ", allow.ws = TRUE)
#' is.comment("! test", "!")
#' is.comment("code", "")
#'
is.comment <- function(line, comm.char = "#", allow.ws = FALSE){
  if(comment.char == "") return(FALSE)
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
#' numeric.extract("1 2  3", " ")
#'
#' # sep = "" allows any amount of white space between numbers
#' numeric.extract("1 2  3", "")
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
#' nicearname(" FLOW RIGHT FACE") # "FlowRightFace"
#'
nicearname <- function(arty){
  str_replace_all(str_to_title(arty), " ", "")
}
