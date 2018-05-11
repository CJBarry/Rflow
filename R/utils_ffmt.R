# Rflow package - number fixed formatters

# for character strings, with right justify as the default
FFa <- function(string, w = 10L, justify = "right", ...)
  format(string, width = w, justify = justify, ...)

# for integers only: pads integer as string with spaces so is exactly w characters
FFI <- function(value, w = 10){
  fmt <- formatC(value, 0L, w, "d")
  if(any(nchar(fmt) > w)) warning({
    "FFI: some numbers were too wide for w, these will be wider than expected\n"
  }, "w: ", w, "; widest value: ", max(nchar(fmt)))
  fmt
}

# for strings like "    1.0360", for value = 1.036, w = 10, d = 4
FFf <- function(value, w = 10, d = 4){
  fmt <- formatC(value, d, w, "f")

  if(any(nchar(value) > w)) warning({
    "FFf: some numbers were too wide for w, these will be wider than expected\n"
  }, "w: ", w, "; widest value: ", max(nchar(fmt)))

  fmt
}

#' Fixed format scientific
#'
#' @param value
#' @param w
#' @param d
#' @param ew
#'
#' @return
#' vector of character strings
#'
#' @importFrom stringr str_dup
#'
FFe <- function(value, w = 12L, d = 4L, ew = 2L){
  if(any(value != 0 & log10(abs(value)) >= 10^ew))
    warning("FFe: inappropriate values of ew used because exponents require more digits")

  # values which are too small for ew are reset to 0 without warning
  value[-log10(abs(value)) >= 10^ew] <- 0

  # simplest: ew = 2 is used by formatC, unless value is too large or small
  fmt <- formatC(value, d, w, "e")

  # for other values of ew
  if(ew == 2L) return(fmt) else if(ew > 2L){
    fmt <- vapply(value, function(v){
      # deal with 0 simply
      if(identical(as.numeric(value), 0)) return({
        paste0(str_dup(" ", w - d - ew - 4L), "0.", str_dup("0", d), "e+", str_dup("0", ew))
      })

      # if exponent would need to be as long as ew, then formatC will return correct result
      # nec.ew means necessary ew; predicts the apparent ew of the result from formatC
      if((nec.ew <- abs(log10(v))) >= ew) return(formatC(v, d, w, "e"))
      if(nec.ew < 2L) nec.ew <- 2L

      fmt.tmp <- formatC(v, d, w - ew + nec.ew, "e")
      return(paste0(substring(fmt.tmp, 1L, w - ew), str_dup("0", ew - nec.ew),
                    substring(fmt.tmp, w - ew + 1L, w - ew + nec.ew)))
    }, character(1L))

    return(fmt)
  }else if(ew == 1L){
    warning("ew = 1 not yet supported; returned value has ew = 2")
    return(fmt)
  }
}

FFlist <- list(i = function(value, w, d) FFI(value, w), f = FFf, e = FFe)

#' Generic fixed formats
#'
#' @param value
#' @param type
#' @param w
#' @param d
#'
#' @return
#' @importFrom stringr str_to_lower
#'
#' @examples
FFgen <- function(value, type, w = 10L, d = 4L){
  type <- str_to_lower(type)
  do.call(FFlist[[type]], list(value, w, d))
}

#' Read fixed width numbers
#'
#' @param text
#' character[];
#' strings to read
#' @param widths
#' integer[];
#' field width of numbers - may be a list of field widths
#'
#' @return
#'
#' @examples
read.fws <- function(text, widths){
  nc <- nchar(text)
  if(any(nc %% sum(widths) != 0)) stop("Rflow:::read.fws: widths do not fit into string lengths of text")
  text <- paste(text, collapse = "")
  widths <- rep(widths, times = nchar(text)/sum(widths))

  ends <- cumsum(widths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  vapply(seq_along(ends), function(i) substr(text, starts[i], ends[i]), character(1L))
}
