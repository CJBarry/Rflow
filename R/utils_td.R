#Christopher Barry, started on 23/01/2015, at University of Birmingham

#This script provides a function 'td(d, m, y)' that returns an integer representing the date supplied, with td returning 2 on 01/01/1900.  In this way, a consistent day numbering system is easily supplied.  It is also designed to be consistent with Microsoft Office Excel, although in this case negative numbers are permitted for representing dates before 1900.  MS Excel has 1900 as a leap year, when in fact it was not (as a multiple of 100 but not 400), and so 2 represents this 01/01/1900 in this function in order to be consistent with MS Excel for most of its range (after February 1900).

#This script also provides a function 'monthdays(monthno, yearno)' which returns the number of days in a month

#The function works from the year 1200 initially, and then subtracts 255668, the number of days between 01/01/1200 eq. and 01/01/1900
#The function first of all determines how many leap years have passed since 1200 (noting that years that are multiples of 100 but not 400 are not leap years and so 2000 is but 1900 isn't a leap year). The number of days in the complete years since 1900 inclusive is then determined.
#The number of complete months in the current year is then determined and days added accordingly.
#The day number is then added.
#255668 is then subtracted as previously mentioned

# a %/% b returns the integer part of the division
# a %% b returns the remainder of the division

nlmd <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
nlpmd <- c(0, nlmd[1:11])
ly <- logical(400)
ly[seq(4, 400, 4)[-c(25, 50, 75)]] <- TRUE
d400y <- 365*length(ly) + sum(ly)

# works for any date, even in BC
td <- function(d = 1L, m = 1L, y){
  # days in preceding years (since 1 AD inclusive)
  lyfrag <- if(y %% 400 != 1) ly[1:((y - 1L) %% 400)] else logical(0)
  dpy <- (y - 1L)%/%400*d400y + 365*length(lyfrag) + sum(lyfrag)
  
  # days within year in preceding months
  dpm <- sum(nlpmd[1:m]) + ifelse(y %% 400 != 0,
                                  as.integer(ly[y %% 400] && m >= 3L),
                                  ifelse(m >= 3, 1, 0))
  
  # correction such that 1 is 31/12/1899, consistent with MS Excel after
  #  March 1900 (and more correct before): 693594
  corr <- 693594
  
  return(dpy + dpm + d - corr)
}

# copied from MFread's cellref.loc
find.ref.td <- function(x, gcs, rev = FALSE){
  refs <- as.integer(if(rev){
    length(gcs) - approx(gcs, 1:length(gcs), x, "constant", f = 0)$y
  }else{
    approx(gcs, 1:length(gcs), x, "constant", f = 0)$y
  })
  refs[x >= max(gcs)] <- NA_integer_
  refs
}

# revert back to character string date in format "d/m/yyyy"
# not yet vectorised
invtd <- function(x){
  x <- x + 693594
  
  # which set of 400 years is the date in?
  # Works for negative numbers.  Note that it returns the integer that is
  #  lower than the result of straight division, rather than the smaller
  #  integer.  Thus, whereas 3%/%2 = 1, (-3)%/%2 = -2
  y400.set <- x%/%d400y
  
  # number of days through the 400 year set
  d.400yset <- x%%d400y
  
  # current year within the 400 year set
  curryr.400yset <- find.ref.td(d.400yset,
                                cumsum(c(0, ifelse(ly, 366, 365))) + .5)
  
  # number of days in the past years within the set
  # - cy4ys.alt created because ly[0] gives a zero-length result, with
  #    confusing (and incorrect) results
  cy4ys.alt <- curryr.400yset
  cy4ys.alt[cy4ys.alt == 1L] <- NA
  #
  ydays.past <- ifelse(curryr.400yset == 1L, 0L, {
    sum(ifelse(ly[1:(cy4ys.alt - 1L)], 366, 365))
  })
  
  # days into current year
  drem <- x - y400.set*d400y - ydays.past
  
  # current month
  cmdays.y <- list(nl = cumsum(c(0, monthdays(1:12))),
                   l = cumsum(c(0, monthdays(1:12, 1996L))))
  currm.y <- mapply(function(d, y){
    find.ref.td(d, cmdays.y[[ly[y] + 1L]] + .5)
  }, drem, curryr.400yset)
  
  # days into current month (i.e. the date)
  # - for cmy, see comment for cy4ys
  cmy <- currm.y
  cmy[cmy == 1L] <- NA
  drem <- drem - ifelse(currm.y == 1L, 0L, {
    cmdays.y[[ly[curryr.400yset] + 1L]][cmy]
  })
  
  paste(floor(drem), currm.y, curryr.400yset + y400.set*400L, sep = "/")
}

# this one takes significantly longer than the first td function
td2 <- function(d, m, y){
  x <- as.Date(paste(y, m, d, sep = "-"))
  return(julian(x, as.Date("1899-12-30")))
}

monthdays <- function(monthno, yearno = 1999L){
  nlmd[monthno] + (ly[(yearno %% 400) + 400*((yearno %% 400) == 0)] & monthno == 2)
}

yeardays <- function(yearno) 365 + (ly[(yearno %% 400) + 400*((yearno %% 400) == 0)])
{}
