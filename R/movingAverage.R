#' movingAverage.
#'
#' This function calculates a moving average.
#' @param x The vector
#' @param n The number of samples
#' @param centered If FALSE, then average current sample and previous (n-1) samples. If TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.) the colnames they should be changed to
#' @return data with fixed colnames
#' @export

movingAverage <- function(x, n=3, centered=FALSE) {
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  # If count less than x, lower n
  if(length(count) < n) { n <- length(count) }
  #
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  # Add the centered data
  new <- x
  # Add to count list wherever there isn't a
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i+1
  }
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    i <- i+1
  }
  # return sum divided by count
  s/count
}
