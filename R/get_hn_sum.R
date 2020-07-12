#' Calculate summaries of the half-normal distribution.
#'
#' @param sigma The scale parameter (>0).
#' @param probs Vector with the quantiles to be produced (default: 0.025, 0.25, 0.5, 0.75, 0.975)
#'
#' @return A named vector with mean, sd, and the quantiles.
#' @export
#'
#' @examples
get_hn_sum <- function(sigma, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){
  quants <- qnorm(0.5 + probs/2, mean = 0, sd = sigma)
  mean <- sigma * sqrt(2) / sqrt(pi)
  sd <- sigma * sqrt((1 - 2/pi))
  
  out <- c(mean = mean, sd = sd, quants)
  return(out)
}
