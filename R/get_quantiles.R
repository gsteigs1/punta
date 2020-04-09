#' Function to extract 95% CrI, quartiles, and median, IQR and CrI.range for a specific node from a jags fit.
#'
#' @param X Output from jags fit (object of class `rjags`).
#' @param node Name of the node that is to be summarized.
#'
#' @return
#' @export
#'
#' @examples
get_quantiles <- function(X, node){
  sd <- X$BUGSoutput$sims.list[[node]]
  qs <- quantile(sd, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  IQR <- qs[4] - qs[2]
  CrIR <- qs[5] - qs[1]
  out <- data.frame(t(qs), IQR, CrIR, check.names = FALSE)
  return(out)
}
