#' Calculate posterior probabilities for a node to be above (list of) threshold.
#'
#' @param X Output of jags fit (object of class `rjags`).
#' @param thresholds Vector of thresholds.
#' @param node.stem Name of the node to extract (e.g. "p").
#' @param node.extras Name extensions (e.g. ".new" if nodes "p.new" monitored and to be extracted).
#' @param study.labs Labels for the extracted node.stem chains (e.g. if p[1], p[2] etc. correspond to studies).
#'
#' @return
#' @export
#'
#' @examples
get_bin_threshold <- function(X, 
                              thresholds,
                              node.stem,
                              node.extras = c(".pop", ".new"),
                              study.labs = NULL){
  
  x_sims <- X$BUGSoutput$sims.list
  node_sims <- NULL
  
  # select nodes of interest
  if(!is.null(node.stem)){
    node_sims <- x_sims[[node.stem]]
  }
  if(!is.null(study.labs)){
    colnames(node_sims) <- study.labs
  }
  
  if(!is.null(node.extras)){
    for(i in node.extras){
      node_sims <- cbind(node_sims, x_sims[[paste(node.stem, i, sep = "")]])
    }
    study.labs <- c(study.labs, paste(node.stem, node.extras, sep = ""))
  }
  
  # calculate probability of being above thresholds
  samples_better <- array(NA, dim = c(dim(node_sims), length(thresholds)))
  for(i in 1:length(thresholds)){
    samples_better[,,i] <- node_sims >= thresholds[i]  
  }
  
  prob_better <- apply(samples_better, MAR = c(2, 3), FUN = mean)
  dimnames(prob_better) <- list(study.labs, paste("ProbAbove", thresholds, sep = ""))
  
  
  out <- data.frame(node = rownames(prob_better),
                    prob_better)
  return(out)
}
