#' Calculate shrinkage factors from a Bayesian MA (posterior summaries from full MCMC of the shrinkage factor).
#'
#' @param X Bayesian MA fit (object of class jags).
#' @param node.pop Name of the node with the population level parameter.
#' @param node.grp Name of the node with the group level parameter.
#' @param id.grp Interger id of the group level node for which to calculate the shrinkage factors.
#' @param summaries The numerical summaries to extract from the posteriors of the shrinkage factors.
#' @return
#' @export
#'
#' @examples
get_shrinkage <- function(X, node.pop, node.grp, id.grp, summaries = c("mean", "median")){

  if(length(id.grp) != 1){
    stop("id.grp must be a single integer")
  }
  
  sims <- X$BUGSoutput$sims.list

  sims.pop <- as.vector(sims[[node.pop]])
  sims.grp <- sims[[node.grp]][,id.grp]
  
  r <- X$model$data()$r[id.grp]
  n <- X$model$data()$n[id.grp]
  obs1 <- boot::logit(r/n)
  obs <- rep(obs1, length(sims.grp))
  
  sims.shrinkage <- (sims.grp - obs) / (sims.pop - obs)
  
  sums <- list()
  for (i in 1:length(summaries)){
    sums_i <- do.call(summaries[i], list(sims.shrinkage))
    sums[[i]] <- sums_i
  }  
  names(sums) <- summaries
  
  out <- data.frame(tissue = id.grp,
                    as.data.frame(sums)) 
  return(out)
}
