get_bin_threshold <- function(X, 
                              thresholds,
                              node.stem,
                              node.extras = c(".pop", ".new"),
                              study.labs = NULL){
  
  # select nodes of interest
  x_sims <- X$BUGSoutput$sims.list
  node_sims <- x_sims[[node.stem]]
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
