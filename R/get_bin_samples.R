get_bin_samples <- function(X, 
                            node.stem,
                            node.extras = c(".pop", ".new"),
                            study.labs = NULL,
                            var.name = NULL,
                            val.name = NULL){
  
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
    colnames(node_sims) <- study.labs
  }

  # provide output in long format (data.frame)
  samples <- melt(node_sims)[2:3] # drop row ids
  samples$Var2 <- as.character(samples$Var2)
  if (!is.null(var.name)){
    colnames(samples)[1] <- var.name  
  }
  if (!is.null(val.name)){
    colnames(samples)[2] <- val.name  
  }
  
  out <- samples  
  return(out)
}
