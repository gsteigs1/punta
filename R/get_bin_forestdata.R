get_bin_forestdata <- function(X, node.stem,
                               study.labs = NULL,
                               node.extras = c(".pop", ".new"), 
                               cols = c("X50.", "X2.5.", "X97.5."),
                               col.labs = c("est", "ci_lo", "ci_up")){

  x_sum <- data.frame(node = rownames(X$BUGSoutput$summary),
                      X$BUGSoutput$summary,
                      stringsAsFactors = FALSE)

  # extract the study specific estimates and the columns of interest
  r_id <- grep(paste(node.stem, "[", sep = ""), x_sum$node, fixed = TRUE)
  c_id <- c("node", cols)
  x_sel <- x_sum[r_id, c_id]
  
  if(!is.null(study.labs)){
    x_sel$node <- study.labs
  }
  
  # add population estimate and new study prediction if requested

  
  if (!is.null(node.extras)){
    r_ext <- paste(node.stem, node.extras, sep = "")
    x_ext <- x_sum[r_ext, c_id]
  }
    
  out <- rbind(x_sel, x_ext)
  colnames(out) <- c("node", col.labs)
  
  return(out)
}
  