#' Produce Markov trace from Weibull survival models for OS and PFS assuming same shape but different scale parameters.
#'
#' @param time Vector wtih time-points at which the trace is produced.
#' @param model.pars List with named elements for the Weibull model, `scale.os, scale.pfs, shape`. 
#' 
#' @details The Weibull survival model uses the parameterization: S(t) = exp(-scale * t ^ shape).
#' 
#' @author Adapted from the example publically available from the DARTH group: https://github.com/DARTH-git/Partitioned-Survival-Analysis
#'
#' @return
#' @export
#'
#' @examples
get_partsurv <- function(time, model.pars){
  # Weibull model
  os  <- exp(-model.pars$scale.os  * time ^ model.pars$shape)
  pfs <- exp(-model.pars$scale.pfs * time ^ model.pars$shape)
  
  # Partitioned survival model
  pre.prog <- pfs             
  post.prog <- os - pfs
  post.prog[post.prog < 0] <- 0
  death <- 1 - os
  trace <- data.frame(pre.prog = pre.prog, post.prog = post.prog, death = death)
  
  return(trace)  
}
