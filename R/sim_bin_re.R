#' Simulate aggregate level data from a random effects model for binary outcomes.
#'
#' @param p.pop True (population) event probability.
#' @param re.sd Random effects standard deviation.
#' @param n.studies Number of studies.
#' @param n.sizes Either a single integer giving the (same) sample size of all studies, or a vector of integers of length n.studies with the different group sizes.
#' 
#' @details The random effects model is put on the logit scale. This means the study specific event rates $p_i$ are related via the exchangeability model $$\logit(p_i)\sim N(\logit(p), re.sd^2)$$
#'
#' @return
#' @export
#'
#' @examples
sim_bin_re <- function(p.pop, re.sd, n.studies, n.sizes, r.seed = NULL){
  
  if(!is.null(r.seed)){
    set.seed(r.seed)
  }

  # population effect
  theta.pop <- boot::logit(p.pop)
  
  # study specific effects
  theta <- rnorm(n.studies, mean = theta.pop, sd = re.sd)
  p <- boot::inv.logit(theta)
  
  # study specific observed events
  if(length(n.sizes) == 1){
    n.sizes <- rep(n.sizes, n.studies)
  }
  x <- rbinom(rep(1, n.studies), size = n.sizes, prob = p)
  
  # output
  out <- list(x = x, 
              n = n.sizes,
              # rates.obs = x / n.sizes,
              p = p, 
              #theta = theta, 
              p.pop = p.pop, 
              # theta.pop = theta.pop, 
              n.studies = n.studies, 
              re.sd = re.sd, 
              r.seed = r.seed)
  return(out)
}