#' Perform simulation/analysis round for binary data.
#'
#' @param par.sim Named list with input parameters for simulation model: p.pop, re.sd, n.studies, n.sizes (see function `sim_bin_re` for details).
#' @param par.an  Named list with input parameters for analysis model.
#' @param r.seed  Optional argument to set the R seed (default is NULL).
#'
#' @return
#' @export
#'
#' @examples
simstd_bin <- function(par.sim, par.an, r.seed = NULL){
  # simulate data
  dsim <- sim_bin_re(p.pop = par.sim$p.pop, re.sd = par.sim$re.sd, n.studies = par.sim$n.studies, n.sizes = par.sim$n.sizes, r.seed = r.seed)
  
  # fit RE MA model (with JAGS)
  djags <- list(ns = dsim$n.studies,
                r = dsim$x,
                n = dsim$n)
  
  if(!is.null(r.seed)){
    set.seed(3 * r.seed)
  }
  
  fit <- R2jags::jags(data = c(djags, par.an$prior),
                      n.chains = par.an$n_chains,
                      n.iter = par.an$n_iter,
                      n.burnin = par.an$n_burnin,
                      n.thin = par.an$n_thin,
                      parameters = par.an$monitor,
                      inits = NULL,
                      model.file = par.an$model.file)
  
  # extract results
  out <- data.frame(node = c("p.new", "re.sd"),
                    rbind(get_quantiles(fit, "p.new"),
                          get_quantiles(fit, "re.sd")),
                    check.names = FALSE)
  rownames(out) <- NULL
  return(out)
}