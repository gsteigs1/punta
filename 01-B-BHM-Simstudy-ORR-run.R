# ***************************************************** #
# Title: 01-B-BHM-Simstudy-ORR-run                      #
# Author: Sandro Gsteiger                               #
# Description: Simulation study to assess impact of:    #
#              - different design options/sample sizes  #
#                (number of groups, size of groups)     #
#              - different RE SD priors                 #
#              on predictive uncertainty.               #
#              Hypothetical setup with binary data.     #
#                                                       #
# This script runs the simulations. The file            #
#  01-B-BHM-Simstudy-ORR-report.Rmd                     #
# reads in the output and reports the results.          #
#                                                       #
# Working directory must be project main directory.     #
# ***************************************************** #

# ======================================================
# *** Setup ***
# ======================================================

rm(list = ls())

start.time <- Sys.time()

library(dplyr)
library(tidyr)

for(s in dir("R")){
  source(file.path(".", "R", s))
}

r_seed <- 8945334

# start all output with this
prefix <- "01-B-" 

# ======================================================
# *** Simulation scenarios ***
# ======================================================

## -----------------------------------------------------
## Scenarios
## -----------------------------------------------------

# True rate and heterogeneity
bin_pars <- list(
  p.pop = 0.5,
  re.sd = c(0.1, 0.3, 0.7)
)

bin_grid <- expand.grid(bin_pars)

# Size of lead group and of subsequent groups: scenarios with larger lead group
size_pars1 <- list(
  lead.grp = c("yes"),
  lead.grp.size = 20,
  subseq.grp.size = c(2, 5, 10, 20)
)

size_grid1 <- expand.grid(size_pars1) %>% 
  filter(lead.grp.size > subseq.grp.size)

# Size of lead group and of subsequent groups: scenarios without larger lead group
size_grid2 <- data.frame(
  lead.grp = c("no"),
  lead.grp.size = c(2, 5, 10, 20),
  subseq.grp.size = c(2, 5, 10, 20),
  stringsAsFactors = FALSE
) 

size_grid <- rbind(size_grid1, size_grid2)

# number of tissue types (groups)
grp_grid <- data.frame(n.grp = c(5, 7, 10, 15, 20))

# heterogeneity priors
prior_grid <- data.frame(
  dist = c("U(0,5)", "HN(0.5)", "HN(1)"),
  stringsAsFactors = FALSE
)

# construct the full grid of all simulation scenarios
full_grid <- bin_grid %>%
  crossing(size_grid) %>%
  crossing(grp_grid) %>%
  crossing(prior_grid)
  
dim(full_grid)
head(full_grid)
tail(full_grid)

# save to include in report
save(bin_grid, size_grid, grp_grid, prior_grid, full_grid, 
     file = paste("outputs/", prefix, "scenario-grid", ".RData", sep = ""))


## -----------------------------------------------------
## Global parameters
## -----------------------------------------------------

global_par <- list(
  n_sim = 500,
  n_chains = 3,
  n_iter =  6000,
  n_burnin = 1000,
  n_thin = 1,
  p.threshold = 0.3
)

# save to include in report
save(global_par, 
     file = paste("outputs/", prefix, "global_par", ".RData", sep = ""))

# ======================================================
# *** Run the simulation ***
# ======================================================

set.seed(r_seed)

## results data frame: for each scenario
##  2 nodes (p.new, re.sd) with 8 and 7 posterior summaries -> 15 items
##  x 6 simulation summary metrics
res <- data.frame(matrix(NA, nrow = 15 * 6 * nrow(full_grid), ncol = 5,
                  dimnames = list(NULL, c("scenario", "node", "var1", "var2", "value"))))

## run simulations
for(i in 1:nrow(full_grid)){ # outer loop: scenarios
  cat("Scenario ", i, " / ", nrow(full_grid), "\n")
  scen_i <- full_grid[i, ]
  
  sim_i <- data.frame(matrix(NA, ncol = 4, nrow = 15 * global_par$n_sim, 
                      dimnames = list(NULL, c("iteration", "node", "var1", "value"))))

  if (get_prior(scen_i$dist)$dist == "HN"){
    model_i <- "models/bin_abs_re_hn.jgmod"
    prior_i <- list(prior.fe.mean = boot::logit(0.3), prior.fe.prec = 1/10,
                    prior.re.mean = get_prior(scen_i$dist)$par1, 
                    prior.re.prec = 1 / (get_prior(scen_i$dist)$par2)^2 )
  }
  if (get_prior(scen_i$dist)$dist == "U"){
    model_i <- "models/bin_abs_re_unif.jgmod"
    prior_i <- list(prior.fe.mean = boot::logit(0.3), prior.fe.prec = 1/10,
                    prior.re.lo = get_prior(scen_i$dist)$par1,
                    prior.re.up = get_prior(scen_i$dist)$par2)
  }
  
  for(j in 1:global_par$n_sim){ # inner loop: n_sim iterations per scenario
    cat(".\n")
    sink(file = "temp.sink.txt") # sink R2jags model summary (as no quiet mode)
    sim_ij <- simstd_bin(par.sim = list(p.pop = scen_i$p.pop, 
                                        re.sd = scen_i$re.sd, 
                                        n.studies = scen_i$n.grp, 
                                        n.sizes = c(scen_i$lead.grp.size, 
                                                    rep(scen_i$subseq.grp.size, 
                                                        scen_i$n.grp - 1))), 
                         par.an = list(n_chains = global_par$n_chains, 
                                       n_iter = global_par$n_iter, 
                                       n_burnin = global_par$n_burnin, 
                                       n_thin = global_par$n_thin,
                                       monitor = c("p.new", "re.sd"),
                                       model.file = model_i,
                                       prior = prior_i,
                                       p.threshold = global_par$p.threshold), 
                         r.seed = NULL)
    sink(file = NULL) # stop
    
    sim_i[(j - 1)*15 + 1:15, ] <- data.frame(iteration = j, sim_ij)

    rm(sim_ij)
  }
  
  sum_i <- sim_i %>% group_by(node, var1) %>%
    summarize(mean = mean(value), sd = sd(value), Q1 = quantile(value, probs = 0.25), 
              Q3 = quantile(value, probs = 0.75), CIlo = quantile(value, probs = 0.025), 
              CIup = quantile(value, probs = 0.975)) %>%
    tidyr::gather(var2, value, 3:8)
  
  res[(i - 1) * (15 * 6) + 1:(15 * 6), ] <- data.frame(scenario = i, sum_i)

  cat("\n")
  rm(scen_i, model_i, prior_i, sim_i, sum_i)
}
dim(res)
head(res)
tail(res)


# save for use in reporting script
save(res, 
     file = paste("outputs/", prefix, "res", ".RData", sep = ""))


# ======================================================
# *** Session info ***
# ======================================================

# THIS FILE:
rstudioapi::getActiveDocumentContext()$path

# TODAY:
date()

# RUN TIME:  
end.time <- Sys.time()
print(paste('Run time =', round(as.numeric(end.time, units = "secs") - as.numeric(start.time, units = "secs"), 2), 'seconds', sep = ' '))


# SYSTEM AND PACKAGES:
sessionInfo()

# End of file
