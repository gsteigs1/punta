#' Calculate total discounted costs and QALYs from Markov traces.
#'
#' @param time Vector with time points in years (used for discounting).
#' @param trace Markov trace (data.frame or matrix).
#' @param costs_per_cycle Vector with costs per cycle.
#' @param utils Vector with utilities per cycle.
#' @param cycle_length Cycle length in years (e.g. 1/12 if cycle length is months).
#' @param dr Discount rate per year (same for costs and utilities).
#' 
#' @author Adapted from the example publically available from the DARTH group: https://github.com/DARTH-git/Partitioned-Survival-Analysis
#'
#' @return
#' @export
#'
#' @examples
get_ee <- function(time, trace, costs_per_cycle, utils, cycle_length, dr){
  
  if(ncol(trace) != length(costs_per_cycle)){
    stop("Number of columns in Markov trace must match number of cost inputs.")
  }
  if(ncol(trace) != length(utils)){
    stop("Number of columns in Markov trace must match number of utility inputs.")
  }
  
  if(!is.matrix(trace)){
    trace <- as.matrix(trace)
  }
  
  # costs and utilities per cycle
  trace_c <- trace %*% costs_per_cycle   
  trace_u <- trace %*% utils
  
  # discounted costs and utilities
  dr_w <- 1 / (1 + dr) ^ (time)  # dr in years -> time must be in years
  tot_c_d <- t(dr_w) %*% trace_c # total discounted costs: sum of dicsounted costs per cycle
  tot_u_d <- t(dr_w) %*% trace_u * cycle_length # total discounted QALYs: AUC (multiply by cycle length)
  
  out <- c(total.costs = tot_c_d,
           total.QALYs = tot_u_d)
  
  return(out)
}
