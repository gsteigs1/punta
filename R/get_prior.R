#' Extract the RE prior distribution and parameter values from character descriptors.
#'
#' @param X Character string defining the prior. Currently supported are: "U(a,b)" and "HN(sigma)".
#'
#' @return
#' @export
#'
#' @examples
get_prior <- function(X){
  split1 <- strsplit(X, split = "(", fixed = TRUE)[[1]]
  split2 <- strsplit(tail(split1, 1), split = ",", fixed = TRUE)[[1]]
  split3 <- strsplit(tail(split2, 1), split = ")", fixed = TRUE)[[1]]
  
  out <- list(dist = split1[1],
              par1 = ifelse(length(split2) == 2,
                            as.numeric(split2[1]),
                            0),
              par2 = as.numeric(split3[1]))
  
  return(out)
}
