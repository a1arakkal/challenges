
#' log_like_cpp - a function to compute the log likelihood for challenge 8
#' @return A scalar 
#' @details
#'    More summary needed here. 
#' 
#' @import Rcpp
#' @useDynLib challenges   
#' @export
log_like_cpp <- function(params, X, y, tpts) {
  .Call('_challenges_log_like_cpp', PACKAGE = 'challenges', params, X, y, tpts)
}

#' log_post_cpp - a function to compute the log posterior for challenge 8
#' @return A scalar 
#' @details
#'    More summary needed here. 
#' 
#' @import Rcpp
#' @useDynLib challenges   
#' @export
log_post_cpp <- function(params, X, y, tpts) {
  .Call('_challenges_log_post_cpp', PACKAGE = 'challenges', params, X, y, tpts)
}

#' MH - a function to run MH for challenge 8
#' @return A scalar for number of accepted draws and also modifies the matrix supplied to the first argument in place
#' @details
#'    More summary needed here. 
#' 
#' @import Rcpp
#' @useDynLib challenges   
#' @export
MH <- function(MH_draws, proposal_cov, X, y, tpts) {
  .Call('_challenges_MH', PACKAGE = 'challenges', MH_draws, proposal_cov, X, y, tpts)
}

