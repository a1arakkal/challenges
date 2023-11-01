
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
