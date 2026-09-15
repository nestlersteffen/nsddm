#' nsddmm: Fit triadic and dyadic round-robin data
#'
#' @description Provides functionality for maximum likelihood and Bayesian 
#'   estimation for the 4-parameter drift-diffusion model.
#'
#' @keywords internal
#' @useDynLib tsrm, .registration = TRUE
#' @useDynLib tsrm_TMBExports, .registration = TRUE
#' @importFrom TMB MakeADFun
#' @importFrom methods as
#' @import Rcpp
#' @import RcppArmadillo
#' @import stats
#' @import utils
#' @import mvtnorm
#' @import numDeriv
#' @import nloptr
#' @import rtdists
#' @import MCMCpack
"_PACKAGE"