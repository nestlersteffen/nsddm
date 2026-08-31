
#pragma once

//// File Name: ddm_mvnorm.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

inline const double log2pi = std::log(2.0 * M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);

arma::vec dmvnrm_arma( arma::mat x,  
                       arma::rowvec mean,  
                       arma::mat sigma, 
                       bool use_log = false );