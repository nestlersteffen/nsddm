
#pragma once

//// File Name: ddm_derivatives.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

double ddm_derivatives_rcpp( const arma::vec& MU, 
    const arma::mat& L, const arma::mat& invSIGMA, 
    const arma::rowvec& alphainvSIGMA, const arma::rowvec& pos, 
    const int type_jj );