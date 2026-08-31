
#pragma once

//// File Name: ddm_bayes_helpers.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

arma::colvec ddm_generate_alphai_rcpp( const arma::colvec& c_alpha, const arma::mat& sigma, 
    const arma::colvec& MU, const arma::mat& SIGMA, const double tau2,
    const double epsilon, const double pi_mix, const std::string type_proposal = "std" );

arma::vec ddm_update_a_rcpp( arma::vec a_d, const double nu0,
    const arma::mat& c_invSIGa, const int K );

arma::mat outer_arma(arma::mat x, arma::rowvec mean );

Rcpp::List ddm_update_sigmaalpha_rcpp( const arma::mat& alpha, 
    const arma::colvec& c_MUa, const int I, const int K,
    const int nu0, const arma::mat& S0, arma::vec a_d,
    const std::string type_proposal );