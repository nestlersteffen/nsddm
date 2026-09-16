
#pragma once

//// File Name: ddm_lan.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

// weights matrices for lan:
inline arma::mat W1, W2, W3, W4;
inline arma::vec b1, b2, b3, b4;

arma::vec ddm4_lanll_dnn_rcpp( const arma::vec& rt, const arma::ivec& x, const double& a, 
     const double& t0, const double& z, const double& v, SEXP dnn );

void lan_load_weights_rcpp( const std::string& path, const std::string& ddm="four" );

SEXP lan_load_dnn_rcpp( const std::string& ddm );

double lan_forward_rcpp( const arma::vec& input );

arma::vec lan_forward_backward_rcpp( const arma::vec& input, const int& idx = 4 );

arma::vec ddm4_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v );

arma::vec ddm4_lanll_weights_batch_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v );

arma::vec ddm4_lanll_weights_grad_batch_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v, 
    const int& idx = 4 );

arma::vec ddm7_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v,
    const double& sv, const double& sz, const double& st0 );