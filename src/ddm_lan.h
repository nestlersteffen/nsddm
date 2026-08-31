
#pragma once

#define EIGEN_DONT_PARALLELIZE

//// File Name: ddm_lan.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

// weights matrices for lan:
inline arma::mat W1, W2, W3, W4;
inline arma::vec b1, b2, b3, b4;

// weights matrices for lan (eigen, for fast batched forward pass):
// inline Eigen::MatrixXd W1_e, W2_e, W3_e, W4_e;
// inline Eigen::RowVectorXd b1_e, b2_e, b3_e, b4_e;

arma::vec ddm4_lanll_dnn_rcpp( const arma::vec& rt, const arma::ivec& x, const double& a, 
     const double& t0, const double& z, const double& v, const int& K, SEXP dnn );

void lan_load_weights_rcpp( const std::string& path, const std::string& ddm="four" );

double lan_forward_rcpp( const arma::vec& input );

//arma::vec lan_forward_batch_rcpp( const arma::mat& X );

arma::vec lan_forward_backward_rcpp( const arma::vec& input, const int& idx = 4 );

arma::vec ddm4_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v );

// arma::vec ddm4_lanll_weights_rcpp_new( const arma::vec& rt, const arma::ivec& x,
//     const double& a, const double& t0, const double& z, const double& v );

arma::vec ddm7_lanll_weights_rcpp( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v,
    const double& sv, const double& sz, const double& st0 );