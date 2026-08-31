
#pragma once

//// File Name: ddm4_density.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>   

inline const int MAX = 2147483647;

double ddm4_dstand_rcpp( const double& t, const double& a, const double& w, const double& v, 
    const double& s2, const int& kmax, const double& delta );

int get_ks( const double& t, const double& err );

int get_kl( const double& t, const double& err );

double ddm4_dnavfuss_rcpp( const double& t, const double& a, const double& w, 
    const double& v, const double& err );

arma::vec ddm4_pdf_rcpp( const arma::vec& rt, const arma::ivec& x, 
  const double& a, const double& t0, const double& w = NA_REAL, const double& z = NA_REAL, 
  const double& v = 0, const std::string& type_ddm = "std", 
  const int& kmax = 5000, const double& delta = 1e-29 );

double ddm4_logLdata_rcpp( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& x, const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std", const int& kmax = 5000, 
    const double& delta = 1e-29, bool use_lan = false );