
#pragma once

//// File Name: ddm4_density.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

double ddm7_large_integrand_rcpp( const double& tx, 
    const double& a, const double& z, const double& v, 
    const double& t0, const double& sv, const int& kmax, 
    const double& delta );

double ddm7_small_integrand_rcpp( const double& tx, 
    const double& a, const double& z, const double& v, 
    const double& t0, const double& sv, const int& kmax, 
    const double& delta );

double ddm7_dao_rcpp( const double& tx, 
    const double& a, const double& z, const double& v, 
    const double& t0, const double& sv, const double& sz,
    const double& st0, const int& kmax, const double& err = 1e-6 );

double ddm7_dstand_rcpp( const double& tx, const int& x, 
    const double& a, const double& t0, const double& z, 
    const double& v, const double& st0, const double& sz, 
    const double& sv, const double& s2, const int& kmax, 
    const double& delta, const double& epsilon, 
    const double& min_RT );

arma::vec ddm7_pdf_rcpp( const arma::vec& rt, const arma::ivec& x, 
  const double& a, const double& t0, const double& z, const double& v, 
  const double& sv, const double& sz, const double& st0, 
  const std::string& type_ddm = "std", const int& kmax = 5000,
  const double& delta = 1e-29 );

double ddm7_logLdata_rcpp( const arma::colvec& alpha, 
    const arma::vec& rt, const arma::ivec& x, 
    const std::string& type_alpha = "dao",  // aktuell noch ungenutzt, aber für einheitliche Signatur
    const std::string& type_ddm = "std",      // aktuell noch ungenutzt, aber für einheitliche Signatur
    const int& kmax = 5000,
    const double& delta = 1e-29,
    bool use_lan = false );