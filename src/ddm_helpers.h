
#pragma once

//// File Name: ddm_helpers.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>    

bool ddm4_parmcheck( const double& a, const double& t0, const double& z, const double& w );

bool ddm7_parmcheck( const double& a, const double& t0, const double& z, 
    const double& sv, const double& st0, const double& sz );

double ddm_sumexplog_rcpp( Rcpp::NumericVector x);