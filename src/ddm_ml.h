
#pragma once

//// File Name: ddm_ml.h
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h> 
#include "ddm_lan.h" 
#include "ddm_mvnorm.h" 
#include "ddm4_density.h"    
#include "ddm7_density.h"
#include "ddm_derivatives.h"

double ddm4_nllfct_export( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs, const std::string& type_alpha,
    const std::string& type_ddm, const int& kmax, const double& delta, 
    bool use_lan );

double ddm4_random_llfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha,
    const std::string& type_ddm,
    const int& kmax, const double& delta,
    bool use_lan);

double ddm4_random_nllfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha,
    const std::string& type_ddm,
    const int& kmax,
    const double& delta,
    bool use_lan );

double ddm7_nllfct_export( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs, const std::string& type_alpha,
    const std::string& type_ddm, const int& kmax, const double& delta, bool use_lan );

double ddm7_random_llfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha,
    const std::string& type_ddm,
    const int& kmax,
    const double& delta,
    bool use_lan );

double ddm7_random_nllfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha,
    const std::string& type_ddm,
    const int& kmax,
    const double& delta,
    bool use_lan );

Rcpp::List ddm_agh_gradfct_singleperson_rcpp( 
    const arma::vec& typevec,
    const arma::mat& posmat,
    const arma::vec& rti, const arma::ivec& xsi,
    const arma::mat& pts, const arma::vec& wgh,
    const arma::vec& MU, 
    const arma::mat& L, 
    const arma::mat& SIGMA,
    const arma::mat& invSIGMA,
    const std::string& type_alpha,
    const std::string& type_ddm,
    const std::string& ddm,
    const int& kmax,
    const double& delta,
    bool use_lan, std::string lan_path );