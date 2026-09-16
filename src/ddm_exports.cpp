//// File Name: ddm_exports.cpp
//// File Version: 0.01

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "ddm4_density.h"
#include "ddm7_density.h"
#include "ddm_lan.h"
#include "ddm_mvnorm.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec ddm4_pdf_export( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& w = NA_REAL, const double& z = NA_REAL,
    const double& v = 0, const std::string& type_ddm = "std", const int& kmax = 5000, 
    const double& delta = 1e-29 ) {
    return ddm4_pdf_rcpp( rt, x, a, t0, w, z, v, type_ddm, kmax, delta );
}

// [[Rcpp::export]]
arma::vec ddm7_pdf_export( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v,
    const double& sv, const double& sz, const double& st0, const std::string& type_ddm = "std", 
    const int& kmax = 5000, const double& delta = 1e-29 ) {
    return ddm7_pdf_rcpp( rt, x, a, t0, z, v, sv, sz, st0, type_ddm, kmax, delta );
}

// [[Rcpp::export]]
void lan_load_weights_export( const std::string& path, const std::string& ddm ) {
    lan_load_weights_rcpp( path, ddm );
}

// [[Rcpp::export]]
SEXP lan_load_dnn_export( const std::string& ddm ) {
    return lan_load_dnn_rcpp( ddm );
}

// [[Rcpp::export]]
arma::vec ddm4_lanll_weights_export( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v ) {
    return ddm4_lanll_weights_rcpp( rt, x, a, t0, z, v );
}

// [[Rcpp::export]]
arma::vec ddm4_lanll_weights_batch_export( const arma::vec& rt, const arma::ivec& x,
    const double& a, const double& t0, const double& z, const double& v ) {
    return ddm4_lanll_weights_batch_rcpp( rt, x, a, t0, z, v );
}

// [[Rcpp::export]]
Rcpp::List ddm4_lanll_grad_weights_rcpp( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs ) 
{
    
    //- transform parameters:
    double v, a, z, t0;
    v  = alpha(0);
    a  = std::exp( alpha(1) ) + std::exp( alpha(2) );
    z  = std::exp( alpha(2) );
    t0 = std::exp( alpha(3) );

    //- get information for the loop:
    int n = rt.size();
    // arma::vec output(5, arma::fill::zeros);
    // arma::vec input(6);
    
    // //- loop:
    // for ( int i = 0; i < n; i++ ) {
    //     input(0) = a;
    //     input(1) = v;
    //     input(2) = t0;
    //     input(3) = z;
    //     input(4) = xs(i);
    //     input(5) = rt(i);
    //     output += lan_forward_backward_rcpp( input, 4 );
    // }
    arma::vec output = ddm4_lanll_weights_grad_batch_rcpp( rt, xs, a, t0, z, v, 4 );

    //- final step:
    arma::mat J( 4, 4, arma::fill::zeros );
    J(0,1) = exp( alpha( 1 ) ); 
    J(0,2) = z; J(1,0) = 1; J(2,3) = t0; J(3,2) = z;
    arma::vec grad_us = J.t()*output.subvec(1,4);

    //- make the list:
    return Rcpp::List::create(
        Rcpp::Named("objective") = -1*output(0),
        Rcpp::Named("gradient") = -1*grad_us
    ); 

}

// [[Rcpp::export]]
Rcpp::List ddm7_lanll_grad_weights_rcpp( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs ) 
{
    
    //- transform parameters:
    double v, a, t0, z, sv, sz, st0;
    v   = alpha(0);
    sv  = std::exp( alpha(4) );
    sz  = std::exp( alpha(5) );
    st0 = std::exp( alpha(6) );
    z   = std::exp( alpha(2) ) + 0.5*sz;
    t0  = std::exp( alpha(3) ) + 0.5*st0;
    a   = std::exp( alpha(1) ) + z + 0.5*sz;

    //- get information for the loop:
    int n = rt.size();
    arma::vec output(8, arma::fill::zeros);
    arma::vec input(9);
    
    //- loop:
    for ( int i = 0; i < n; i++ ) {
        input(0) = a;
        input(1) = v;
        input(2) = t0;
        input(3) = z;
        input(4) = sv;   
        input(5) = sz;   
        input(6) = st0;  
        input(7) = xs(i);
        input(8) = rt(i);
        output += lan_forward_backward_rcpp( input, 7 );
    }

    //- final step:
    arma::mat J( 7, 7, arma::fill::zeros );
    J(0,1) = exp( alpha(1) ); J(0,2) = exp( alpha(2) ); J(0,5) = sz; 
    J(1,0) = 1; J(2,3) = exp( alpha(3) ); J(2,6) = 0.5*st0; 
    J(3,2) = exp( alpha(2) ); J(3,5) = 0.5*sz; 
    J(4,4) = sv; J(5,5) = sz; J(6,6) = st0;
    arma::vec grad_us = J.t()*output.subvec(1,7);

    //- make the list:
    return Rcpp::List::create(
        Rcpp::Named("objective") = -1*output(0),
        Rcpp::Named("gradient") = -1*grad_us
    ); 

}

// [[Rcpp::export]]
Rcpp::List ddm_lan_nllfct_gradient_rcpp( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs, const arma::vec& MU, const arma::mat& SIGMA, 
    const arma::mat& invSIGMA, const std::string ddm="four" ) 
{
    
    //- get gradient of data with backward pass:
    Rcpp::List logData;
    if ( ddm == "four" ) {
        logData = ddm4_lanll_grad_weights_rcpp(alpha,rt,xs);
    } else {
        logData = ddm7_lanll_grad_weights_rcpp(alpha,rt,xs);
    }
    double logData_val = Rcpp::as<double>( logData["objective"] );
    arma::vec logData_grad = Rcpp::as<arma::vec>( logData["gradient"] );

    //- add value and gradient of the prior:
    arma::mat alpha_mat( 1, alpha.size() );
    alpha_mat.row(0) = alpha.t();
    double log_prior = arma::as_scalar( dmvnrm_arma( alpha_mat, MU.t(), SIGMA, true ) );
    arma::vec prior_grad = invSIGMA*(alpha-MU); 

    //- make the list:
    return Rcpp::List::create(
        Rcpp::Named("objective") = logData_val-log_prior,
        Rcpp::Named("gradient") =  logData_grad+prior_grad
    ); 

}