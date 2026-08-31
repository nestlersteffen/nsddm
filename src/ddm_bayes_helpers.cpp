
#include "ddm_bayes_helpers.h"   

using namespace Rcpp;
using namespace arma;

arma::colvec ddm_generate_alphai_rcpp( const arma::colvec& c_alpha, const arma::mat& sigma, 
    const arma::colvec& MU, const arma::mat& SIGMA, const double tau2,
    const double epsilon, const double pi_mix, const std::string type_proposal ) 
{
    arma::colvec new_alpha;
    if ( type_proposal == "mixture" ) {
        if ( R::runif(0,1) < pi_mix ) {
            new_alpha = arma::mvnrnd( c_alpha, epsilon * SIGMA, 1 );
        } else {
            new_alpha = arma::mvnrnd( MU, SIGMA, 1 );
        }
    } else {
        new_alpha = arma::mvnrnd( c_alpha, tau2 * sigma, 1 );
    }
    return new_alpha;
}

arma::vec ddm_update_a_rcpp( arma::vec a_d, const double nu0,
    const arma::mat& c_invSIGa, const int K )
{
    double alpha_ig = ( nu0 + K ) / 2.0;
    
    for ( int d = 0; d < K; d++ ) {
        double beta_ig = nu0 * c_invSIGa(d,d) + 1; // 1.0 / ( a_d(d) * a_d(d) );
        //- sample from IG = 1/Gamma:
        a_d(d) = 1.0 / R::rgamma( alpha_ig, 1.0/beta_ig );
    }
    return a_d;
}

arma::mat outer_arma(arma::mat x, arma::rowvec mean ) 
{ 
    int n = x.n_rows; 
    int p = x.n_cols; 
    arma::rowvec z(p);
    arma::mat out(p,p, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        out += z.t()*z;     
    }  
    return out;
}

Rcpp::List ddm_update_sigmaalpha_rcpp( const arma::mat& alpha, 
    const arma::colvec& c_MUa, const int I, const int K,
    const int nu0, const arma::mat& S0, arma::vec a_d,
    const std::string type_proposal )
{
    //- compute covariance matrix of alpha:
    arma::mat alphacov = outer_arma( alpha, c_MUa.t() );
    
    //- compute degrees of freedom and scale matrix:
    double k_alpha;
    arma::mat B_alpha;
    if ( type_proposal == "pmwg" ) {
        k_alpha = nu0 + K - 1 + I;
        B_alpha = 2.0 * nu0 * arma::diagmat( 1.0/a_d ) + alphacov;
    } else {
        k_alpha = nu0 + K + I;
        B_alpha = S0 + alphacov;
    }
    
    //- make the draw:
    arma::mat c_invSIGa = arma::wishrnd( arma::inv_sympd(B_alpha), k_alpha );
    arma::mat c_SIGa    = arma::inv_sympd( c_invSIGa );
    
    //- update a_d:
    if ( type_proposal == "pmwg" ) {
        a_d = ddm_update_a_rcpp( a_d, nu0, c_invSIGa, K );
    }
    
    //- output:
    return Rcpp::List::create(
        Rcpp::Named("c_SIGa")    = c_SIGa,
        Rcpp::Named("c_invSIGa") = c_invSIGa,
        Rcpp::Named("a_d")       = a_d
    );
}