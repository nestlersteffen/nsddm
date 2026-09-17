//// File Name: ddm_bayes.cpp
//// File Version: 0.01

// [[Rcpp__interfaces(r, cpp)]]  substitute "__" by "::"
// for exporting Rcpp functions for use in another package

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>    
#include "ddm_helpers.h" 
#include "ddm_bayes_helpers.h" 
#include "ddm_lan.h" 
#include "ddm_mvnorm.h" 
#include "ddm4_density.h"    
#include "ddm7_density.h"
#include "ddm_derivatives.h"

using namespace Rcpp;
using namespace arma;

// define the function to be used:
typedef double (*LogLFun)(
    const arma::colvec&, const arma::vec&, const arma::ivec&,
    const std::string&, const std::string&, const int&, const double&, bool
);

arma::colvec ddm_pmwg_step_rcpp( const arma::colvec& c_alpha, const double c_logData,
    const arma::vec& rti, const arma::ivec& xsi, const int K,
    const arma::colvec& c_MUa, const arma::mat& c_SIGa,
    const arma::colvec& mu_hat, const arma::mat& sigma_hat, const bool use_adapt,
    const int R, const double epsilon, const double pi_mix,
    const double pi_mix1, const double pi_mix2, const double pi_mix3,
    const std::string type_alpha, const std::string method,
    bool& accepted, double& logData_new, LogLFun logL_fun, 
    const int& kmax = 5000, const double& delta = 1e-29, bool use_lan = false )
{
    //- step 1: fix first particle as current alpha:
    arma::mat particles( R, K, arma::fill::zeros );
    particles.row(0) = c_alpha.t();

    //- step 2: draw R-1 new particles from proposal:
    for ( int r = 1; r < R; r++ ) {
        double u = R::runif(0,1);
        if ( use_adapt ) {
            if ( u < pi_mix1 ) {
                particles.row(r) = arma::mvnrnd( mu_hat, sigma_hat, 1 ).t();
            } else if ( u < pi_mix1 + pi_mix2 ) {
                particles.row(r) = arma::mvnrnd( c_alpha, sigma_hat, 1 ).t();
            } else {
                particles.row(r) = arma::mvnrnd( c_MUa, c_SIGa, 1 ).t();
            }
        } else {
            if ( u < pi_mix ) {
                particles.row(r) = arma::mvnrnd( c_alpha, epsilon*c_SIGa, 1 ).t();
            } else {
                particles.row(r) = arma::mvnrnd( c_MUa, c_SIGa, 1 ).t();
            }
        }
    }

    //- step 3: compute log-weights:
    arma::vec log_liks( R, arma::fill::zeros );
    arma::vec log_weights( R, arma::fill::zeros );
    arma::vec log_prior = dmvnrm_arma( particles, c_MUa.t(), c_SIGa, true );

    for ( int r = 0; r < R; r++ ) {

        //- log-likelihood:
        if ( r == 0 ) {
            log_liks(r) = c_logData;
        } else {
            log_liks(r) = logL_fun( particles.row(r).t(), rti, xsi,
                type_alpha, method, kmax, delta, use_lan );
        }

        //- log-proposal:
        arma::mat p_mat = particles.row(r);
        arma::mat c_mat( 1, K );
        c_mat.row(0) = c_alpha.t();
        double log_prop;
        if ( use_adapt ) {
            log_prop = std::log(
                pi_mix1 * arma::as_scalar( dmvnrm_arma( p_mat, mu_hat.t(), sigma_hat, false ) ) +
                pi_mix2 * arma::as_scalar( dmvnrm_arma( p_mat, c_alpha.t(), sigma_hat, false ) ) +
                pi_mix3 * arma::as_scalar( dmvnrm_arma( p_mat, c_MUa.t(), c_SIGa, false ) ) );
        } else {
            log_prop = std::log(
                pi_mix * arma::as_scalar( dmvnrm_arma( p_mat, c_alpha.t(), epsilon*c_SIGa, false ) ) +
                (1-pi_mix) * arma::as_scalar( dmvnrm_arma( p_mat, c_MUa.t(), c_SIGa, false ) ) );
        }
        log_weights(r) = log_liks(r) + log_prior(r) - log_prop;
    }

    //- step 4: normalize weights (log-sum-exp):
    log_weights -= log_weights.max();
    arma::vec weights = arma::exp( log_weights );
    weights /= arma::sum( weights );

    //- step 5: sample new index:
    arma::vec cumweights = arma::cumsum( weights );
    double u = R::runif(0,1);
    int k = 0;
    while ( k < R-1 && u > cumweights(k) ) k++;

    //- output:
    accepted     = ( k != 0 );
    logData_new  = log_liks(k);
    return particles.row(k).t();
}

arma::colvec ddm_standard_step_rcpp( const arma::colvec& c_alpha, const double c_logData,
    const arma::vec& rti, const arma::ivec& xsi, const int K,
    const arma::colvec& c_MUa, const arma::mat& c_SIGa,
    const arma::mat& sigma, const double tau2,
    const double epsilon, const double pi_mix, 
    const std::string type_proposal, const std::string type_alpha, const std::string method,
    bool& accepted, double& logData_new, LogLFun logL_fun,
    const int& kmax = 5000, const double& delta = 1e-29, bool use_lan = false )
{
    //- initialize output:
    accepted    = false;
    logData_new = c_logData;
    arma::colvec alpha_new = c_alpha;

    //- step 1: logPrior for current alpha:
    arma::mat c_mat( 1, K );
    c_mat.row(0) = c_alpha.t();
    double c_logPrior = arma::as_scalar( dmvnrm_arma( c_mat, c_MUa.t(), c_SIGa, true ) );

    //- step 2: generate proposal:
    arma::colvec p_alpha = ddm_generate_alphai_rcpp( c_alpha, sigma, c_MUa, c_SIGa,
        tau2, epsilon, pi_mix, type_proposal );

    //- step 3: logPrior and logData for proposal:
    arma::mat p_mat( 1, K );
    p_mat.row(0) = p_alpha.t();
    double p_logPrior = arma::as_scalar( dmvnrm_arma( p_mat, c_MUa.t(), c_SIGa, true ) );
    double p_logData  = logL_fun( p_alpha, rti, xsi, type_alpha, method, kmax, delta, use_lan );

    //- step 4: compute log-ratio:
    double log_llratio;
    if ( type_proposal == "mixture" ) {
        double log_q_forward = std::log(
            pi_mix * arma::as_scalar( dmvnrm_arma( p_mat, c_alpha.t(), epsilon*c_SIGa, false ) ) +
            (1-pi_mix) * arma::as_scalar( dmvnrm_arma( p_mat, c_MUa.t(), c_SIGa, false ) ) );
        double log_q_backward = std::log(
            pi_mix * arma::as_scalar( dmvnrm_arma( c_mat, p_alpha.t(), epsilon*c_SIGa, false ) ) +
            (1-pi_mix) * arma::as_scalar( dmvnrm_arma( c_mat, c_MUa.t(), c_SIGa, false ) ) );
        log_llratio = ( p_logData + p_logPrior - log_q_forward ) -
                      ( c_logData + c_logPrior - log_q_backward );
    } else {
        log_llratio = ( p_logData + p_logPrior ) - ( c_logData + c_logPrior );
    }

    //- step 5: accept or reject:
    double llratio = std::exp( std::min( log_llratio, 0.0 ) );
    if ( R::runif(0,1) < llratio ) {
        alpha_new   = p_alpha;
        accepted    = true;
        logData_new = p_logData;
    }

    return alpha_new;
}

// [[Rcpp::export]]

Rcpp::List ddm_chain_rcpp( const arma::vec& rts, const arma::ivec& xs,
    const arma::mat& info, const int I, const int K, arma::mat alpha,
    arma::colvec c_MUa, arma::mat c_SIGa, const arma::colvec& m, 
    const arma::mat& invM, const arma::mat& S0, const int nu0, arma::vec a_d,
    const int biter, const int burnin, bool use_adapt_dao,
    const double tau2, const double epsilon, 
    const double pi_mix, const double pi_mix1, const double pi_mix2, const double pi_mix3, const int R,
    const std::string method, const std::string type_alpha, const std::string type_proposal,
    const std::string type_sigma_prior, const std::string ddm, bool verbose, 
    const int& kmax = 5000, const double& delta = 1e-29,
    bool use_lan = false, std::string lan_path = "" )
{

    //- should we load the lan-weights?
    if ( use_lan ) {
        lan_load_weights_rcpp( lan_path );
    }

    //- load the correct ddm - function
    LogLFun logL_fun = (ddm == "four") 
        ? ddm4_logLdata_rcpp 
        : ddm7_logLdata_rcpp;

    // get indices to extract sigma:
    arma::uvec lower_indices = arma::find(arma::trimatl(arma::ones(K,K)));

    //- make matrices
    arma::mat MUa(K,biter);
    MUa.col(0) = c_MUa;
    arma::mat SIGa(K*(K+1)/2,biter);
    SIGa.col(0) = c_SIGa( lower_indices );

    //- array to store alpha draws after burnin for adaptation:
    arma::cube alpha_store;
    arma::mat mu_hats;
    arma::cube sigma_hats;
    arma::uvec has_adapted( I, arma::fill::zeros );
    if ( type_proposal == "pmwg" && use_adapt_dao ) {
        alpha_store.zeros( biter - burnin, I, K );
        mu_hats.zeros( I, K );
        sigma_hats.zeros( K, K, I );
        // initialize sigma_hats as identity matrices:
        for ( int i = 0; i < I; i++ ) {
            sigma_hats.slice(i) = arma::eye( K, K );
        }
    }

    //- pre-step: compute logData of current us:
    arma::vec c_logData( I, arma::fill::zeros );
    int idx = 0;
    for ( int i = 0; i < I; i++ ) {
        // get the data:
        int ni = info(i, 1);
        arma::vec rti = rts.subvec(idx, idx + ni - 1);
        arma::ivec xsi = xs.subvec(idx, idx + ni - 1);
        // compute logData - values
        c_logData(i) = logL_fun( alpha.row(i).t(), rti, xsi, 
            type_alpha, method, kmax, delta, use_lan );
        // increase idx:
        idx += ni;
    }

    // let's go:
    int accept = 0;
    int accept_post = 0;
    arma::mat c_invSIGa = arma::inv_sympd( c_SIGa );

    for ( int nn = 0; nn < biter; nn++ ) {

        if ( verbose & (nn > 0) & (nn % 10 == 0)) {
            Rcpp::Rcout << "Iteration: " << nn << std::endl;
        }
        
        // * update MUalpha:
        arma::mat covMUa = arma::inv_sympd( invM + I*c_invSIGa );
        arma::rowvec alphasum = arma::sum( alpha, 0 );
        arma::colvec eMUa = covMUa*( invM*m + c_invSIGa*alphasum.t() );
        arma::colvec c_MUa = arma::mvnrnd( eMUa, covMUa, 1 );
        MUa.col(nn) = c_MUa;

        // %% update SIGMAalpha    
        Rcpp::List tmp_SIGa = ddm_update_sigmaalpha_rcpp( alpha, c_MUa, I, K,
            nu0, S0, a_d, type_sigma_prior );
        c_invSIGa = Rcpp::as<arma::mat>( tmp_SIGa["c_invSIGa"] );
        c_SIGa = Rcpp::as<arma::mat>( tmp_SIGa["c_SIGa"] );
        a_d = Rcpp::as<arma::vec>( tmp_SIGa["a_d"] );
        SIGa.col(nn) = c_SIGa( lower_indices );

        // -----------------
        // * update alpha

        //- check wether to use adaptive algo and whether in adaptation phase:
        int n_adapt = 0;
        bool in_adapt = false;
        if ( use_adapt_dao ) {
            n_adapt  = nn - burnin;
            in_adapt = n_adapt > 0;
        }
        
        idx = 0;
        for ( int i = 0; i < I; i++ ) {
            
            // get the data:
            int ni = info(i, 1);
            arma::vec rti = rts.subvec(idx, idx + ni - 1);
            arma::ivec xsi = xs.subvec(idx, idx + ni - 1);
            
            // make the step:
            bool accepted;
            double logData_new;
            arma::colvec alpha_new( K, arma::fill::zeros );
            if ( type_proposal == "pmwg" ) {

                // compute info for mixture:
                arma::colvec mu_hat( K, arma::fill::zeros );
                arma::mat sigma_hat( K, K, arma::fill::eye );
                
                // relevant in case adaptation is used:
                bool use_adapt = false;

                if ( use_adapt_dao && in_adapt && n_adapt > K + 1 && n_adapt % 20 == 0 && n_adapt < 5000 ) {
                    
                    //- extract draws for person i: slice is (n_adapt x K)
                    arma::mat draws_i = alpha_store.tube( 
                        arma::span(0, n_adapt-1), arma::span(i, i) );
                    draws_i.reshape( n_adapt, K );
                    arma::colvec mu_tmp    = arma::mean( draws_i, 0 ).t();
                    arma::mat    sigma_tmp = arma::cov( draws_i );
                    
                    //- make cholesky decomposition test:
                    arma::mat chol_test;
                    if ( arma::chol( chol_test, sigma_tmp ) ) {
                        mu_hats.row(i)      = mu_tmp.t(); 
                        sigma_hats.slice(i) = sigma_tmp;  
                        has_adapted(i)      = 1;
                    }
                    // is test fails, we save nothing
                }

                if ( use_adapt_dao && has_adapted(i) ) {  
                    mu_hat    = mu_hats.row(i).t();
                    sigma_hat = sigma_hats.slice(i);
                    use_adapt = true;
                }
                
                // get new alpha
                alpha_new = ddm_pmwg_step_rcpp( 
                    alpha.row(i).t(), c_logData(i), rti, xsi, K, 
                    c_MUa, c_SIGa, 
                    mu_hat, sigma_hat, use_adapt,
                    R, epsilon, pi_mix, pi_mix1, pi_mix2, pi_mix3,
                    type_alpha, method, accepted, logData_new, 
                    logL_fun, kmax, delta, use_lan );
            
                // store draw for adaptation after burnin:
                if ( use_adapt_dao && in_adapt ) {
                    for ( int k = 0; k < K; k++ ) {
                        alpha_store( n_adapt - 1, i, k ) = alpha_new(k);
                    }
                }

            } else {
                alpha_new = ddm_standard_step_rcpp( alpha.row(i).t(), 
                    c_logData(i), rti, xsi, K, c_MUa, c_SIGa, 
                    arma::eye(K, K),tau2,epsilon, 
                    pi_mix, type_proposal, type_alpha, method, accepted, 
                    logData_new, logL_fun, kmax, delta, use_lan );
            }

            // change everything:
            alpha.row(i) = alpha_new.t();
            c_logData(i) = logData_new;
            if ( accepted ) {
                accept += 1;
                if ( nn > burnin ) {
                    accept_post += 1;
                }
            }

            // increase idx:
            idx += ni;

        }
        // -----------

    } // biter

    return Rcpp::List::create(
        Rcpp::Named("MUa") = MUa.t(),
        Rcpp::Named("SIGa") = SIGa.t(),
        Rcpp::Named("alpha") = alpha,
        Rcpp::Named("accept") = accept,
        Rcpp::Named("accept_post") = accept_post
    );  

}