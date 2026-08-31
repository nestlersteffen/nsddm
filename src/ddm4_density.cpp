
#include "ddm4_density.h"
#include "ddm_helpers.h"
#include "ddm_lan.h"

using namespace Rcpp;
using namespace arma;

///*********************************************
///**  4-parameter DDM - Tuerlinckx (2004)
///*********************************************

double ddm4_dstand_rcpp( const double& t, const double& a, const double& w, const double& v, 
    const double& s2, const int& kmax, const double& delta ) 
{
  
    // Rcpp::Rcout << "kmax: " << kmax << std::endl;
    // Rcpp::Rcout << "delta: " << delta << std::endl;

    //* compute multiplicative term:
    double mult = (s2*M_PI)/(a*a) * std::exp( ( -w*a*v - 0.5*v*v*t )/s2 );
    
    //* approximate sum:
    double sum_hist[3] = {0.0, 0.0, 0.0}; 
    for ( int k = 1; k < kmax; k++ ){
        
        //* for things to save Note:
        sum_hist[0] = sum_hist[1];
        sum_hist[1] = sum_hist[2];
    
        //* compute elements of the sum:
        double pt_1 = std::sin( M_PI*k*w );
        double pt_2 = std::exp( (-0.5*M_PI*M_PI*k*k*s2*t)/(a*a) );
        
        //* compute the sum and save:
        sum_hist[2] = sum_hist[1] + k*pt_1*pt_2;

        //* check:
        if ( ( std::abs( sum_hist[0] - sum_hist[1] ) < delta ) & 
            ( std::abs( sum_hist[1] - sum_hist[2] ) < delta ) & 
            ( sum_hist[2] > 0 ) ) {
        break;
        }

    } // end k - loop

    //* final check
    if ( sum_hist[2] <= 0 ) {
        return 0;
    }
    
    //* output:
    return sum_hist[2]*mult;
}

///*********************************************
///**  4-parameter DDM - Navarro & Fuss (2009)
///*********************************************

int get_ks( const double& t, const double& err ) //** small-time k
{ 
    // see Navarro & Fuss (2009), p.29 
    if (err * 2 * sqrt(2 * M_PI * t) < 1) { // if error threshold is set low enough
        double ks = 2 + sqrt(-2 * t * log(2 * err * sqrt(2 * M_PI * t)));
        double bc = sqrt(t) + 1; // boundary conditions
        if (ks > MAX || bc > MAX) return MAX;
        return std::ceil( std::max( ks, bc ) ); 
  }
  return 2;
}

int get_kl( const double& t, const double& err ) //** large-time k
{ 
    // see Navarro & Fuss (2009), p.29 
    double bc = 1 / (M_PI * sqrt(t) ); // boundary conditions 
    if ( bc > MAX) return MAX;
    if (err * M_PI * t < 1) { // error threshold is low enough
        double kl = sqrt(-2 * log( M_PI * t * err) / ( M_PI*M_PI*t));
        if ( kl > MAX ) return MAX;
        return std::ceil( std::max( kl, bc ) );
    }
    return std::ceil( bc ); 
}

double ddm4_dnavfuss_rcpp( const double& t, const double& a, const double& w, 
    const double& v, const double& err )
{

    //* adapt reaction time:
    double tt = t/(a*a);

    //* get ks and kl:
    int ks = get_ks( t/(a*a), err );
    int kl = get_kl( t/(a*a), err );

    //* compute multiplicative term:
    double mult = std::exp( -w*a*v - 0.5*v*v*t )/(a*a);

    //* approximate sum:
    double dens = 0;
    if ( ks < kl ) { // small-time is better
        double gamma = -1/(2*tt);
        double tmp_sum = 0.0;
        int k_lower = -(ks-1)/2;
        int k_upper = ks/2;    
        for (int k = k_lower; k <= k_upper; k++) {
            tmp_sum += (w + 2*k) * exp( gamma*(w + 2*k)*(w + 2*k) );
        }
        dens = tmp_sum/sqrt(2*M_PI*tt*tt*tt);
    } else { // large-time is better
        double gamma = (-0.5*M_PI*M_PI*tt);
        double tmp_sum = 0.0;
        for (int k = 1; k <= kl; k++) {
            tmp_sum += sin(M_PI*k*w) * exp( gamma*k*k ) * k;
        }
        dens = tmp_sum*M_PI;
    }
    return dens*mult;
}

arma::vec ddm4_pdf_rcpp( const arma::vec& rt, const arma::ivec& x, 
  const double& a, const double& t0, const double& w, const double& z, 
  const double& v, const std::string& type_ddm, const int& kmax, const double& delta ) 
{
  
    //* how many densities to compute?
    int nn = rt.size();
    arma::vec densities(nn, arma::fill::zeros);

    //* parametrization?
    double w_, z_;
    if( !R_IsNA(w) ) {
        w_ = w;
        z_ = a * w_;
    } else if( !R_IsNA(z) ) {
        z_ = z;
        w_ = z / a;
    } else {
        // Rcpp::stop("Either w or z have to be supplied.");
        return densities;
    }
    
    //* check parameters once (gelten für alle i)
    bool check = ddm4_parmcheck( a, t0, z_, w_ );
    if( !check ) return densities;

    //* let's go:
    for ( int i = 0; i < nn; i++ ) {

        double t = rt[i] - t0;
        if( t <= 0 ) continue;

        if ( type_ddm == "std" ) {    
            if ( x[i] == 0 ) { // response is "lower" or 0 
                densities[i] = ddm4_dstand_rcpp( t, a, w_, v, 1, kmax, delta );
            } else { // response is "upper" or 1
                densities[i] = ddm4_dstand_rcpp( t, a, 1-w_, -v, 1, kmax, delta );
            }
    
        } else if ( type_ddm == "navfuss" ) {
            if ( x[i] == 0 ) { // response is "lower" or 0 
                densities[i] = ddm4_dnavfuss_rcpp( t, a, w_, v, 0.000001 );
            } else { // response is "upper" or 1
                densities[i] = ddm4_dnavfuss_rcpp( t, a, 1-w_, -v, 0.000001 );
            }
        } 
    }

    //* output:
    return densities;
}

double ddm4_logLdata_rcpp( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& x, const std::string& type_alpha, const std::string& type_ddm, 
    const int& kmax, const double& delta, bool use_lan )
{
    
    //* get DDM parameters from alpha:
    double v, a, t0, w, z;
    v  = alpha(0);
    t0 = std::exp( alpha(3) );
    
    if ( type_alpha == "dao" ) {
        a  = std::exp( alpha(1) ) + std::exp( alpha(2) );
        z  = std::exp( alpha(2) );
        w  = z / a;
    } else { // logit-type
        a  = std::exp( alpha(1) );
        w  = 1.0 / ( 1.0 + std::exp( -alpha(2) ) ); // plogis
        z  = a * w;
    }

    //* compute densities:
    arma::vec log_densities;
    if ( use_lan ) {
        log_densities = ddm4_lanll_weights_rcpp( rt, x, a, t0, z, v );
    } else {
        // z immer übergeben, w weglassen (NA_REAL ist Default)
        arma::vec densities = ddm4_pdf_rcpp( rt, x, a, t0, w, NA_REAL, v, type_ddm, kmax, delta );
        densities.replace( 0.0, 1e-29 );
        densities( arma::find_nonfinite(densities) ).fill( 1e-29 );
        log_densities = arma::log( densities );
    }

    //* compute log-likelihood:
    double nll = arma::sum( log_densities );
    return nll;
}