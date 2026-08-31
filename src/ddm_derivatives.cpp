
#include "ddm_derivatives.h"  

using namespace Rcpp;
using namespace arma;

double ddm_derivatives_rcpp( const arma::vec& MU, 
    const arma::mat& L, const arma::mat& invSIGMA, 
    const arma::rowvec& alphainvSIGMA, const arma::rowvec& pos, 
    const int type_jj )
{
    // set derivative to zero:
    double tmp = 0;
    
    // make a string vector: 
    if ( type_jj == 0) { // "MU" 
        // 1_ij - matrix:
        int p = MU.n_rows;
        arma::colvec ONE_ij = zeros( p ); 
        ONE_ij( pos(0,0) - 1, 0 ) = 1;
        // derivative: 
        arma::mat res = alphainvSIGMA * ONE_ij;
        tmp = res(0,0); 
    }
    
    if ( type_jj == 1 ) { // "L"
        // 1_ij - matrix
        int p = invSIGMA.n_rows;
        arma::mat ONE_ij = zeros( p, p );
        ONE_ij( pos(0,0) - 1 , pos(0,1) - 1 ) = 1;
        // derivative: 
        arma::mat dSIGMA = ONE_ij*L.t() + L*ONE_ij.t();
        // true ll derivative:
        double pt1 = -0.5*trace( invSIGMA * dSIGMA );
        arma::mat pt2 = 0.5*( alphainvSIGMA * dSIGMA * alphainvSIGMA.t() );
        tmp = pt1 + pt2(0,0); 
    }
    
    if ( type_jj == 2 ) { // "SIGMA"
        // 1_ij - matrix
        int p = invSIGMA.n_rows;
        arma::mat ONE_ij = zeros( p, p );
        ONE_ij( pos(0,0) - 1 , pos(0,1) - 1 ) = 1;
        // derivative: 
        arma::mat dSIGMA = ONE_ij + ONE_ij.t() - ONE_ij*ONE_ij;
        // true ll derivative:
        double pt1 = -0.5*trace( invSIGMA * dSIGMA );
        arma::mat pt2 = 0.5*( alphainvSIGMA * dSIGMA * alphainvSIGMA.t() );
        tmp = pt1 + pt2(0,0); 
    }
    
    // output:
    return tmp;
} 