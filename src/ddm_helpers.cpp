
#include "ddm_helpers.h"    

using namespace Rcpp;
using namespace arma;

bool ddm4_parmcheck( const double& a, const double& t0, const double& z, const double& w )  
{
    bool out = TRUE;
      
    if ( z < 0 || z > a || a <= 0 || t0 <= 0 ) {
        out = FALSE;
    }
    if ( w <= 0 || w >= 1 ) {
        out = FALSE;
    }
    if ( t0 <= 0 ) {
        out = FALSE;
    }
    return out;
}

bool ddm7_parmcheck( const double& a, const double& t0, const double& z, 
    const double& sv, const double& st0, const double& sz )
{
    bool out = TRUE;
       
    if ( a <= 0 || sz > a || z + 0.5*sz > a ) {
        out = FALSE;
    }
    if ( z < 0 || sz < 0 || z - 0.5*sz < 0 || sz > 2 ) {
        out = FALSE;
    }
    if ( t0 < 0 || st0 < 0 || st0 > 1 ) {
        out = FALSE;
    }
    if ( sv < 0 || sv > 3 ) {
        out = FALSE;
    }
    return out; 
}

double ddm_sumexplog_rcpp( Rcpp::NumericVector x) 
{
    int n = x.size();
    //- get maximimum:
    double xmax = max( x );
    //- compute sum:
    double xsum = 0;
    for( int i = 0; i < n; ++i) {
        xsum += std::exp( x(i) - xmax );
    }
    //- result:
    double out = xmax + std::log( xsum );
    return out;  
}