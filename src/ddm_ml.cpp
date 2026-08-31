
#include "ddm_ml.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double ddm4_nllfct_export( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs, const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std", const int& kmax = 5000,
    const double& delta = 1e-29, bool use_lan = false )
{
    
    return -1*ddm4_logLdata_rcpp( alpha, rt, xs, type_alpha, type_ddm, kmax, delta, use_lan );
}

// [[Rcpp::export]]
double ddm4_random_llfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std",
    const int& kmax = 5000, const double& delta = 1e-29,
    bool use_lan = false )
{

    //- get ll value of the data:
    double log_data = ddm4_logLdata_rcpp(alpha, rti, xsi, type_alpha, type_ddm, kmax, delta, use_lan );

    //- get ll value of the prior:
    arma::mat alpha_mat( 1, alpha.size() );
    alpha_mat.row(0) = alpha.t();
    double log_prior = arma::as_scalar( dmvnrm_arma( alpha_mat, MU.t(), SIGMA, true ) );
    
    //- output:
    return log_data + log_prior;

}

// [[Rcpp::export]]
double ddm4_random_nllfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std",
    const int& kmax = 5000, const double& delta = 1e-29,
    bool use_lan = false )
{
    return -1*ddm4_random_llfct_rcpp(alpha, rti, xsi, MU, SIGMA, type_alpha, type_ddm, kmax, delta, use_lan); 
}

//[[Rcpp::export]]
double ddm7_nllfct_export( const arma::colvec& alpha, const arma::vec& rt, 
    const arma::ivec& xs, const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std", const int& kmax = 5000,
    const double& delta = 1e-29, bool use_lan = false )
{
    
    return -1*ddm7_logLdata_rcpp( alpha, rt, xs, type_alpha, type_ddm, kmax, delta, use_lan );
}

// [[Rcpp::export]]
double ddm7_random_llfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std",
    const int& kmax = 5000,
    const double& delta = 1e-29,
    bool use_lan = false )
{

    //- get ll value of the data:
    double log_data = ddm7_logLdata_rcpp(alpha, rti, xsi, type_alpha, type_ddm, kmax, delta, use_lan );

    //- get ll value of the prior:
    arma::mat alpha_mat( 1, alpha.size() );
    alpha_mat.row(0) = alpha.t();
    double log_prior = arma::as_scalar( dmvnrm_arma( alpha_mat, MU.t(), SIGMA, true ) );
    
    //- output:
    return log_data + log_prior;

}

// [[Rcpp::export]]
double ddm7_random_nllfct_rcpp( const arma::colvec& alpha, 
    const arma::vec& rti, const arma::ivec& xsi, 
    const arma::vec& MU, const arma::mat& SIGMA,
    const std::string& type_alpha = "dao",
    const std::string& type_ddm = "std",
    const int& kmax = 5000,
    const double& delta = 1e-29,
    bool use_lan = false )
{
    return -1*ddm7_random_llfct_rcpp(alpha, rti, xsi, MU, SIGMA, type_alpha, type_ddm, kmax, delta, use_lan); 
}

// [[Rcpp::export]]
Rcpp::List ddm_agh_gradfct_singleperson_rcpp( 
    const arma::vec& typevec,
    const arma::mat& posmat,
    const arma::vec& rti, const arma::ivec& xsi,
    const arma::mat& pts, const arma::vec& wgh,
    const arma::vec& MU, 
    const arma::mat& L, 
    const arma::mat& SIGMA,
    const arma::mat& invSIGMA,
    const std::string& type_alpha="dao",
    const std::string& type_ddm="std",
    const std::string& ddm="four",
    const int& kmax = 5000,
    const double& delta = 1e-29,
    bool use_lan = false, std::string lan_path = "" )
{
    
    //- should we load the lan-weights?
    if ( use_lan ) {
        lan_load_weights_rcpp( lan_path );
    }

    // set pointer to correct ddm-function
    typedef double (*DdmFun)(
        const arma::colvec&, const arma::vec&, const arma::ivec&,
        const arma::vec&, const arma::mat&,
        const std::string&, const std::string&, 
        const int&, const double&, bool
    );
    DdmFun ddm_fun = (ddm == "four") 
        ? ddm4_random_llfct_rcpp 
        : ddm7_random_llfct_rcpp;

    //- matrix to save the results:
    int npts = pts.n_rows;
    int lgt_tab = posmat.n_rows;
    arma::vec outl = zeros( npts );
    arma::mat outg = zeros( npts, lgt_tab );

    //- iterate...
    for ( int pp = 0; pp < npts; pp++ ) {
        
        //- get the point:
        arma::rowvec Ppp = pts.row(pp);

        //- compute ll of response data:
        double fpp = ddm_fun( Ppp.t(), rti, xsi, MU, SIGMA, type_alpha, type_ddm, kmax, delta, use_lan );
        double fii = std::exp( fpp )*wgh(pp);
        outl( pp ) = fii;

        //- some precomputations for the gradient:
        arma::colvec resPpp = Ppp.t() - MU;
        arma::rowvec alphainvSIGMA = resPpp.t()*invSIGMA;
        
        //- now we compute the gradient:
        for (int jj=0; jj<lgt_tab; jj++) {
            
            //- get information for derivative
            int type_jj = typevec( jj );
            arma::rowvec pos = posmat.row( jj );
            
            //- compute derivative:
            double der = 0;
            if ( type_jj >= 1 ) {
                der = ddm_derivatives_rcpp( MU, L, invSIGMA, alphainvSIGMA, pos, type_jj ); 
            } else {
                der = alphainvSIGMA( pos( 0, 0 ) - 1 );
            }

            //- save results:
            outg( pp, jj ) = der*fii;
        
        } // jj   
  
    } // pp
  
    // final computations:
    double lli = sum( outl );
    arma::rowvec gradi = sum( outg, 0);
    
    // make list:
    return Rcpp::List::create(
        Rcpp::Named("lli") = lli, 
        Rcpp::Named("gri") = gradi );   
}