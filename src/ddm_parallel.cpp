
#ifdef _OPENMP
  #include <omp.h>
#endif

#include "ddm_ml.h"
#include "ddm_derivatives.h"
#include "ddm_lan.h"

// [[Rcpp::export]]
Rcpp::List ddm_agh_gradfct_allpersons_rcpp( 
    const arma::vec& typevec, const arma::mat& posmat,
    const arma::vec& rts, const arma::ivec& xs,
    const arma::imat& info, const int I,
    const arma::cube& pts_array, const arma::mat& wgh_mat,
    const arma::vec& MU, 
    const arma::mat& L, 
    const arma::mat& SIGMA,
    const arma::mat& invSIGMA,
    const std::string& type_ddm="std",
    const std::string& type_alpha="dao",
    const std::string& ddm="four",
    const int& kmax = 5000,
    const double& delta = 1e-29,
    bool use_lan = false, std::string lan_path = "",
    int n_threads = 0 )
{
    
    //- should we load the lan-weights?
    if ( use_lan ) {
        lan_load_weights_rcpp( lan_path );
    }

    //- define number of threads
    #ifdef _OPENMP
        int max_threads = omp_get_max_threads();
        int use_threads;
        if ( n_threads > 0 ) {
            use_threads = std::min( n_threads, I );
        } else {
            use_threads = std::min( max_threads, I );
        }
        if ( use_threads < 1 ) use_threads = 1;
        omp_set_num_threads( use_threads );
    #endif

    //- things for the output:
    int lgt_tab = posmat.n_rows;
    arma::vec lli_out( I );
    arma::mat gri_out( I, lgt_tab );

    //- iterate...
    #pragma omp parallel for
    for ( int i = 0; i < I; i++ ) {
        
        //- get the relevant information:
        int ni   = info(i, 1);
        int idx2 = info(i, 2) - 1;       // 0-basiert
        int idx1 = idx2 - ni + 1;
        arma::vec  rti = rts.subvec( idx1, idx2 );
        arma::ivec xsi = xs.subvec( idx1, idx2 );

        // ddm-Funktionspointer (wie im Original):
        typedef double (*DdmFun)(
            const arma::colvec&, const arma::vec&, const arma::ivec&,
            const arma::vec&, const arma::mat&,
            const std::string&, const std::string&, 
            const int&, const double&, bool
        );
        DdmFun ddm_fun = (ddm == "four")
            ? ddm4_random_llfct_rcpp
            : ddm7_random_llfct_rcpp;

        //- get the points:
        const arma::mat pts = pts_array.slice(i);
        const arma::vec wgh = wgh_mat.row(i);
        int npts = pts.n_rows;
        
        //- person-specific components:
        arma::vec outl = arma::zeros( npts );
        arma::mat outg = arma::zeros( npts, lgt_tab );

        // AGH-Quadratur ueber die Punkte:
        for ( int pp = 0; pp < npts; pp++ ) {
            
            arma::rowvec Ppp = pts.row( pp );
            double fpp = ddm_fun( Ppp.t(), rti, xsi, MU, SIGMA,
                                   type_alpha, type_ddm, kmax, 
                                   delta, use_lan );
            double fii = std::exp( fpp ) * wgh( pp );
            outl( pp ) = fii;
            
            arma::colvec resPpp = Ppp.t() - MU;
            arma::rowvec alphainvSIGMA = resPpp.t() * invSIGMA;
            
            for ( int jj = 0; jj < lgt_tab; jj++ ) {
                int type_jj = typevec( jj );
                arma::rowvec pos = posmat.row( jj );
                double der = 0;
                if ( type_jj >= 1 ) {
                    der = ddm_derivatives_rcpp( MU, L, invSIGMA, 
                                                 alphainvSIGMA, pos, type_jj );
                } else {
                    der = alphainvSIGMA( pos(0, 0) - 1 );
                }
                outg( pp, jj ) = der * fii;
            }
        }
        
        // Ergebnisse der Person in thread-lokale Ausgabezeilen schreiben:
        lli_out( i ) = arma::sum( outl );
        gri_out.row( i ) = arma::sum( outg, 0 );

    }

    // make list:
    return Rcpp::List::create(
        Rcpp::Named("lli") = lli_out, 
        Rcpp::Named("gri") = gri_out );   
}