
#include "ddm7_density.h"
#include "ddm4_density.h"
#include "ddm_helpers.h"
#include "ddm_lan.h"

using namespace Rcpp;
using namespace arma;

double ddm7_large_integrand_rcpp( const double& tx, 
    const double& a, const double& z, const double& v, 
    const double& t0, const double& sv, const int& kmax, 
    const double& delta ) 
{
  
    //* compute true reaction time:
    double td = tx - t0;

    //* compute multiplicative term:
    double tdsv2 = 1 + td*sv*sv;
    double log_mult1 = log( M_PI) - 2*log(a) - 0.5*log( tdsv2 );
    double log_mult2 = -0.5*( v*v*td + 2*v*z - z*z*sv*sv )/tdsv2; 
    
    //* approximate sum:
    double sum_hist[3] = {0.0, 0.0, 0.0}; 
    for ( int k = 1; k < kmax; k++ ){
        
        //* for things to save Note:
        sum_hist[0] = sum_hist[1];
        sum_hist[1] = sum_hist[2];
    
        //* compute elements of the sum:
        double pt_1 = std::sin( ( M_PI*k*z )/a );
        double log_pt_2 = (-0.5*M_PI*M_PI*k*k*td)/(a*a);

        //* compute log_hist:
        double log_hist = log( k ) + log_pt_2;

        //* compute the sum and save:
        sum_hist[2] = sum_hist[1] + pt_1*std::exp( log_hist );

        //* check:
        if ( ( std::abs( sum_hist[0] - sum_hist[1] ) < delta ) && 
             ( std::abs( sum_hist[1] - sum_hist[2] ) < delta ) && 
             ( sum_hist[2] > 0 ) ) {
        break;
        }

    } // end k - loop

    //* final check
    if ( sum_hist[2] <= 0 ) {
        return 0;
    }
    
    //* output:
    return sum_hist[2]*exp( log_mult1 )*exp( log_mult2);
}

double ddm7_small_integrand_rcpp( const double& tx,
    const double& a, const double& z, const double& v,
    const double& t0, const double& sv, const int& kmax,
    const double& delta )
{
    //* compute true reaction time:
    double td = tx - t0;

    //* multiplicative terms:
    double sv2 = sv * sv;
    double A   = 1.0 + td * sv2;
    double mult1 = a * std::pow( td, -1.5 ) / std::sqrt( 2.0 * M_PI * A );

    //* log_mult2 - numerisch stabil:
    double log_mult2;
    if ( sv < 1e-6 ) {
        log_mult2 = -0.5 * ( v*v*td + 2.0*v*z );
    } else {
        double B = v - z * sv2;
        log_mult2 = -0.5 * v*v / sv2 + 0.5 * B*B / ( sv2 * A );
    }
    double mult2 = std::exp( log_mult2 );

    //* infinite sum over k = 0, +1, -1, +2, -2, ...:
    auto term_in_the_sum = [&]( int k ) -> double {
        double coeff = z / a + 2.0 * k;
        return coeff * std::exp( -( z + 2.0*k*a )*( z + 2.0*k*a ) / ( 2.0*td ) );
    };

    double sum_hist[3] = {0.0, 0.0, 0.0};
    sum_hist[2] = term_in_the_sum( 0 );

    for ( int m = 1; m < kmax; m++ ) {
        sum_hist[0] = sum_hist[1];
        sum_hist[1] = sum_hist[2];
        sum_hist[2] = sum_hist[1] + term_in_the_sum( m ) + term_in_the_sum( -m );

        if ( m > 1 &&
             std::abs( sum_hist[0] - sum_hist[1] ) < delta &&
             std::abs( sum_hist[1] - sum_hist[2] ) < delta ) break;
    }

    double result = mult1 * mult2 * sum_hist[2];
    if ( !std::isfinite( result ) || result <= 0.0 ) return 0.0;
    return result;
}

double ddm7_dao_rcpp( const double& tx,
    const double& a, const double& z, const double& v,
    const double& t0, const double& sv, const double& sz,
    const double& st0, const int& kmax, const double& err )
{
    //* GL points and weights:
    double pts[10] = {-0.973906528517172,-0.865063366688985,
        -0.679409568299024,-0.433395394129247,-0.148874338981631,
         0.148874338981631, 0.433395394129247, 0.679409568299024,
         0.865063366688985, 0.973906528517172};
    double wgh[10] = {0.066671344308688, 0.149451349150581,
        0.219086362515982, 0.269266719309996, 0.295524224714753,
        0.295524224714753, 0.269266719309996, 0.219086362515982,
        0.149451349150581, 0.066671344308688};

    double total = 0.0;
    for ( int i = 0; i < 10; i++ ) {
        for ( int j = 0; j < 10; j++ ) {

            //* get current quadrature points:
            double z_ij  = z  + 0.5*sz  * pts[i];
            double t0_ij = t0 + 0.5*st0 * pts[j];

            //* check conditions:
            if ( z_ij  <= 0.0 || z_ij >= a ) continue;
            if ( t0_ij <= 0.0 ) continue;
            if ( tx - t0_ij <= 0.0 ) continue;

            //* Navarro & Fuss: get ks and kl for this quadrature point:
            double td_ij = tx - t0_ij;
            double tt_ij = td_ij / ( a*a );
            int ks = get_ks( tt_ij, err );
            int kl = get_kl( tt_ij, err );

            //* choose integrand:
            double f_ij;
            if ( ks < kl ) {
                f_ij = ddm7_small_integrand_rcpp( tx, a, z_ij, v, t0_ij, sv, ks, 1e-29 );
            } else {
                f_ij = ddm7_large_integrand_rcpp( tx, a, z_ij, v, t0_ij, sv, kl, 1e-29 );
            }

            total += wgh[i] * wgh[j] * f_ij;
        }
    }

    return 0.25 * total;
}

double ddm7_dstand_rcpp( const double& tx, const int& x, 
    const double& a, const double& t0, const double& z, 
    const double& v, const double& st0, const double& sz, 
    const double& sv, const double& s2, const int& kmax, 
    const double& delta, const double& epsilon, 
    const double& min_RT ) 
{
    //* points and weights for Gauss-Hermite (u-integration):
    double pts_us[10] = {-3.43615911883774,-2.53273167423279,
        -1.75668364929988,-1.03661082978951,-0.34290132722370,
         0.34290132722370, 1.03661082978951, 1.75668364929988,
         2.53273167423279, 3.43615911883774};
    double wgh_us[10] = {0.00000764043286, 0.00134364574678,
        0.03387439445548, 0.24013861108230, 0.61086263373530,
        0.61086263373530, 0.24013861108230, 0.03387439445548,
        0.00134364574678, 0.00000764043286};

    //* points and weights for Gauss-Legendre (z-integration):
    double pts_zs[10] = {-0.973906528517172,-0.865063366688985,
        -0.679409568299024,-0.433395394129247,-0.148874338981631,
         0.148874338981631, 0.433395394129247, 0.679409568299024,
         0.865063366688985, 0.973906528517172};
    double wgh_zs[10] = {0.066671344308688, 0.149451349150581,
        0.219086362515982, 0.269266719309996, 0.295524224714753,
        0.295524224714753, 0.269266719309996, 0.219086362515982,
        0.149451349150581, 0.066671344308688};

    // Rcpp::Rcout << "kmax: " << kmax << std::endl;
    // Rcpp::Rcout << "delta: " << delta << std::endl;

    //* collect parameters:
    double a2 = a*a;
    double Zu = z + 0.5*sz;
    double Zl = z - 0.5*sz;
    double Tl = t0 - 0.5*st0;
    double Tu = std::min( tx, t0 + 0.5*st0 );

    //* prepare points and transform them:
    double us[10];
    double w_us[10];
    double zs[10];
    double w_zs[10];
    int nr_us = 10;
    int nr_zs = 10;
    double sqrt2  = std::sqrt(2.0);
    double sqrtpi = std::sqrt(M_PI);
    for ( int i = 0; i < 10; i++ ) {
        us[i]   = sqrt2*pts_us[i]*sv + v;
        w_us[i] = wgh_us[i]/sqrtpi;
        zs[i]   = pts_zs[i]*0.5*sz + z;
        w_zs[i] = wgh_zs[i]*0.5*sz;
    }

    //* change variable:
    if ( x == 1 ) {
        Zu = ( a - z ) + 0.5*sz;
        Zl = ( a - z ) - 0.5*sz;
        for ( int i = 0; i < 10; i++ ) {
            us[i] = -us[i];
        }
    }

    //* compute constant things:
    double tTu = -0.5*(tx - Tu);
    double tTl = -0.5*(tx - Tl);

    //* early exit:
    if ( ( tx - t0 + 0.5*st0 ) <= min_RT ) {
        return 0.0;
    }

    double out = 0.0;

    //* standard case:
    if ( tx > t0 + 0.5*st0 ) {

        double sum_hist[3] = {0.0, 0.0, 0.0};
        for ( int k = 1; k < kmax; k++ ) {

            //* for things to save:
            sum_hist[0] = sum_hist[1];
            sum_hist[1] = sum_hist[2];
            double sum_us = 0.0;

            //* some further constant things:
            double pika  = (M_PI*k)/a;
            double pikas = (M_PI*M_PI*k*k*s2)/a2;

            //* gauss quadrature over u:
            for ( int uu = 0; uu < nr_us; uu++ ) {
                //* some tmps:
                double u2s2pikas = us[uu]*us[uu]/s2 + pikas;
                //* start computations:
                double t1   = 1.0/(u2s2pikas*u2s2pikas);
                double t2_1 = std::exp(-us[uu]*Zu/s2)*
                    (-us[uu]*std::sin(pika*Zu)/s2 - pika*std::cos(pika*Zu));
                double t2_2 = std::exp(-us[uu]*Zl/s2)*
                    (-us[uu]*std::sin(pika*Zl)/s2 - pika*std::cos(pika*Zl));
                double t3   = std::exp(u2s2pikas*tTu) - std::exp(u2s2pikas*tTl);
                sum_us += t1*(t2_1 - t2_2)*t3*w_us[uu];
            } //* end uu - loop

            sum_hist[2] = sum_hist[1] + k*sum_us;

            if ( ( std::abs( sum_hist[0] - sum_hist[1] ) < delta ) &&
                 ( std::abs( sum_hist[1] - sum_hist[2] ) < delta ) &&
                 ( sum_hist[2] > 0.0 ) ) {
                break;
            }
        } //* end k - loop

        out = 2.0*sum_hist[2]*s2*s2*M_PI/(a2*sz*st0);

    } else { //* tx <= t0 + st0/2

        double sum_us = 0.0;
        for ( int uu = 0; uu < nr_us; uu++ ) {

            double sum_zs = 0.0;

            if ( std::abs(us[uu]) > epsilon ) {

                for ( int zz = 0; zz < nr_zs; zz++ ) {

                    //* approximate first sum (see Equation 22):
                    double zzz1 = (a - zs[zz])*x + zs[zz]*(1 - x);
                    double zzz2 = (a - zs[zz])*(1 - x) + zs[zz]*x;
                    double sum1 = (a2/(M_PI*s2))*
                        std::sinh(zzz2*us[uu]/s2)/std::sinh((us[uu]*a)/s2);

                    //* get the second sum:
                    double sum_hist[3] = {0.0, 0.0, 0.0};
                    for ( int k = 1; k < kmax; k++ ) {
                        sum_hist[0] = sum_hist[1];
                        sum_hist[1] = sum_hist[2];
                        double pika = (M_PI*k)/a;
                        double u2s2pikas = us[uu]*us[uu]/s2 + (M_PI*M_PI*k*k*s2)/a2;
                        //* tmp computations:
                        double t1 = 1.0/u2s2pikas;
                        double t2 = std::sin(pika*zzz1);
                        double t3 = std::exp(u2s2pikas*tTl);
                        sum_hist[2] = sum_hist[1] + k*t1*t2*t3;
                        if ( ( std::abs( sum_hist[0] - sum_hist[1] ) < delta ) &&
                             ( std::abs( sum_hist[1] - sum_hist[2] ) < delta ) &&
                             ( sum_hist[2] > 0.0 ) ) {
                            break;
                        }
                    } //* end k - loop

                    sum_zs += w_zs[zz]/sz*(sum1 - 2.0*sum_hist[2])*
                        (M_PI*s2)/(a2*st0)*std::exp(-zzz1*us[uu]/s2);

                } //* end zs - loop

            } else { //* abs(us[uu]) < epsilon: degenerate v ~ 0 case

                double sum_hist[3] = {0.0, 0.0, 0.0};
                double sl = -(Zl*M_PI*M_PI)/(2.0*a) + (Zl*Zl*M_PI*M_PI)/(4.0*a2);
                double su = -(Zu*M_PI*M_PI)/(2.0*a) + (Zu*Zu*M_PI*M_PI)/(4.0*a2);
                for ( int k = 1; k < kmax; k++ ) {
                    sum_hist[0] = sum_hist[1];
                    sum_hist[1] = sum_hist[2];
                    double pika  = (M_PI*k)/a;
                    double pika2 = (M_PI*M_PI*k*k*s2)/a2;
                    //* tmp computations:
                    double t1 = std::cos(pika*Zl);
                    double t2 = std::cos(pika*Zu);
                    double t3 = std::exp(pika2*tTl);
                    sum_hist[2] = sum_hist[1] + (1.0/((double)k*k))*(t1 - t2)*t3;
                    if ( ( std::abs( sum_hist[0] - sum_hist[1] ) < delta ) &&
                         ( std::abs( sum_hist[1] - sum_hist[2] ) < delta ) &&
                         ( sum_hist[2] > 0.0 ) ) {
                        break;
                    }
                } //* end k - loop

                sum_zs = ((2.0*a)/(st0*sz*s2*M_PI*M_PI))*(sl - su - sum_hist[2]);

            } //* end if/else on us[uu]

            sum_us += sum_zs*w_us[uu];

        } //* end us - loop

        out = sum_us;
    }

    if ( out <= 0.0 || !std::isfinite( out ) ) {
        return 0.0;
    }

    //* output:
    return out;
}

arma::vec ddm7_pdf_rcpp( const arma::vec& rt, const arma::ivec& x, 
  const double& a, const double& t0, const double& z, const double& v, 
  const double& sv, const double& sz, const double& st0, 
  const std::string& type_ddm, const int& kmax, const double& delta ) 
{
  
    //* how many densities to compute?
    int nn = rt.size();
    arma::vec densities(nn, arma::fill::zeros);

    //* check parameters once (gelten für alle i)
    bool check = ddm7_parmcheck( a, t0, z, sv, st0, sz );
    if( !check ) return densities;

    //* let's go:
    for ( int i = 0; i < nn; i++ ) {

        if ( type_ddm == "std" ) {
            densities[i] = ddm7_dstand_rcpp( rt[i], x[i], a, t0, z, v, st0, sz, sv, 1, kmax, delta, 
                1e-7, 0.001 );
        } else if ( type_ddm == "dao" ) {
            if ( x[i] == 0 ) { // response is "lower" or 0 
                //densities[i] = ddm7_large_density_rcpp( rt[i], a, z, v, t0, sv, sz, st0, kmax, delta );
                densities[i] = ddm7_dao_rcpp( rt[i], a, z, v, t0, sv, sz, st0, kmax, 1e-6 );
            } else { // response is "upper" or 1
                //densities[i] = ddm7_large_density_rcpp( rt[i], a, a-z, -v, t0, sv, sz, st0, kmax, delta );
                densities[i] = ddm7_dao_rcpp( rt[i], a, a-z, -v, t0, sv, sz, st0, kmax, 1e-6 );
            }
        }
    }

    //* output:
    return densities;
}

double ddm7_logLdata_rcpp( const arma::colvec& alpha, 
    const arma::vec& rt, const arma::ivec& x, 
    const std::string& type_alpha,  // aktuell noch ungenutzt, aber für einheitliche Signatur
    const std::string& type_ddm,      // aktuell noch ungenutzt, aber für einheitliche Signatur
    const int& kmax, const double& delta,
    bool use_lan )
{
    
    //* get DDM parameters from alpha:
    double v, a, t0, z, sv, sz, st0;
    v   = alpha(0);
    sv  = std::exp( alpha(4) );
    sz  = std::exp( alpha(5) );
    st0 = std::exp( alpha(6) );
    z   = std::exp( alpha(2) ) + 0.5*sz;
    t0  = std::exp( alpha(3) ) + 0.5*st0;
    a   = std::exp( alpha(1) ) + z + 0.5*sz;

    //* compute densities:
    arma::vec log_densities;
    if ( use_lan ) {
        // if ( use_tf ) {
        //     Rcpp::Rcout << "Dnn" << std::endl;
        //     // log_densities = ddm7_lanll_dnn_rcpp( rt, x, a, t0, z, v, dnn ); 
        // } else {
            Rcpp::Rcout << "Lan" << std::endl;
            log_densities = ddm7_lanll_weights_rcpp( rt, x, a, t0, z, v, sv, sz, st0 );
        // }
    } else {
        // z immer übergeben, w weglassen (NA_REAL ist Default)
        arma::vec densities = ddm7_pdf_rcpp( rt, x, a, t0, z, v, sv, sz, st0, type_ddm, kmax, delta ) ;
        densities( arma::find( densities <= 0.0 ) ).fill( 1e-29 );
        densities( arma::find_nonfinite(densities) ).fill( 1e-29 );
        log_densities = arma::log( densities );
    }

    //* compute log-likelihood:
    double nll = arma::sum( log_densities );
    return nll;
}