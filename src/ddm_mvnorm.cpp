
#include "ddm_mvnorm.h"

using namespace Rcpp;
using namespace arma;

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat)
{
    arma::uword const n = trimat.n_cols;
    for(unsigned j = n; j-- > 0;) {
        double tmp(0.);
        for(unsigned i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x[i];
        x[j] = tmp;
    }
}

arma::vec dmvnrm_arma( arma::mat x,  
                       arma::rowvec mean,  
                       arma::mat sigma, 
                       bool use_log ) 
{ 
    using arma::uword;
    uword const n = x.n_rows, 
    xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;

    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
     }  
    if (use_log)
        return out;
    return exp(out);
}