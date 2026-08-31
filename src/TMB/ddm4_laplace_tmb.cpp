
#include <TMB.hpp>
#include "ddmTMB.h"
using namespace ddmTMB;

// Joint negative log-likelihood for DDM with absolute starting point
template<class Type>
struct joint_nll_ddm {
    
    // Data for DDM:
    vector<Type> rti;     // reaction times
    vector<Type> xsi;      // responses (0/1)
    
    // Fixed effects (population means):
    vector<Type> MU;      // [a, nu, t0, z]
    matrix<Type> SIGMA;   // covariance matrix for random effects
    
    // DDM specific parameters:
    int method;
    int type_alpha;
    
    // Constructor:
    joint_nll_ddm(vector<Type> rti_,
                  vector<Type> xsi_,
                  vector<Type> MU_,
                  matrix<Type> SIGMA_,
                  int method_,
                  int type_alpha_) :
        rti(rti_), xsi(xsi_), MU(MU_), SIGMA(SIGMA_), 
        method(method_), type_alpha(type_alpha_) {}
    
    // Evaluate the negative joint log-likelihood
    template <typename T>
    T operator()(vector<T> u) {
        
        T res = T( 0 );
        
        // Cast matrices for template compatibility:
        matrix<T> tSIGMA = SIGMA.template cast<T>();
        vector<T> tMU = MU.template cast<T>();
        
        // Transform random effects to DDM parameters
        // vector<T> total_effects = u + MU.template cast<T>();
        
        // Extract parameters - assumed order: [v, a, z, t0 ]
        T v  = u(0);
        T t0 = exp( u(3) );
        T a  = exp( u(1) );
        T z  = 0; 
        T w  = 0;
        if ( type_alpha == 0 ) { // dao version
            z = exp( u(2) );
            a += z;
            w = z/a; 
        } else {
            w = invlogit( u(2) );
            z = w*a; 
        }
        
        // Random effects contribution (multivariate normal prior):
        density::MVNORM_t<T> neg_log_density(tSIGMA);
        vector<T> eu = u - tMU;
        res = neg_log_density( eu );
        
        // DDM likelihood for all trials of this person:
        for( int i = 0; i < rti.size(); i++ ) {
            
            // get data:
            T rt = T( rti(i) );
            T x = T( xsi(i) );
            
            // change reaction time to process-based component
            T t = rt - t0;

            // now, compute the density
            T tmp = T( 0 );
            if ( method == 0 ) {  // Navarro & Fuss
                if ( x == 0 ) tmp = ddm4_dnavfuss_tmb( t, a, w, v, 0.000001 );
                else          tmp = ddm4_dnavfuss_tmb( t, a, 1-w, -v, 0.000001 );
            } else if ( method == 1) {  // standard
                if ( x == 0 ) tmp = ddm4_dstand_tmb( t, a, w, v, 1, 100 ); // 50
                else          tmp = ddm4_dstand_tmb( t, a, 1-w, -v, 1, 100 );
            }
            tmp = CppAD::CondExpGt( t, T(0), tmp, T(1e-29) );
            res -= log( tmp );
        }
        
        // output
        return res;
    }
};

template<class Type>
Type objective_function<Type>::operator() ()
{
    
    using namespace density;

    // Data:
    DATA_VECTOR(rts);       
    DATA_VECTOR(xs);  
    DATA_MATRIX(ustart);      
    DATA_IMATRIX(info);        // infos about level 2 units [person_id, n_trials]
    DATA_INTEGER(K);           // no. of random effects (should be 4 for a,nu,t0,z)
    DATA_INTEGER(N);           // no. of persons
    DATA_INTEGER(method);
    DATA_INTEGER(type_alpha);
    DATA_INTEGER(niter);       // no. of Newton iterations
    
    // Parameters:
    PARAMETER_VECTOR(MU);    // Fixed effects: [log(a), nu, log(t0), logit(z)]
    PARAMETER_VECTOR(vecL);  // Lower triangular elements of Cholesky decomposition
    
    matrix<Type> L(K, K);
    L.setZero();

    int k = 0;
    for(int j = 0; j < K; j++) {
        for(int i = j; i < K; i++) { 
            L(i, j) = vecL(k); 
            k++;
        }
    }
    matrix<Type> SIGMA = L * L.transpose();
    
    matrix<Type> usave(N,K);
    usave.setZero();
    vector<Type> res(N); 
    
    // Iterate over persons:
    Type nll = 0;
    for(int n = 0; n < N; n++) {
        // PARALLEL_REGION {
        
        // get the data of person ii:
        int idx1 = info(n, 2) - info(n, 1);
        int idx2 = info(n, 1);
        vector<Type> rti = rts.segment( idx1, idx2 );
        vector<Type> xsi = xs.segment( idx1, idx2 );
        
        // Construct joint negative log-likelihood for person ii:
        joint_nll_ddm<Type> jnll( rti, xsi, MU, SIGMA, method, type_alpha );
        
        // Initialize random effects:
        vector<Type> u(K);
        u = ustart.row(n);
        
        // Calculate Laplace approximation:
        Type tmp = laplace(jnll, u, niter);
        res(n) = tmp;
        usave.row(n) = u;
        nll += tmp;
        
    }
    
    REPORT( res );
    REPORT( usave );

    return nll;
}
