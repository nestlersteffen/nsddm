#include <TMB.hpp>
#include "ddmTMB.h"
using namespace ddmTMB;


template<class Type, class Functor>
matrix<Type> eva_hessian(Functor f, vector<Type> &u) {
    return autodiff::hessian(f, u);
}

// ll_data that we use in EVA

template<class Type>
struct joint_nll {
    
    // Data for DDM:
    vector<Type> rti; // reaction times person i
    vector<Type> xsi; // responses (0/1) person i
    
    // DDM specific parameters:
    int method;       // which approximation to use
    int type_alpha;
    
    // Constructor:
    joint_nll(  vector<Type> rti_,
                vector<Type> xsi_,
                int method_,
                int type_alpha_ ) :
        rti(rti_), xsi(xsi_), method(method_), type_alpha(type_alpha_) {}
    
    // Evaluate the negative joint log-likelihood with random effects
    template <typename T>
    T operator()(vector<T> u) {
        
        T res = 0;

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
                if ( x == 0 ) tmp = ddm4_dstand_tmb( t, a, w, v, 1, 50 );
                else          tmp = ddm4_dstand_tmb( t, a, 1-w, -v, 1, 50 ); //50
            }
            tmp = CppAD::CondExpGt( t, T(0), tmp, T(1e-29) );
            res += log( tmp );
        }
        
        return res;
    }
};

// main function:

template<class Type>
Type objective_function<Type>::operator() () 
{
	using namespace density;

	// -------------
    //  parameters
    PARAMETER_VECTOR(MU);    
    PARAMETER_VECTOR(vecL); 
    PARAMETER_VECTOR(u);
    
    // ------
    //  data
    DATA_VECTOR(rts);       
    DATA_VECTOR(xs);  
    DATA_IMATRIX(info);        // infos about level 2 units [person_id, n_trials]
    DATA_INTEGER(K);           // no. of random effects (should be 4 for a,nu,t0,z)
    DATA_INTEGER(I);           // no. of persons
    DATA_INTEGER(method);
    DATA_INTEGER(type_alpha);
  
    // Construct covariance matrix from Cholesky decomposition:
    // matrix<Type> L(K, K);
    // L.setZero();
    // int k = 0;
    // for(int j = 0; j < K; j++) {
    //     for(int i = j; i < K; i++) {
    //         if(i == j) {
    //             L(i, j) = exp(vecL(k));  // Diagonale immer positiv
    //         } else {
    //             L(i, j) = vecL(k);       // Nebendiagonale frei
    //         }
    //         k++;
    //     }
    // }
    // matrix<Type> Lt = L;
    // matrix<Type> SIGMA = L * L.transpose();
    // matrix<Type> invSIGMA = SIGMA.inverse();
    

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
    matrix<Type> invSIGMA = SIGMA.inverse();

    // ----------------------------------------
    //  log-likelihood-computations

    Type app2 = 0.5*atomic::logdet( invSIGMA );
    
    vector<Type> res( I );
    matrix<Type> check(I, 3);
    for ( int n = 0; n < I; n++ ) {
      	
      	// get the data of person ii:
        int idx1 = info(n, 2) - info(n, 1);
        int idx2 = info(n, 1);
        // int ni = info(n,1);
        vector<Type> rti = rts.segment( idx1, idx2 );
        vector<Type> xsi = xs.segment( idx1,  idx2 );
        
        // get variational parameters:
        vector<Type> ui(K);
        ui.setZero();
        for ( int kk = 0; kk < K; kk++ ) {
          ui(kk) = u( n*K + kk );
        }
        vector<Type> rui = ui - MU;
        matrix<Type> ruitrui = ddmTMB::outer_product( rui );

        // now, we compute the likelihood:
        Type app4 = -0.5*(  ruitrui*invSIGMA ).trace();

        // construct eva approximation:
        joint_nll<Type> jnll( rti, xsi, method, type_alpha );

        // calculate approximation
        Type app13 = ddmTMB::eva( jnll, ui, invSIGMA );
        
        // add all things together:
        res(n) = app13 + app2 + app4;
        check(n,0) = app13;
        check(n,1) = app2;
        check(n,2) = app4; 

        // Zusätzlich für Debugging (nur Person 0):
        if (n == 47) {
            matrix<Type> H = eva_hessian(jnll, ui);
            matrix<Type> invUistar = invSIGMA - H;
            Type fval = jnll(ui);
            REPORT(fval);
            REPORT(ui);
            REPORT(H);
            REPORT(invUistar);
            REPORT(invSIGMA);
        }

    }

    REPORT( res );
    REPORT( check );

    // finally...
    Type nll = -1*res.sum();
    return nll;

}
