/// @file ddm4_tmb.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ddm4_tmb( objective_function<Type>* obj )
{
    // data:
    DATA_VECTOR(rts);
    DATA_IVECTOR(xs);
    DATA_INTEGER(method);  // 0 = navfuss, 1 = standard
    DATA_INTEGER(type_alpha);
    int n = rts.size();
    
    // parameters:
    PARAMETER_VECTOR(parm);
    
    // transform parms depending on type_alpha:
    Type v  = parm(0); 
    Type t0 = exp( parm(3) );
    Type a  = exp( parm(1) );
    Type w = Type(0); Type z = Type(0);
    if ( type_alpha == 0 ) { // dao version
        z = exp( parm(2) );
        a += z;
        w = z/a; 
    } else {
        w = invlogit( parm(2) );
        z = w*a; 
    }
    
    // initialize nll:
    Type nll = Type(0);

    for ( int i = 0; i < n; i++ ) {
        
        // get single trial reaction time and response
        Type rt = rts(i);
        Type x = xs(i);
        
        // change reaction time to process-based component
        Type t = rt - t0;

        // now, compute the density
        Type tmp = Type(0);
        if ( method == 0 ) {  // Navarro & Fuss
            if ( x == 0 ) tmp = ddmTMB::ddm4_dnavfuss_tmb( t, a, w, v, 0.000001 );
            else          tmp = ddmTMB::ddm4_dnavfuss_tmb( t, a, 1-w, -v, 0.000001 );
        } else if ( method == 1) {  // standard
            if ( x == 0 ) tmp = ddmTMB::ddm4_dstand_tmb( t, a, w, v, 1, 100 );
            else          tmp = ddmTMB::ddm4_dstand_tmb( t, a, 1-w, -v, 1, 100 );
        }
        tmp = CppAD::CondExpGt( t, Type(0), tmp, Type(1e-29) );
        nll -= log( tmp );
    }

    ADREPORT( a );
    ADREPORT( v );
    ADREPORT( z );
    ADREPORT( w );
    ADREPORT( t0 );
    
    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this