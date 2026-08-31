/// @file ddm7_tmb.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ddm7_tmb( objective_function<Type>* obj )
{
    // data:
    DATA_VECTOR(rts);
    DATA_IVECTOR(xs);
    DATA_VECTOR(pts);
    DATA_VECTOR(wgh);
    int n = rts.size();
    
    // parameters:
    PARAMETER_VECTOR(parm);
    
    // transform parms depending on type_alpha:
    Type v   = parm(0); 
    Type sv  = exp( parm(4) );
    Type sz  = exp( parm(5) );
    Type st0 = exp(parm(6) );
    Type z  = exp( parm(2) ) + 0.5*sz;
    Type t0 = exp( parm(3) ) + 0.5*st0;
    Type a  = exp( parm(1) ) + z + 0.5*sz;
    
    // initialize nll:
    Type nll = Type(0);

    for ( int i = 0; i < n; i++ ) {
        
        // get single trial reaction time and response
        Type rt = rts(i);
        Type x = xs(i);
        
        Type tmp = Type(0);
        if ( x == 0 ) tmp = ddmTMB::ddm7_density_dao2_tmb( rt, a, z, v, t0, sv, sz, st0, pts, wgh, 50 );
        else          tmp = ddmTMB::ddm7_density_dao2_tmb( rt, a, a-z, -v, t0, sv, sz, st0, pts, wgh, 50 );
        
        nll -= log( tmp );
    
    }

    ADREPORT( a );
    ADREPORT( v );
    ADREPORT( z );
    ADREPORT( t0 );
    ADREPORT( sv );
    ADREPORT( sz );
    ADREPORT( st0 );
    
    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this