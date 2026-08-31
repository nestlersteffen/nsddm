/// @file ddm7_getmodes_tmb.hpp

// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type ddm7_getmodes_tmb( objective_function<Type>* obj )
{
    
    using namespace density;

    // data:
    DATA_VECTOR(rts);
    DATA_IVECTOR(xs);
    DATA_VECTOR(pts);
    DATA_VECTOR(wgh);
    DATA_VECTOR(MU);
    DATA_MATRIX(SIGMA);
    int n = rts.size();
    
    // parameters:
    PARAMETER_VECTOR( us ); 
    
    // compute random effects density:
    vector<Type> eus = us - MU;
    Type nll = MVNORM( SIGMA )( eus ); 
    
    // transform parms depending on type_alpha:
    Type v   = us(0); 
    Type sv  = exp( us(4) );
    Type sz  = exp( us(5) );
    Type st0 = exp( us(6) );
    Type z   = exp( us(2) ) + 0.5*sz;
    Type t0  = exp( us(3) ) + 0.5*st0;
    Type a   = exp( us(1) ) + z + 0.5*sz;
    

    for ( int i = 0; i < n; i++ ) {
        
        // get single trial reaction time and response
        Type rt = rts(i);
        Type x = xs(i);
        
        Type tmp = Type(1e-27);
        if ( x == 0 ) tmp = ddmTMB::ddm7_density_dao2_tmb( rt, a, z, v, t0, sv, sz, st0, pts, wgh, 50 );
        else          tmp = ddmTMB::ddm7_density_dao2_tmb( rt, a, a-z, -v, t0, sv, sz, st0, pts, wgh, 50 );
        
        nll -= log( tmp );
    
    }
    
    return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this