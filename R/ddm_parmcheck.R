
#---- function to check parms

ddm4_parmcheck <- function( a = NULL, t0 = NULL, z = NULL, w = NULL )  
{
	if ( z < 0 || z > a || a <= 0 || t0 <= 0 ) {
	    return( FALSE )
	}
	if ( w <= 0 || w >= 1 ) {
	    return( FALSE )
	}
	if ( t0 <= 0 ) {
	    return( FALSE )
	}
	return( TRUE )
}

ddm7_parmcheck <- function( a = NULL, t0 = NULL, z = NULL, sv = NULL, sz=NULL,st0 = NULL )  
{
	if ( a <= 0 || sz > a || z + 0.5*sz > a ) {
        return( FALSE )
    }
    if ( z < 0 || sz < 0 || z - 0.5*sz < 0 || sz > 2 ) {
        return( FALSE )
    }
    if ( t0 < 0 || st0 < 0 || st0 > 1 ) {
        return( FALSE )
    }
    if ( sv < 0 || sv > 3 ) {
        return( FALSE )
    }
	return( TRUE )
}