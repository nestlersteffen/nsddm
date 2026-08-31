
#--- function to compute SIGMA from lower cholesky factor:

ddm_compute_sigma <- function( parm_list=NULL, type_sigma="cholesky" ) 
{
    #--- get information:
    L <- parm_list$L 

    #--- compute SIGMA:
    SIGMA <- parm_list$SIGMA
    if ( type_sigma == "cholesky" ) {
        SIGMA <- tcrossprod( L )
    }
    
    #--- add results to parm_list:
    parm_list[["L"]] <- L
    parm_list[["SIGMA"]] <- SIGMA
    return( parm_list )
}