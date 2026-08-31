
ddm_change_to_sigma <- function( parm_list=NULL, parm_table=NULL, args=NULL )
{
    
    #- insert current parm:
    parm_list <- ddm_include_free_parameters( parm=parm_table$est, parm_list=parm_list, 
        parm_table=parm_table )
  
    #- get parameter matrices:
    parm_list <- ddm_compute_sigma( parm_list=parm_list, type_sigma=args$type_sigma )
    SIGMA 	  <- parm_list$SIGMA
    SIGMA_elements <- SIGMA[ lower.tri( SIGMA, diag = TRUE ) ]
    
	#- change L part:
    idx <- which( parm_table$type == "L" )
    parm_table[idx,"type"]  <- "SIGMA"
	parm_table[idx,"est"]   <- SIGMA_elements
	parm_table[idx,"ntype"] <- 2

	#- return result
    return( parm_table = parm_table )

}