
ddm_agh_llfct <- function( parm=NULL, data_list=NULL, parm_list=NULL, 
    parm_table=NULL, pts_array=NULL, wgh=NULL, args=NULL, weights=NULL, 
    ddm="four", add_terms=TRUE )
{
  
    #- insert parm:
    parm_list <- ddm_include_free_parameters( parm=parm, parm_list=parm_list, 
        parm_table=parm_table )
  
    #- get parameter matrices:
    parm_list <- ddm_compute_sigma( parm_list=parm_list, type_sigma=args$type_sigma )
    MU        <- parm_list[["MU"]]
    L         <- parm_list[["L"]]
    SIGMA     <- parm_list[["SIGMA"]]

    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    K    <- data_list[["K"]]
    info <- data_list[["info"]]
  
    #- select function depending on ddm:
    ddm_fun <- switch( ddm,
        "four"  = ddm4_random_nllfct,
        "seven" = ddm7_random_nllfct,
        stop( paste( "The 4- or 7-parameter DDM can be estimated." ) )
    )

    #- make vector to save ll-values
    res <- rep( 0, I )
  
    #- compute loglik for each person
    for ( i in seq( I ) ) {
    
        #- collect person data:
        idx1 <- info[i,3] - info[i,2] + 1 
        idx2 <- info[i,3]
        rti  <- rts[idx1:idx2]
        xsi  <- xs[idx1:idx2]
        
        #- collect mode info:
        ranef <- modes_list[[ i ]]
        
        #- for each random point we compute the Likelihood
        npts <- nrow( pts_array[,,i] )
        fapp <- rep( 0, npts )
          
        for ( ptx in seq( npts ) ) {
            
            fapp[ ptx ] <- -1*ddm_fun( us = pts_array[ptx,,i], 
                rt=rti, xs=xsi, MU=MU, SIGMA=SIGMA, args=args, 
                weights=weights )
        
        }

        #- final computations:
        res[ i ] <- ddm_sumexplog( fapp + log( wgh[,i] ) ) 

    }

    #- add approximation terms:
    if ( add_terms ) {
        if ( args$maxlik_list$method == "AGH" ) {
            sqrtDET <- do.call("c", lapply( modes_list, function(x) as.numeric( sqrt( x$DET ) ) ) )
            res <- res + log( ( 2*pi )^(K/2) * sqrtDET )
        } else {
            res <- res - log( npts )
        }
    }
  
    #- output:
    out <- -1*sum( res )
    return( out )
}