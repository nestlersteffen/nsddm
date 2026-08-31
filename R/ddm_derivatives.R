
ddm_derivatives <- function( MU=NULL, L=NULL, 
    invSIGMA=NULL, pos=NULL, type_jj=NULL, 
	alphainvSIGMA=NULL ) 
{
  
    #- set derivative to zero:
    res <- 0
      
    #- formulas for the derivatives:
    if ( type_jj == "MU" ) { 
        #- 1_ij - matrix:
        p <- nrow( MU )
        ONE_ij <- matrix( 0, nrow = p, ncol = 1 ) 
        ONE_ij[ pos[1], 1 ] <- 1
        #- derivative: 
        res <- alphainvSIGMA%*%ONE_ij
    }

    if ( type_jj == "L" ) {
	    #- 1_ij - matrix
	    p <- nrow( invSIGMA )
	    ONE_ij <- diag( 0, p ) 
	    ONE_ij[ pos[1], pos[2] ] <- 1  
	    #- derivative: 
	    dSIGMA <- ONE_ij %*% t( L ) + L %*% t( ONE_ij )
	    #- true ll derivative:
	    pt1 <- -0.5*sum( invSIGMA * dSIGMA )
	    pt2 <-  0.5*( alphainvSIGMA%*%( dSIGMA%*%t(alphainvSIGMA) ) )
	    res <-  pt1 + pt2   
  	}

    if ( type_jj == "SIGMA" ) {
	    #- 1_ij - matrix
	    p <- nrow( invSIGMA )
	    ONE_ij <- diag( 0, p ) 
	    ONE_ij[ pos[1], pos[2] ] <- 1  
	    #- derivative: 
	    dSIGMA <- ONE_ij + t( ONE_ij ) - ONE_ij %*% ONE_ij
	    #- true ll derivative:
	    pt1 <- -0.5*sum( invSIGMA * dSIGMA )
	    pt2 <-  0.5*( alphainvSIGMA%*%( dSIGMA%*%t(alphainvSIGMA) ) )
	    res <-  pt1 + pt2   
  	}
  
    # if ( mat == "S" ) {
    #     #- 1_ij - matrix
    #     p <- nrow( iSIGMA )
    #     ONE_ij <- diag( 0, p ) 
    #     ONE_ij[ pos[1], pos[2] ] <- 1 
    #     #- transform R in case type_sigma == "TMB":
    #     tR <- R
    #     if ( control$type_sigma == "TMB" ) {
    #         tR <- mptmem_compute_r_tmb( R  )
    #     } 
    #     #- derivative: 
    #     dSIGMA <- t( ONE_ij ) %*% tR %*% S + t( S ) %*% tR %*% ONE_ij
    #     #- true ll derivative:
    #     pt1 <- -0.5*sum( iSIGMA * dSIGMA )
    #     pt2 <-  0.5*( usiSIGMA%*%( dSIGMA%*%t(usiSIGMA) ) )
    #     res <-  pt1 + pt2   
    # }

    # if ( mat == "R" ) {
    #     #- 1_ij - matrix
    #     p  <- nrow( iSIGMA )
    #     dR <- diag( 0, p ) 
    #     if ( control$type_sigma == "COR" ) {
    #         dR[ pos[1], pos[2] ] <- 1  
    #         dR[ pos[2], pos[1] ] <- 1 
            
    #     } else {
    #         dR <- mptmem_compute_dr_tmb( R, pos  )
    #     }    
    #     #- derivative: 
    #     dSIGMA <- t( S ) %*% dR %*% S
    #     #- true ll derivative:
    #     pt1 <- -0.5*sum( iSIGMA * dSIGMA )
    #     pt2 <-  0.5*( usiSIGMA%*%( dSIGMA%*%t(usiSIGMA) ) )
    #     res <-  pt1 + pt2   
    # }

    #- output
    return( res )

} 