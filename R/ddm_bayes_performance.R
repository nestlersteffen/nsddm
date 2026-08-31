
ddm_iact <- function( chain, L_max = 1000 ) 
{
    #- chain length:
    M <- length( chain )

    #- compute autocorrelations:
    ac <- acf( chain, lag.max = min( L_max, M-1 ), plot = FALSE )$acf[,,1]
    
    #- determine cutoff L_M:
    se_ac <- 2 / sqrt( M )
    L_M <- min( which( abs( ac ) < se_ac ) - 1, L_max )
    if ( is.na( L_M ) || L_M < 1 ) L_M <- 1
    
    #- compute IACT:
    iact <- 1 + 2 * sum( ac[2:(L_M+1)] )
    
    return( iact )
}

ddm_rhat <- function( chains )
{
    #- chains ist ein Array mit Dimensionen c(biter, K, nchain)
    #- oder eine Liste von Matrizen mit Dimensionen c(biter, K)
    
    #- get dimensions:
    biter  <- dim(chains)[1]
    K      <- dim(chains)[2]
    nchain <- dim(chains)[3]
    
    #- compute Rhat for each parameter:
    rhat <- rep( 0, K )
    for ( k in 1:K ) {
        
        #- get chains for parameter k:
        chains_k <- chains[,k,]  # biter x nchain Matrix
        
        #- within-chain variance (W):
        s2_j <- apply( chains_k, 2, var )  # Varianz pro Kette
        W    <- mean( s2_j )
        
        #- between-chain variance (B):
        mu_j  <- apply( chains_k, 2, mean )  # Mittelwert pro Kette
        mu    <- mean( mu_j )
        B     <- biter / (nchain - 1) * sum( (mu_j - mu)^2 )
        
        #- pooled variance estimate:
        V_hat <- (1 - 1/biter) * W + (1/biter) * B
        
        #- Rhat:
        rhat[k] <- sqrt( V_hat / W )
    }
    
    return( rhat )
}