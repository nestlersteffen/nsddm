
ddm_inits_mu <- function( rts=NULL, type_alpha=NULL, ddm="four" ) 
{
    MUa <- NULL
    if ( ddm == "four" ) {
        v  <- 0
        a  <- runif(1,0.5,2)
        z  <- a/2
        t0 <- min( rts ) + 0.1*runif(1,0,0.1)
        if ( type_alpha == "dao" ) {
            MUa <- c( v, log( a - z ), log( z ), log( t0 )  )
        } else if ( type_alpha == "logit" ) {
            MUa <- c( v, log( a ), qlogis( z/a ), log( t0 )  )
        }
    } else if ( ddm == "seven" ) { 
        #- generate a suitable vector:
        repeat {
            v   <- 0
            sv  <- runif(1, 0, 2)
            a   <- runif(1, 0.5, 2)
            z   <- a/2
            sz  <- runif(1, 0, 0.5)
            st0 <- runif(1, 0, 0.1)
            t0  <- min(rts) + 0.5*st0 + 0.1*runif(1, 0, 0.1)
            #- check constraints:
            if ( (a - z - 0.5*sz) > 0 &
                 (z - 0.5*sz) > 0 &
                (t0 - 0.5*st0) > 0 ) {
                break
            }
        }
        MUa <- c( v, log(a - z - 0.5*sz), log(z - 0.5*sz), log(t0 - 0.5*st0), 
            log( sv ), log(sz), log(st0) )
    }
    return( MUa )
}

ddm_inits_bayes <- function( rts = NULL, I=NULL, K=NULL, args=NULL, ddm=NULL )
{
        
    #- initialize prior information:
    m <- args$bayes_list$m
    M <- args$bayes_list$M
    if ( is.null( args$bayes_list$m) ) m <- rep(0,K)
    if ( is.null( args$bayes_list$M) ) M <- diag(1,K) 
    
    S0  <- args$bayes_list$S0
    nu0 <- args$bayes_list$nu0
    a_d <- 1 / rgamma( K, 0.5, 1)
    if ( args$bayes_list$type_sigma_prior == "huang_wand" ) {
        if (is.null(args$bayes_list$nu0)) nu0 <- 2
        if (is.null(args$bayes_list$S0))  S0  <- 4*diag(1/a_d, K)
    } else {
        if (is.null(args$bayes_list$nu0)) nu0 <- K + 2   
        if (is.null(args$bayes_list$S0))  S0  <- diag(K)
    }
    
    #- get no. of chains:
    nchain <- args$bayes_list$nchain
    
    #- now generate initial values for each chain:
    alpha <- array( 0, dim = c( I, K, nchain ) )
    MUa   <- array( 0, dim = c( K, nchain ) )
    SIGa  <- array( 0, dim = c( K, K, nchain ) ) 
    for ( nc in 1:nchain ) {
        MUa[,nc]    <- ddm_inits_mu( rts=rts, type_alpha=args$type_alpha, ddm=ddm )
        SIGa[,,nc]  <- 4*diag( 1/MCMCpack::rinvgamma(K, 0.5, 1), K )
        alpha[,,nc] <- mvtnorm::rmvnorm( I, MUa[,nc], 0.1*diag(1,K) )
    }

    return( list( alpha=alpha, m=m, M=M, S0=S0, nu0=nu0, a_d=a_d, MUa=MUa, SIGa=SIGa ) )
}