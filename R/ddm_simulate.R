
#------ functions to simulate some data:

ddm4_simulate <- function( ni=100, alpha=NULL, type_alpha=NULL )
{
    #- transform parameters
    v  <- alpha[1]
    t0 <- base::exp( alpha[4] )
    if ( type_alpha == "dao" ) {
        a <- base::exp( alpha[2] ) + base::exp( alpha[3] )
        z <- base::exp( alpha[3] )
    } else if ( type_alpha == "logit" ) { 
        a <- base::exp( alpha[2] );
        w <- 1 / ( 1 + base::exp( -alpha[3] ) )
        z <- a * w
    }
    #- now simulate the data:
    sim_sample <- rtdists::rdiffusion( ni, a=a, v=v, t0=t0, z=z )
    sim_sample$xs <- ifelse( sim_sample$response == "lower", 0, 1 )
    sim_sample <- sim_sample[,c( "rt", "xs" )]
    return( sim_sample )
}

ddm7_simulate <- function( ni=100, alpha=NULL, type_alpha=NULL  )
{
    #- transform parameters
    v   <- alpha[1]
    sv  <- exp(alpha[5])
    sz  <- exp(alpha[6])
    st0 <- exp(alpha[7])
    z   <- exp(alpha[3]) + 0.5*sz
    t0  <- exp(alpha[4]) + 0.5*st0
    a   <- exp(alpha[2]) + z + 0.5*sz
    #- now simulate the data:
    sim_sample <- rtdists::rdiffusion( n = ni, a=a, v=v, z=z, sz=sz, sv=sv, t0=t0-0.5*st0,
        st0=st0, s=1, precision=3 )
    sim_sample$xs <- ifelse( sim_sample$response == "lower", 0, 1 )
    sim_sample <- sim_sample[,c( "rt", "xs" )]
    return( sim_sample )
}

ddm_simulate <- function( I=100, ni=100, MU=NULL, SIGMA=NULL, 
    type_alpha="dao", ddm="four" )
{
    #- draw random effects:
    alpha <- mvtnorm::rmvnorm( I, MU, SIGMA )
    
    #- which function to use:
    ddm_fun <- switch( ddm,
        "four"  = ddm4_simulate,
        "seven" = ddm7_simulate,
        stop( paste( "Only the 4- or 7-parameter DDM can be simulated." ) )
    )

    #- now simulate the data:
    sim_list <- vector("list", I)
    for ( i in seq(I) ) {
        #- get the random effect:
        alphai <- alpha[i,]
        #- alphai is on an unconstrained scale, we have to transform it:
        tmp_sample <- ddm_fun( ni=ni, alpha=alphai, type_alpha=type_alpha )
        tmp_sample$id <- i
        sim_list[[i]] <- tmp_sample
    }
    sim_sample <- do.call(rbind, sim_list)
    return( sim_sample )
}

rddm7 <- function( n, a, v, z, t0, sv=0, sz=0, st0=0, s=1, dt=0.001, 
                   max_time=10 ) 
{
    #- storage:
    rt       <- numeric(n)
    response <- integer(n)
    
    #- sqrt(dt) for the Wiener increments:
    sqrt_dt <- sqrt(dt)
    
    for ( i in seq_len(n) ) {
        
        #- draw trial-specific parameters:
        v_i  <- rnorm(1, mean=v, sd=sv)
        z_i  <- runif(1, min=z - sz/2, max=z + sz/2)
        t0_i <- runif(1, min=t0 - st0/2, max=t0 + st0/2)
        
        #- simulate the diffusion process:
        x <- z_i
        t <- 0
        while ( x > 0 && x < a && t < max_time ) {
            x <- x + v_i*dt + s*sqrt_dt*rnorm(1)
            t <- t + dt
        }
        
        #- record outcome:
        if ( t >= max_time ) {
            #- no boundary hit: mark as NA (rare in sensible parameter ranges)
            rt[i]       <- NA
            response[i] <- NA
        } else {
            rt[i]       <- t + t0_i
            response[i] <- ifelse( x >= a, 1, 0 )   # 1 = upper, 0 = lower
        }
    }
    
    return( data.frame( rt=rt, xs=response ) )
}

rddm7_v2 <- function( n, a, v, z, t0, sv=0, sz=0, st0=0, s=1, dt=0.001, 
                   max_time=10 ) 
{
    rt       <- numeric(n)
    response <- integer(n)
    
    sqrt_dt <- sqrt(dt)
    s2_dt   <- s*s*dt
    
    for ( i in seq_len(n) ) {
        
        #- draw trial-specific parameters:
        v_i  <- rnorm(1, mean=v, sd=sv)
        z_i  <- runif(1, min=z - sz/2, max=z + sz/2)
        t0_i <- runif(1, min=t0 - st0/2, max=t0 + st0/2)
        
        #- simulate the diffusion process:
        x_prev <- z_i
        t      <- 0
        hit    <- NA_integer_
        
        while ( t < max_time ) {
            x_new <- x_prev + v_i*dt + s*sqrt_dt*rnorm(1)
            t     <- t + dt
            
            #- direct boundary hit?
            if ( x_new >= a ) {
                hit <- 1L
                break
            }
            if ( x_new <= 0 ) {
                hit <- 0L
                break
            }
            
            #- brownian bridge correction: did we cross in between?
            #  (only when both endpoints are strictly inside)
            p_upper <- exp( -2*(a - x_prev)*(a - x_new)/s2_dt )
            p_lower <- exp( -2*x_prev*x_new/s2_dt )
            
            u <- runif(1)
            if ( u < p_upper ) {
                hit <- 1L
                break
            }
            if ( u < p_upper + p_lower ) {
                hit <- 0L
                break
            }
            
            x_prev <- x_new
        }
        
        if ( is.na(hit) ) {
            rt[i]       <- NA
            response[i] <- NA
        } else {
            rt[i]       <- t + t0_i
            response[i] <- hit
        }
    }
    
    return( data.frame( rt=rt, xs=response ) )
}