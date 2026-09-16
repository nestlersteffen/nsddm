
#--- --------------------------------------------------------------
#---  standard four-parameter DDM density (see Tuerlinckx, 2004 )
#--- --------------------------------------------------------------

ddm4_dstand <- function( t = NULL, a = NULL, w = NULL, v = NULL, 
    s2=NULL, kmax=NULL, delta=NULL, err=NULL ) 
{
  
    #- compute multiplicative term:
    mult <- (s2*pi)/(a*a)*base::exp( ( -w*a*v - 0.5*v*v*t )/s2 )
    
    #- approximate the infinite sum:
    sum_hist <- rep(0,3)
    for ( k in seq( kmax ) ) {
        
        # for things to save
        sum_hist[1] <- sum_hist[2]
        sum_hist[2] <- sum_hist[3]
    
        # compute elements of the sum:
        pt_1 <- base::sin( pi*k*w )
        pt_2 <- base::exp( (-0.5*pi*pi*k*k*s2*t)/(a*a) )
        
        # compute the sum and save:
        sum_hist[3] <- sum_hist[2] + k*pt_1*pt_2
        
        # check:
        if ( ( base::abs( sum_hist[1] - sum_hist[2] ) < delta ) & 
             ( base::abs( sum_hist[2] - sum_hist[3] ) < delta ) & 
             ( sum_hist[3] > 0 ) ) {
                break;
        }
    
    } # end k - loop

    #- a final check:
    if ( sum_hist[3] <= 0 ) {
        return( 0 )
    }

    #- output:
    return( sum_hist[3]*mult );
}

#--- ----------------------------------------------------------------------
#--- Navarro & Fuss four-parameter DDM density (see Navarro & Fuss, 2009 )
#--- ----------------------------------------------------------------------

MAX <- .Machine$integer.max

#--- get the sum limits for small and large time

get_ks <- function( t = NULL, err = NULL ) 
{
    if (err * 2 * base::sqrt(2 * pi * t) < 1) {
        ks <- 2 + base::sqrt(-2 * t * base::log(2 * err * base::sqrt(2 * pi * t)))
        bc <- base::sqrt(t) + 1
        if (ks > MAX || bc > MAX) return(MAX)
        return( base::ceiling( base::max( ks, bc ) ) )
    }
    return( 2 )
}

get_kl <- function( t = NULL, err = NULL ) {
    bc <- 1 / (pi * base::sqrt(t))
    if (bc > MAX) return(MAX)
    if (err * pi * t < 1) {
        kl <- base::sqrt(-2 * base::log(pi * t * err) / (pi^2 * t))
        if (kl > MAX) return(MAX)
        return( base::ceiling( base::max(kl, bc) ) )
    }
    return( base::ceiling( bc ) )
}

#--- now the density:

ddm4_dnavfuss <- function(t = NULL, a = NULL, w = NULL, v = NULL, s2=NULL, kmax=NULL, delta=NULL, err=NULL ) 
{
    
    #- transform time:
    tt <- t / (a^2)
      
    #- get ks and kl:
    ks <- get_ks( t=tt, err=err)
    kl <- get_kl( t=tt, err=err)
      
    #- compute multiplicative term
    mult <- base::exp(-w * a * v - 0.5 * v^2 * t) / (a^2)
      
    #- approximate sum
    if (ks < kl) {  # Small-time expansion
        gamma   <- -1 / (2 * tt)
        k_lower <- -base::floor((ks - 1) / 2)
        k_upper <-  base::floor(ks / 2)
        ks_seq  <- k_lower:k_upper
        tmp_sum <- base::sum( (w + 2 * ks_seq) * base::exp(gamma * (w + 2 * ks_seq)^2) )
        dens    <- tmp_sum / base::sqrt(2 * pi * tt^3)
    } else {        # Large-time expansion
        gamma   <- -0.5 * pi^2 * tt
        kl_seq  <- 1:kl
        tmp_sum <- base::sum( base::sin(pi * kl_seq * w) * base::exp(gamma * kl_seq^2) * kl_seq)
        dens    <- tmp_sum * pi
    }

    return(dens * mult)

}

#--- --------------------------------------------------------------
#---  wrapper function for a vector of reaction times and choices
#--- --------------------------------------------------------------

ddm4_pdf <- function( rt = NULL, x = NULL, 
    a = NULL, t0 = NULL, w = NULL, z = NULL, v = NULL, 
    type_ddm="std", s2=1, kmax=5000, delta=1e-29, err=0.000001 )
{

    #- how many densities to compute?
    nn <- length( rt )
    densities <- rep( 0, length = nn )
    
    #- parametrization?
    if ( is.null( w ) & !is.null( z ) ) { # z parametrization
        w <- z/a
    } else {
        z <- a*w
    }

    #- check parms?
    check <- ddm4_parmcheck( a=a, t0=t0, z=z, w=w )  
    if ( !check ) return( densities )

    #- select the function:
    ddm4_fun <- switch( type_ddm,
        "std"     = ddm4_dstand,
        "navfuss" = ddm4_dnavfuss,
        stop( paste( "Unknown type of ddm density calculation method:", type_ddm ) )
    )

    #- compute ts:
    ts <- rt - t0 

    #- let's go:
    for ( i in seq( nn ) ) {
        if ( ts[i] <= 0 ) next
        if ( x[i] == 0 ) {
            densities[i] <- ddm4_fun( t=ts[i], a=a, w=w, v=v, s2=s2, kmax=kmax, delta=delta, err=err )
        } else { # response is "upper" or 1
            densities[i] <- ddm4_fun( t=ts[i], a=a, w=1-w, v=-v, s2=s2, kmax=kmax, delta=delta, err=err )
        }
    }

    return( densities )
} 

#--- --------------------------------------------------------------
#---      log-likelihood functions given a parameter vector u
#--- --------------------------------------------------------------

#- standard data-llfct:

ddm4_logLdata <- function( us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL, dnn=NULL ) 
{ 
    #- Step 1: transform parameters
    if ( args$type_alpha == "dao") {
        v  <- us[1]
        z  <- exp(us[3])
        t0 <- exp(us[4])
        a  <- exp(us[2]) + z
    } else if (args$type_alpha == "logit") {
        v  <- us[1]
        w  <- plogis(us[3])
        t0 <- exp(us[4])
        a  <- exp(us[2])
        z  <- w * a
    }
      
    #- step 2: compute likelihood
    if ( args$use_lan ) {
        if ( args$use_tf ) {
            print("lan with tf")
            nll <- ddm4_lan_llfct_dnn( dnn=dnn, rt=rt, xs=xs, v=v, a=a, z=z, t0=t0 )    
        } else {
            print("lan with weights")
            nll <- ddm4_lanll_weights( weights=weights, a=a, v=v, t0=t0, z=z, rt=rt, xs=xs )
        }
    } else {
        ll <- ddm4_pdf(rt=rt, x=xs, a=a, t0=t0, z=z, v=v, type_ddm=args$type_ddm, kmax=args$kmax, delta=args$delta )
        ll[ll == 0 | ll == Inf] <- 1e-29
        nll <- sum(log(ll))
    }
    
    return( nll )
}

#- standard negative data-llfct:

ddm4_nllfct <- function(us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL, dnn=NULL) {
  -ddm4_logLdata(us=us, rt=rt, x=xs, args=args, weights=weights, dnn=dnn)
}

#- standard negative data-llfct including random effects

ddm4_random_nllfct <- function( us=NULL, rt=NULL, xs=NULL, 
    MU=NULL, SIGMA=NULL, args=NULL, weights=NULL, dnn=NULL ) 
{
    #- data likelihood:
    log_pData <- ddm4_nllfct( us=us, rt=rt, xs=xs, args=args, weights=weights, dnn=dnn ) 
    #- log prior:
    log_prior <- mvtnorm::dmvnorm( us, MU, SIGMA, log = TRUE )
    return( log_pData - log_prior )
}