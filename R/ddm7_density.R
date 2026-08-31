
#--- ------------------------------------------------------------------
#---  standard seven-parameter DDM density following Dao et al., 2025 
#--- -----------------------------------------------------------------

ddm7_large_integrand <- function( tx=NULL, a=NULL, z=NULL, v=NULL, 
    t0=NULL, sv=NULL, kmax=NULL, delta=NULL ) 
{
  
    #- compute "true" reaction time:
    td <- tx - t0

    #- compute multiplicative terms:
    mult1 <- pi/( a*a*sqrt(1 + td*sv*sv) )
    mult2 <- base::exp(-0.5*( ( v*v*td + 2*v*z - z*z*sv*sv )/(td*sv*sv + 1) ) )
    
    #- approximate the infinite sum:
    sum_hist <- rep(0,3)
    for ( k in seq( kmax ) ) {
        
        # for things to save
        sum_hist[1] <- sum_hist[2]
        sum_hist[2] <- sum_hist[3]
    
        # compute elements of the sum:
        pt_1 <- base::sin( (pi*k*z)/a )
        pt_2 <- base::exp( (-0.5*pi*pi*k*k*td)/(a*a) )
        
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

    #- final result:
    return( sum_hist[3]*mult1*mult2 )
}

ddm7_small_integrand <- function( tx=NULL, a=NULL, z=NULL, v=NULL,
    t0=NULL, sv=NULL, kmax=NULL, delta=NULL)
{
    #- compute "true" reaction time:
    td <- tx - t0

    #-multiplicative terms ( see Dao et al., page 7 )
    mult1 <- (a*td^(-3/2))/sqrt(2*pi*(1 + td*sv*sv))
    log_mult2 <- ( -0.5*(v*v)/(sv*sv) ) + ( (v-z*sv*sv)^2 / (2*sv*sv*(1 + td*sv*sv) ) )
    mult2 <- exp( log_mult2 )

    #- compute the infinite sum in pairs of (-k, k)
    term_in_the_sum <- function( k ) {
        coeff <- z / a + 2 * k
        coeff * exp(-(z + 2 * k * a)^2 / (2 * td) )
    }
    
    sum_hist    <- rep(0, 3)
    sum_hist[3] <- term_in_the_sum(k=0)  # k = 0
    for ( m in seq( kmax ) ) {
        
        # for things to save:
        sum_hist[1] <- sum_hist[2]
        sum_hist[2] <- sum_hist[3]
        
        # compute the sum and save:
        sum_hist[3] <- sum_hist[2] + term_in_the_sum(k=m) + term_in_the_sum(k=-m)
    
        # check:
        if ( abs(sum_hist[1] - sum_hist[2]) < delta &&
             abs(sum_hist[2] - sum_hist[3]) < delta) break
    }

    result <- mult1 * mult2 * sum_hist[3]
    if ( !is.finite (result ) || result <= 0 ) return(0)
    return(result)
}

ddm7_dao <- function( tx=NULL, x=NULL,a=NULL, z=NULL, v=NULL, t0=NULL, sv=NULL, 
    sz=NULL, st0=NULL, kmax=NULL, delta=NULL, err=0.000001 ) 
{
  
    #- get the Gauss-Legendre point:
    pts <- c(-0.973906528517172, -0.865063366688985, -0.679409568299024,
             -0.433395394129247, -0.148874338981631,  0.148874338981631,
              0.433395394129247,  0.679409568299024,  0.865063366688985, 0.973906528517172)
    wgh <- c(0.066671344308688, 0.149451349150581, 0.219086362515982,
             0.269266719309996, 0.295524224714753, 0.295524224714753,
             0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688)
  
    #- change variable:
    if ( x == 1 ) {
        z <- ( a - z )
        v <- -v
    } 

    # transform points for t0 and z
    z_pts  <- z  + (sz/2)  * pts   # 10 Punkte im z-Intervall
    t0_pts <- t0 + (st0/2) * pts   # 10 Punkte im tau-Intervall
  
    total <- 0
    for (i in seq_along(pts)) {
        
        for (j in seq_along(pts)) {
        
            #- get current z und tau:
            z_ij  <- z_pts[i]
            t0_ij <- t0_pts[j]  
        
            #- check conditions:
            if (z_ij <= 0 || z_ij >= a) next
            if (t0_ij <= 0) next
            if (tx - t0_ij <= 0) next  

            #- determine optimal k:
            td_ij <- tx - t0_ij
            tt_ij <- td_ij / (a^2)
            ks <- get_ks(t=tt_ij, err=err)
            kl <- get_kl(t=tt_ij, err=err)
        
            if (ks < kl) {
                f_ij <- ddm7_small_integrand(tx=tx, a=a, z=z_ij, v=v, 
                    t0=t0_ij, sv=sv, kmax=ks, delta=delta)
            } else {
                f_ij <- ddm7_large_integrand(tx=tx, a=a, z=z_ij, v=v, 
                    t0=t0_ij, sv=sv, kmax=kl, delta=delta)
            }   
        
            #- compute integrand:
            total <- total + wgh[i] * wgh[j] * f_ij
        
        }
    }
  
    # Faktor 1/4 (die s_z und s_tau kürzen sich heraus)
    return(total / 4)
}

#--- ------------------------------------------------------------------
#---  standard seven-parameter DDM density following Tuerlinckx (2004) 
#--- -----------------------------------------------------------------

ddm7_dstand <- function( tx=NULL, x=NULL, a=NULL, t0=NULL, z=NULL, v=NULL, 
    st0=NULL, sz=NULL, sv=NULL, s2=1, kmax=NULL, delta=NULL, epsilon=1e-7, 
    min_RT=0.001)#, err=NULL ) 
{
  
    #- points and weights for the distribution function:
    pts_us <- c(-3.43615911883774,-2.53273167423279,-1.75668364929988,-1.03661082978951,
        -0.34290132722370, 0.34290132722370, 1.03661082978951, 1.75668364929988, 
        2.53273167423279, 3.43615911883774 )
    wgh_us <- c(0.00000764043286, 0.00134364574678, 0.03387439445548, 0.24013861108230,
        0.61086263373530, 0.61086263373530, 0.24013861108230, 0.03387439445548, 
        0.00134364574678, 0.00000764043286 )

    pts_zs <- c(-0.973906528517172,-0.865063366688985,-0.679409568299024,-0.433395394129247,
      -0.148874338981631, 0.148874338981631, 0.433395394129247, 0.679409568299024,
      0.865063366688985, 0.973906528517172 )
    wgh_zs <- c( 0.066671344308688,0.149451349150581, 0.219086362515982,0.269266719309996, 
        0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982,
        0.149451349150581, 0.066671344308688 )

    #- collect parameters:
    a2 <- a*a
    Zu <- z + sz/2
    Zl <- z - sz/2
    Tl <- t0 - st0/2 
    Tu <- min(tx, t0 + st0/2 )
      
    #- prepare points and transform them
    us    <- sqrt(2)*pts_us*sv + v
    w_us  <- (1/sqrt(pi))*wgh_us
    nr_us <- 10
      
    zs    <- pts_zs*(sz/2) + z 
    w_zs  <- wgh_zs*(sz/2) 
    nr_zs <- 10 

    #- change variable:
    if ( x == 1 ) {
        Zu   <- ( a - z ) + sz/2
        Zl   <- ( a - z ) - sz/2 
        us   <- -1*us
    } 

    #- compute constant things:
    tTu  <- -0.5*(tx - Tu)
    tTl  <- -0.5*(tx - Tl)

    if ( ( tx - t0 + st0/2 ) <= min_RT ) {
        return( 0 )
    } 

    # standard case:
    if ( tx > t0 + st0/2 ) { 
        
        sum_hist <- rep(0,3)
        for (k in 1:kmax) {
                
            # for things to save
            sum_hist[1:2] <- sum_hist[2:3]
            sum_us <- 0 # sum across quadrature points
                
            # some further constant things:
            pika  <- (pi*k)/a 
            pikas <- (pi*pi*k*k*s2)/a2
            
            # gauss quadrature:
            for ( uu in 1:nr_us ) {
                #- some tmps:
                u2s2pikas <- us[uu]^2/s2 + pikas
                #- start computations:
                t1     <- ( u2s2pikas )^(-2) 
                t2.1   <- exp( -us[uu]*Zu/s2 )*( -us[uu]*sin(pika*Zu)/s2 - pika*cos(pika*Zu) )
                t2.2   <- exp( -us[uu]*Zl/s2 )*( -us[uu]*sin(pika*Zl)/s2 - pika*cos(pika*Zl) )
                t3     <- exp( u2s2pikas*tTu ) - exp( u2s2pikas*tTl )
                sum_us <- sum_us + t1*( t2.1 - t2.2 )*t3*w_us[uu]
            } # end uu - loop
                
            sum_hist[3] <- sum_hist[2] + k*sum_us
            if ( ( abs( sum_hist[1] - sum_hist[2] ) < delta ) & 
                 ( abs( sum_hist[2] - sum_hist[3] ) < delta ) & 
                 ( sum_hist[3] > 0 ) ) {
                break
            }
        } # end k - loop
        out <- 2*sum_hist[3]*s2*s2*pi/(a2*sz*st0)
    
    } else if ( tx <= t0 + st0/2 ) {
      
        sum_us <- 0
        for ( uu in 1:nr_us ) {
            
            if ( abs( us[uu] ) > epsilon ) {
                
                sum_zs <- 0
                for ( zz in 1:nr_zs ) {
                    
                    # approximate first sum (see Equation 22):
                    zzz1 <- ( a - zs[zz] )*x + zs[zz]*(1 - x)
                    zzz2 <- ( a - zs[zz] )*(1-x) + zs[zz]*x  
                    sum1 <- (a2/(pi*s2))*sinh( zzz2*us[uu]/s2 )/sinh( (us[uu]*a )/s2 )
                    
                    # get the second sum:
                    sum_hist <- rep(0,3)
                    for (k in 1:kmax) {
                        sum_hist[1:2] <- sum_hist[2:3]
                        pika  <- (pi*k)/a 
                        u2s2pikas <- us[uu]^2/s2 + (pi*pi*k*k*s2)/a2
                        #- tmp computations:
                        t1 <- ( u2s2pikas )^(-1) 
                        t2 <- sin( pika*zzz1 )
                        t3 <- exp( u2s2pikas*tTl )
                        sum_hist[3] <- sum_hist[2] + k*t1*t2*t3
                        if ( ( abs( sum_hist[1] - sum_hist[2] ) < delta ) & 
                             ( abs( sum_hist[2] - sum_hist[3] ) < delta ) & 
                             ( sum_hist[3] > 0 ) ) {
                            break
                        }
                    } # for loop k 
                    sum_zs <- sum_zs + w_zs[zz]/sz*( sum1 - 2*sum_hist[3] )*(pi*s2)/(a2*st0)*exp((-zzz1*us[uu])/s2);
                
                } # end zs - loop
            
            } else if ( abs( us[uu] ) < epsilon ) {
                
                sum_hist <- rep(0,3)
                sl  <- -(Zl*pi^2)/(2*a) + (Zl^2*pi^2)/(4*a2)
                su  <- -(Zu*pi^2)/(2*a) + (Zu^2*pi^2)/(4*a2)
                for (k in 1:kmax) {
                    sum_hist[1:2] <- sum_hist[2:3]
                    pika  <- (pi*k)/a 
                    pika2 <- (pi*pi*k*k*s2)/a2;
                    #- tmp computations:
                    t1 <- cos( pika*Zl )
                    t2 <- cos( pika*Zu )
                    t3 <- exp( pika2*tTl )
                    sum_hist[3] <- sum_hist[2] + (k^(-2))*(t1-t2)*t3
                    if ( ( abs( sum_hist[1] - sum_hist[2] ) < delta ) & 
                         ( abs( sum_hist[2] - sum_hist[3] ) < delta ) & 
                         ( sum_hist[3] > 0 ) ) {
                        break
                    }
                } # end k - loop
            
                sum_zs <- ((2*a)/(st0*sz*s2*pi^2))*( sl - su - sum_hist[3] )
            
            } # end zs - loop
            
            sum_us <- sum_us + sum_zs*w_us[uu]
          
        } # end us - loop
        
        out <- sum_us
    }

    #- a final check:
    if ( out <= 0 ) {
        return( 0 )
    }
    
    # output:
    return( out )
}

#--- --------------------------------------------------------------
#---  wrapper function for a vector of reaction times and choices
#--- --------------------------------------------------------------

ddm7_pdf <- function( rt=NULL, x=NULL, a=NULL, t0=NULL, z=NULL, v=NULL, sv=NULL, 
    sz=NULL, st0=NULL, type_ddm="std", kmax=5000, delta=1e-29 )
{

    #- how many densities to compute?
    nn <- length( rt )
    densities <- rep( 0, length = nn )

    #- check parms?
    check <- ddm7_parmcheck( a=a, t0=t0, z=z, sv=sv, sz=sz, st0=st0 )  
    if ( !check ) return( densities )

    #- select the function:
    ddm7_fun <- switch( type_ddm,
        "std" = ddm7_dstand,
        "dao" = ddm7_dao,
        stop( paste( "Unknown type of ddm density calculation method:", type_ddm ) )
    )
    
    #- let's go:
    for ( i in seq( nn ) ) {
        # if ( ( rt[i] - t0 ) <= 0 ) next
        densities[i] <- ddm7_fun( tx=rt[i], x=x[i], a=a, z=z, v=v, t0=t0, sv=sv, 
            sz=sz, st0=st0, kmax=kmax, delta=delta )
    }

    return( densities )
} 

#--- --------------------------------------------------------------
#---      log-likelihood functions given a parameter vector u
#--- --------------------------------------------------------------

#- standard data-llfct::

ddm7_logLdata <- function( us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL ) 
{ 
    #- Step 1: transform parameters
    v   <- us[1]
    sv  <- exp(us[5])
    sz  <- exp(us[6])
    st0 <- exp(us[7])
    z   <- exp(us[3]) + 0.5*sz
    t0  <- exp(us[4]) + 0.5*st0
    a   <- exp(us[2]) + z + 0.5*sz
    
    #- Step 2: compute density values
    if ( args$use_lan ) {
        nll <- sum( ddm7_lanll_weights( weights, a, v, t0, z, sv, sz, st0, rt, xs, K = length(us) ) )
    } else {
        ll  <- ddm7_pdf(rt=rt, x=xs, a=a, t0=t0, z=z, v=v, sv=sv, sz=sz, st0=st0, 
            type_ddm=args$type_ddm, kmax=args$kmax, delta=args$delta )
        ll[ ll == 0 | ll == Inf ] <- 1e-29
        nll <- sum(log(ll))
    }
    return( nll )
}

#- standard negative data-llfct::

ddm7_nllfct <- function(us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL ) {
  -ddm7_logLdata(us=us, rt=rt, x=xs, args=args, weights=weights )
}

#-  standard negative data-llfct including random effects

ddm7_random_nllfct <- function( us=NULL, rt=NULL, xs=NULL, 
    MU=NULL, SIGMA=NULL, args=NULL, dnn=NULL, weights=NULL ) 
{
    #- data likelihood:
    log_pData <- ddm7_nllfct( us=us, rt=rt, xs=xs, args=args, weights=weights ) 
    #- log prior:
    log_prior <- mvtnorm::dmvnorm( us, MU, SIGMA, log = TRUE )
    return( log_pData - log_prior )
}