
#--- function to update person-specific alphas 

ddm_pmwg_step <- function( c_alpha=NULL, c_logData=NULL, rti=NULL, xsi=NULL, K=NULL,
    c_MUa=NULL, c_SIGa=NULL, mu_hat=NULL, sigma_hat=NULL, args=NULL, weights=NULL, ddm_fun=NULL )
{
    
    #- collect some args:
    R       <- args$bayes_list$R
    pi_mix  <- args$bayes_list$pi_mix
    pi_mix1 <- args$bayes_list$pi_mix1
    pi_mix2 <- args$bayes_list$pi_mix2
    pi_mix3 <- args$bayes_list$pi_mix3
    epsilon <- args$bayes_list$epsilon

    #- check if we use the 3-component mixture:
    use_adapt <- !is.null(mu_hat) && !is.null(sigma_hat)

    #- fix first particle as current alpha:
    particles <- matrix( 0, nrow = R, ncol = K )
    particles[1,] <- c_alpha
    
    #- step 2: draw R-1 new particles from proposal:
    for ( r in 2:R ) {
        if ( use_adapt ) {
            #- 3-component mixture (S11 in Dao et al.):
            u <- runif(1)
            if ( u < pi_mix1 ) {
                particles[r,] <- mvtnorm::rmvnorm( 1, mu_hat, sigma_hat )
            } else if ( u < pi_mix1 + pi_mix2 ) {
                particles[r,] <- mvtnorm::rmvnorm( 1, c_alpha, sigma_hat )
            } else {
                particles[r,] <- mvtnorm::rmvnorm( 1, c_MUa, c_SIGa )
            }
        } else {
            #- 2-component mixture (burn-in phase, top of p. 28):
            if ( runif(1) < pi_mix ) {
                particles[r,] <- mvtnorm::rmvnorm( 1, c_alpha, epsilon * c_SIGa )
            } else {
                particles[r,] <- mvtnorm::rmvnorm( 1, c_MUa, c_SIGa )
            }
        }
    }
    
    #- step 3: compute unnormalized log-weights:
    log_liks    <- rep( 0, R )
    log_weights <- rep( 0, R )
    for ( r in 1:R ) {
        
        #- log-likelihood:
        if ( r == 1 ) {
            log_liks[r] <- c_logData
        } else {
            log_liks[r] <- ddm_fun( us=particles[r,], rt=rti, xs=xsi, args=args, weights=weights )
        }
        
        #- log-prior:
        log_prior <- mvtnorm::dmvnorm( particles[r,], c_MUa, c_SIGa, log=TRUE )
        
        #- log-proposal:
        if ( use_adapt ) {
            log_prop <- log(
                pi_mix1 * mvtnorm::dmvnorm( particles[r,], mu_hat, sigma_hat ) +
                pi_mix2 * mvtnorm::dmvnorm( particles[r,], c_alpha, sigma_hat ) +
                pi_mix3 * mvtnorm::dmvnorm( particles[r,], c_MUa, c_SIGa ) )
        } else {
            log_prop <- log( 
                pi_mix * mvtnorm::dmvnorm( particles[r,], c_alpha, epsilon*c_SIGa ) +
                (1-pi_mix) * mvtnorm::dmvnorm( particles[r,], c_MUa, c_SIGa ) )
        }

        log_prop <- log( 
            pi_mix * mvtnorm::dmvnorm( particles[r,], c_alpha, epsilon*c_SIGa ) +
            (1-pi_mix) * mvtnorm::dmvnorm( particles[r,], c_MUa, c_SIGa ) )
        
        log_weights[r] <- log_liks[r] + log_prior - log_prop
    }
    
    #- step 4: normalize weights (log-sum-exp trick für numerische Stabilität)
    log_weights <- log_weights - max(log_weights)
    particle_weights <- exp( log_weights )
    particle_weights <- particle_weights / sum( particle_weights )
    
    #- step 5: sample new index and return particle:
    k <- sample( 1:R, size=1, prob=particle_weights )
    
    return( list( 
        alpha_new = particles[k,],
        accepted  = (k != 1),  # TRUE wenn nicht das erste Particle gewählt wurde
        logData_new = log_liks[k]
    ) )
}

#- helper function for ddm_standard_step

ddm_logPrior <- function( alphai = NULL, MU = NULL, SIGMA = NULL ) 
{
    #- compute log-prior for MU:
    res <- mvtnorm::dmvnorm( alphai, MU, SIGMA, log = TRUE )
    return( sum( res ) )
}

ddm_generate_alphai <- function( c_alpha=NULL, sigma=NULL, MU=NULL, SIGMA=NULL, args=NULL) 
{
    
    if ( args$bayes_list$type_proposal == "mixture" ) {
        if ( runif(1) < args$pi_mix ) {
            new_alpha <- mvtnorm::rmvnorm( 1, c_alpha, args$epsilon * SIGMA )
        } else {
            new_alpha <- mvtnorm::rmvnorm( 1, MU, SIGMA )
        }
    } else {
        new_alpha <- mvtnorm::rmvnorm( 1, c_alpha, args$tau2 * sigma )
    }
    return( as.vector( new_alpha ) )
}

ddm_standard_step <- function( c_alpha=NULL, c_logData=NULL, rti=NULL, xsi=NULL, 
    K=NULL, c_MUa=NULL, c_SIGa=NULL, args=NULL, weights=NULL, ddm_fun=NULL )
{

    #- collect some args:
    pi_mix  <- args$bayes_list$pi_mix
    epsilon <- args$bayes_list$epsilon

    #- initialize output:
    out <- list( alpha_new = c_alpha, accepted = FALSE, logData_new = c_logData )

    #- step 1: compute logPrior for current alpha:
    c_logPrior <- ddm_logPrior( alphai=c_alpha, MU=c_MUa, SIGMA=c_SIGa )
    
    #- step 2: generate a proposal 
    p_alpha <- ddm_generate_alphai( c_alpha=c_alpha, sigma=diag(1,K),
        MU=c_MUa, SIGMA=c_SIGa, args=args )

    #- step 3: compute all things for posterior of the proposal:
    p_logPrior <- ddm_logPrior( alphai=p_alpha, MU=c_MUa, SIGMA=c_SIGa )
    p_logData  <- ddm_fun( us=p_alpha, rt=rti, xs=xsi, args=args, weights=weights )
    
    #- accept proposal?
    if ( args$bayes_list$type_proposal == "mixture" ) {
        log_q_forward  <- log( pi_mix * mvtnorm::dmvnorm( p_alpha, c_alpha, epsilon*c_SIGa) +
            (1-pi_mix) * mvtnorm::dmvnorm( p_alpha, c_MUa, c_SIGa) )
        log_q_backward <- log( pi_mix * mvtnorm::dmvnorm( c_alpha, p_alpha, epsilon*c_SIGa) +
            (1-pi_mix) * mvtnorm::dmvnorm( c_alpha, c_MUa, c_SIGa) )
        log_llratio <- ( p_logData + p_logPrior - log_q_forward  ) - 
                       ( c_logData + c_logPrior - log_q_backward )
    } else { # symmetrische Proposal, kürzt sich raus
        log_llratio <- ( p_logData + p_logPrior ) - ( c_logData + c_logPrior )
    }
    llratio = exp( min( log_llratio, 0 ) )
    if ( runif(1) < llratio ) {
        out$alpha_new   <- p_alpha
        out$accepted    <- TRUE
        out$logData_new <- p_logData 
    } 
    return( out )
}

#--- functions to update a and SIGAM_alpha

ddm_update_a <- function( a_d=NULL, nu0=NULL, c_invSIGa=NULL, K=NULL ) 
{
    #- update each a_d separately:
    for ( d in 1:K ) {
        #- shape and rate for IG:
        alpha_ig <- ( nu0 + K ) / 2
        beta_ig  <- nu0 * c_invSIGa[d,d] + 1 # 1/( a_d[d]^2 )
        #- sample from IG:
        a_d[d] <- 1 / rgamma( 1, shape=alpha_ig, rate=beta_ig )
    }
    return( a_d )
}

ddm_update_sigmaalpha <- function( alpha=NULL, MUa=NULL, I=NULL, K=NULL, 
    nu0=NULL, S0=NULL, a_d=NULL, type_proposal=NULL )
{

    #- compute covariance matrix of alpha:
    matMUa   <- matrix( MUa, I, K, byrow = TRUE )
    alphacov <- t( alpha - matMUa )%*%( alpha - matMUa )

    #- compute mean and scale for wishart distribution:
    if ( type_proposal == "pmwg" ) {
        k_alpha <- nu0 + K - 1 + I
        B_alpha <- 2*nu0*diag( 1/a_d, K ) + alphacov
    } else if ( type_proposal == "mixture" || type_proposal == "std" ) {
        k_alpha <- nu0 + K + I
        B_alpha <- S0 + alphacov
    }

    #- make the draw:
    c_invSIGa <- MCMCpack::rwish( k_alpha, base::solve( B_alpha ) )
    c_SIGa    <- base::solve( c_invSIGa )

    #- update a_d:
    if ( type_proposal == "pmwg" ) {
        a_d <- ddm_update_a( a_d=a_d, nu0=nu0, c_invSIGa=c_invSIGa, K=K )
    }
                
    #- output:
    out <- list( c_SIGa=c_SIGa, c_invSIGa=c_invSIGa, a_d=a_d )
    return( out )
}