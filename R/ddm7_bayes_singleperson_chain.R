
#- -------------------------- MAIN

ddm7_bayes_singleperson_chain <- function( inits=NULL, rt=NULL,xs=NULL, args=NULL, 
    muPrior_sp=rep(0,7), sdPrior_sp=rep(1,7), verbose=TRUE )
{
      
    #- load weights in case we need them:
    if ( args$use_lan ) {
        ddm_load_weights_rcpp(ddm="seven")
    }

    #- set args for the MCMC chain:
    accept_A <- accept_B <- 0
    accept_post_A <- accept_post_B <- 0
    biter  <- args$bayes_list$biter
    nchain <- args$bayes_list$nchain
    burnin <- args$bayes_list$burnin

    #- make parms-matrix to save results
    parms    <- matrix( 0, ncol=7, nrow=biter)
    parms[1,] <- inits

    #- compute all things for current estimate
    c_logPrior <- sum( dnorm( parms[1,], muPrior_sp, sdPrior_sp, log = TRUE ) )
    c_logData  <- -1*ddm7_nllfct_export( parms[1,], rt, xs, args$type_alpha, args$type_ddm, 
        args$kmax, args$delta, args$use_lan )

    #- sampling:
    for ( nn in 2:biter ) {
  
        #- life signal:
        if ( verbose & nn%%100 == 0 ) {
            print( paste0( "Iteration: ", nn ) )
        }

        #- current parm:
        current  <- parms[nn-1,]

        # ----- BLOCK A: the standard parms
        proposal_A <- current
        proposal_A[1:4] <- as.numeric( mvtnorm::rmvnorm( 1, current[1:4], diag( args$bayes_list$tau2_A, 4 ) ) )
                            
        #- compute logPrior and logData for proposal:
        p_logPrior_A <- sum( dnorm( proposal_A, muPrior_sp, sdPrior_sp, log = TRUE ) )
        p_logData_A  <- -1*ddm7_nllfct_export( proposal_A, rt, xs, args$type_alpha, args$type_ddm, 
        args$kmax, args$delta, args$use_lan )
            
        #- accept proposal?
        log_llratio_A <- ( p_logData_A + p_logPrior_A ) - ( c_logData + c_logPrior )
        llratio_A <- exp( min( log_llratio_A, 0 ) )
        if ( runif(1) < llratio_A ) {
            current  <- proposal_A
            accept_A <- accept_A + 1 
            if ( nn > burnin ) accept_post_A <- accept_post_A + 1
            c_logPrior <- p_logPrior_A
            c_logData  <- p_logData_A
        } 
        
        # ---- BLOCK B: (sv, sz, sτ) = indices 5:7 ----
        proposal_B <- current
        proposal_B[5:7] <- as.numeric( mvtnorm::rmvnorm( 1, current[5:7], diag( args$bayes_list$tau2_B, 3 ) ) )
            
        #- compute logPrior and logData for proposal:
        p_logPrior_B <- sum( dnorm( proposal_B, muPrior_sp, sdPrior_sp, log = TRUE ) )
        p_logData_B  <- -1*ddm7_nllfct_export( proposal_B, rt, xs, args$type_alpha, args$type_ddm, 
            args$kmax, args$delta, args$use_lan )
            
        #- accept proposal?
        log_llratio_B <- ( p_logData_B + p_logPrior_B ) - ( c_logData + c_logPrior )
        llratio_B <- exp( min( log_llratio_B, 0 ) )
        if ( runif(1) < llratio_B  ) {
            current    <- proposal_B
            c_logPrior <- p_logPrior_B
            c_logData  <- p_logData_B
            accept_B   <- accept_B + 1
            if ( nn > burnin ) accept_post_B <- accept_post_B + 1
        }

        #- store current state (after both block updates):
        parms[nn,] <- current

    }

    return( list( parms=parms, accept_A=accept_A, accept_B=accept_B,
        accept_post_A=accept_post_A, accept_post_B = accept_post_B ) )

}

# ---------- end 