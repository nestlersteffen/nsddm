
#- -------------------------- MAIN

ddm4_bayes_singleperson_chain <- function( inits=NULL, rt=NULL,xs=NULL,args=NULL,
    muPrior_sp =c(0,0,0,0), sdPrior_sp=c(1,1,1,1), verbose=TRUE )
{
  
    #- load weights in case we need them:
    if ( args$use_lan ) {
        ddm_load_weights_rcpp(ddm="four")
    }

    #- set args for the MCMC chain:
    accept <- accept_post <- 0
    biter  <- args$bayes_list$biter
    nchain <- args$bayes_list$nchain
    burnin <- args$bayes_list$burnin

    #- initialize chains:
    parms     <- matrix( 0, ncol=4, nrow=biter)
    parms[1,] <- inits
    
    #- compute all things for the first estimate
    c_logPrior <- sum( dnorm( parms[1,], muPrior_sp, sdPrior_sp, log = TRUE ) )
    c_logData  <- -1*ddm4_nllfct_export( parms[1,], rt, xs, args$type_alpha, args$type_ddm, 
        args$kmax, args$delta, args$use_lan )

    #- sampling:
    for ( nn in 2:biter ) {
  
        #- life signal:
        if ( verbose & nn%%100 == 0 ) {
            print( paste0( "Iteration: ", nn ) )
        }

        #- generate a proposal:
        current  <- parms[nn-1,]
        proposal <- as.numeric( mvtnorm::rmvnorm( 1, current, diag( args$bayes_list$tau2_A, 4 ) ) )
                            
        #- compute logPrior and logData for proposal:
        p_logPrior <- sum( dnorm( proposal, muPrior_sp, sdPrior_sp, log = TRUE ) )
        p_logData  <- -1*ddm4_nllfct_export( proposal, rt, xs, args$type_alpha, args$type_ddm, 
        args$kmax, args$delta, args$use_lan )
            
        #- accept proposal?
        log_llratio <- ( p_logData + p_logPrior ) - ( c_logData + c_logPrior )
        llratio <- exp( min( log_llratio, 0 ) )
        if ( runif(1) < llratio ) {
            parms[nn,] <- proposal
            accept       <- accept + 1 
            if ( nn > burnin ) accept_post <- accept_post + 1
            c_logPrior   <- p_logPrior
            c_logData    <- p_logData
        } else {
            parms[nn,] <- current
        }

    }

    #- output:
    out <- list( parms=parms, accept=accept, accept_post=accept_post )
    return( out )

}

#--------- end