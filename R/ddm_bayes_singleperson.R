
#----------------------- MAIN

ddm_bayes_singleperson <- function( rt=NULL, xs=NULL, parms=NULL, args=NULL, ddm="four",
	muPrior_sp=NULL, sdPrior_sp=NULL, inits=NULL, cl=NULL, verbose=TRUE ) 
{

    #- get K:
    if ( ddm=="four") K <- 4 else K <- 7

    #- some checks and initialize priors:
    if ( !is.null( muPrior_sp ) ) {
        if ( length( muPrior_sp ) != K ) stop( "Length of init vector have to match ddm")
    } else {
        muPrior_sp <- rep(0,K)
    }
    if ( !is.null( sdPrior_sp ) ) {
        if ( length( sdPrior_sp ) != K ) stop( "Length of init vector have to match ddm")
    } else {
        sdPrior_sp <- rep(1,K)
    }
    if ( !is.null(inits) ) {
        if ( length(inits) != K ) stop( "Length of init vector have to match ddm")
    }

    #- collect chain args:
    nchain <- args$bayes_list$nchain
    biter  <- args$bayes_list$biter
    burnin <- args$bayes_list$burnin

    #- initialize chains:
    init_mat <- matrix( 0, ncol=K, nrow=nchain )
    for ( nc in 1:nchain ) {
        if ( is.null( inits ) ) {
            init_mat[nc,] <- ddm_inits_mu( rts=rt, type_alpha=args$type_alpha, ddm=ddm )
        } else {
            init_mat[nc,] <- inits + rnorm(K,0,0.01) # add a little bit of noise...
        }
    }

    #- select function depending on ddm:
    ddm_fun <- switch( ddm,
        "four"  = ddm4_bayes_singleperson_chain,
        "seven" = ddm7_bayes_singleperson_chain,
        stop( paste( "The 4- or 7-parameter DDM can be estimated." ) )
    )

    #- start to sample:
    if ( nchain == 1 ) {

    	results <- ddm_fun( inits=init_mat[1,], rt=rt, xs=xs, args=args, verbose=verbose, 
            muPrior_sp=muPrior_sp, sdPrior_sp=sdPrior_sp )
    	results <- list( results )

    } else {

    	#- backend setup
	    if (is.null(cl)) {
	        cl <- parallel::makeCluster(nchain)
	        doParallel::registerDoParallel(cl)
	        on.exit(parallel::stopCluster(cl), add = TRUE)
	       
            parallel::clusterEvalQ(cl, {
                RhpcBLASctl::blas_set_num_threads(1)
                RhpcBLASctl::omp_set_num_threads(1)
            })

        } else {
	        doParallel::registerDoParallel(cl)
	    }

	    #- run chains in parallel:
	    nc <- NULL
	    results <- foreach::foreach( nc = seq(1, nchain, 1 ), .packages = "nsddm" ) %dopar% {
        
        	ddm_fun( inits=init_mat[nc,], rt=rt, xs=xs, args=args, 
    			verbose=FALSE, muPrior_sp=muPrior_sp, sdPrior_sp=sdPrior_sp )
    	} 
        
    }

    #- fill parms and make a cube with deleted burin:
    parms <- array( 0, dim = c( biter, K, nchain ) )
    parms_final <- array( 0, dim = c( biter-burnin, K, nchain ) )
    for ( nc in 1:nchain ) {
    	parms[,,nc]       <- results[[nc]]$parms
        parms_final[,,nc] <- parms[-c(1:burnin),,nc]
    }

    #- compute convergence diagnostics:
    iact <- matrix( 0, nrow = nchain, ncol = K )
    for ( nc in 1:nchain ) {
        iact[nc,]  <- apply( parms_final[,,nc], 2, ddm_iact )
    }
    iact <- colMeans( iact )

    rhat <- NULL
    if ( nchain > 1 ) {
        rhat <- ddm_rhat( parms_final )
    } 

	#- compute acceptance probs:
    if ( K == 4 ) {
        p_accept_all  <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept) / (biter) )
        p_accept_post <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept_post) / ((biter-burnin)) )
        p_accept <- c(p_accept_all, p_accept_post)
        names( p_accept ) <- c("all","post_burnin")
    } else {
        p_accept_all_A  <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept_A) / (biter) )
        p_accept_post_A <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept_post_A) / ((biter-burnin)) )
        p_accept_all_B  <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept_B) / (biter) )
        p_accept_post_B <- mean( sapply(1:nchain, function(nc) results[[nc]]$accept_post_B) / ((biter-burnin)) )
        p_accept <- c( p_accept_all_A, p_accept_post_A, p_accept_all_B, p_accept_post_B )
        names( p_accept ) <- c("all_A","post_burnin_A", "all_B", "post_burnin_B" )
    }
        
    #- now combine the chains and make a nice table:
    parms <- do.call( rbind, lapply( 1:nchain, function(nc) parms_final[,,nc] ) )
    postM  <- apply( parms, 2, mean )
    postMd <- apply( parms, 2, median )
    postSD <- apply( parms, 2, sd )
    qu.parms <- apply( parms, 2, quantile, probs = c( 0.025,0.975 ) )
    tab <- data.frame( Par=paste0("alpha",c(1:K)), 
        M=postM, Md=postMd, SD=postSD, lCI=qu.parms[1,], uCI=qu.parms[2,],
        iact = iact )
    if ( !is.null( rhat ) ) tab$rhat <- rhat
    
    #- output:
    out <- list( parm_table=tab, p_accept=p_accept, chains=parms )
    return( out )

}

