
#- -------------------------- MAIN

ddm_random_bayes <- function( data_list=NULL, parm_table=NULL, args=NULL, 
    ddm="four",verbose=TRUE )
{
    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    K    <- data_list[["K"]]
    info <- data_list[["info"]]
      
    #- get initial values for the chains:
    inits  <- ddm_inits_bayes( rts=rts, I=I, K=K, args=args, ddm=ddm ) 
    m      <- inits$m 
    invM   <- base::solve( inits$M )
    nu0    <- inits$nu0 
    S0     <- inits$S0
    a_d    <- inits$a_d

    #- set args for the MCMC chain:
    accept <- accept_post <- 0
    biter  <- args$bayes_list$biter
    nchain <- args$bayes_list$nchain
    burnin <- args$bayes_list$burnin

    #- sampling:
    if ( nchain == 1 ) {

        if ( !args$use_rcpp ) {

            #- load weights in case we need them:
            if ( args$use_lan ) {
                if ( args$use_tf ) {
                    dnn <- ddm_load_dnn( ddm=ddm ) 
                } else {
                    weights <- ddm_load_weights( ddm=ddm )
                }
            }

            #- load function for the correct ddm:
            ddm_fun <- switch( ddm,
                "four"  = ddm4_logLdata,
                "seven" = ddm7_logLdata,
                stop( paste( "The 4- or 7-parameter DDM can be estimated." ) )
            )

            #- load function to make the step:
            ddm_step <- switch( args$bayes_list$type_proposal,
                "pmwg"     = ddm_pmwg_step,
                "mixture"  = ddm_standard_step,
                stop( paste( "Unknown type of proposal generation." ) )
            )    

            #- indices to save results:
            Sinds  <- which( lower.tri( S0, diag =TRUE ) )

            #- matrices to collect the results:
            MUa       <- matrix( 0, ncol = K, nrow = biter )
            MUa[1,]   <- inits$MUa[,1]
            SIGa      <- matrix( 0, ncol = K*(K+1)/2, nrow = biter )
            c_SIGa    <- inits$SIGa[,,1]
            SIGa[1,]  <- c_SIGa[Sinds]
            c_invSIGa <- base::solve( c_SIGa )

            #- matrix with initial random effects:
            alpha <- inits$alpha[,,1]

            #- matrix to store alpha draws after burnin for adaptation:
            if ( args$bayes_list$use_adapt_dao ) {
                #- matrix to store alpha draw
                alpha_store <- array( 0, dim=c( biter - burnin, I, K) )
                #- matrices to store current mu_hat and sigma_hat per person:
                mu_hats    <- matrix( 0, nrow = I, ncol = K )
                sigma_hats <- array( diag(K), dim = c( K, K, I ) )
                use_adapt  <- FALSE
            } else {
                mu_hat    <- NULL
                sigma_hat <- NULL
            }
            
            #- we pre-compute single log-lik values for current alphas:
            c_logData <- rep( 0, I )
            idx <- 0
            for ( i in seq( I ) ) {
                #- get data:
                ni  <- info[i,2]
                rti <- rts[(idx+1):(idx+ni)]
                xsi <- xs[(idx+1):(idx+ni)]
                #- compute data densities:
                c_logData[i] <- ddm_fun( us=alpha[i,], rt=rti, xs=xsi, args=args, weights=weights, dnn=dnn )
                #- increase idx:
                idx <- idx + ni
            }

            #- let's go and sample:
            for ( nn in 2:biter ) {
              
                #- life signal:
                if ( verbose & nn%%10 == 0 ) {
                    print( paste0( "Iteration: ", nn ) )
                }

                #- update MUalpha:
                cMUa <- base::solve( invM + I*c_invSIGa )
                eMUa <- cMUa %*% ( invM%*%m + c_invSIGa%*%colSums( alpha ) )
                MUa[nn,] <- mvtnorm::rmvnorm( 1, eMUa, cMUa ) 

                #- update SIGMAalpha: 
                tmp_SIGa <- ddm_update_sigmaalpha( alpha=alpha, MUa=MUa[nn,], I=I, K=K, 
                    nu0=nu0, S0=S0, a_d=a_d, type_sigma_prior=args$bayes_list$type_sigma_prior )
                a_d       <- tmp_SIGa$a_d
                c_SIGa    <- tmp_SIGa$c_SIGa
                c_invSIGa <- tmp_SIGa$c_invSIGa
                SIGa[nn,] <- c_SIGa[ Sinds ]
                
                if ( args$bayes_list$use_adapt_dao ) {
                    #- check if we are in adaptation phase:
                    n_adapt  <- nn - burnin
                    in_adapt <- n_adapt > 0
                    #- update mu_hat and sigma_hat every 20 iterations, but stop after 5000:
                    if ( in_adapt && n_adapt > K + 1 && n_adapt %% 20 == 0 && nn < 5000 ) {
                        for ( i in seq(I) ) {
                            mu_hats[i,]      <- colMeans( alpha_store[1:n_adapt, i, ] )
                            sigma_hats[,,i]  <- cov( alpha_store[1:n_adapt, i, ] )
                        }
                        use_adapt <- TRUE
                    }
                }

                #- update alpha:
                idx <- 0
                for ( i in seq( I ) ) {
                
                    #- get the data:
                    ni  <- info[i,2]
                    rti <- rts[(idx+1):(idx+ni)]
                    xsi <- xs[(idx+1):(idx+ni)]

                    #- get current mu_hat and sigma_hat for person i:
                    if ( args$bayes_list$use_adapt_dao ) {
                        if ( use_adapt ) {
                            mu_hat    <- mu_hats[i,]
                            sigma_hat <- sigma_hats[,,i]
                        } else {
                            mu_hat    <- NULL
                            sigma_hat <- NULL
                        }
                    }
          
                    #- make a standard-step or a pmwg-step
                    res_i <- ddm_step( c_alpha=alpha[i,], c_logData=c_logData[i], 
                        rti=rti, xs=xsi, K=K, c_MUa=MUa[nn,], c_SIGa=c_SIGa, 
                        mu_hat=mu_hat, sigma_hat=sigma_hat,
                        args=args, weights=weights, dnn=dnn, ddm_fun=ddm_fun )

                    #- store draw for adaptation after burnin:
                    if ( args$bayes_list$use_adapt_dao ) {
                        if ( in_adapt ) {
                            alpha_store[n_adapt, i, ] <- res_i$alpha_new
                        }
                    }
                    
                    #- save results:
                    alpha[i,]    <- res_i$alpha_new
                    c_logData[i] <- res_i$logData_new
                    if ( res_i$accepted ) { 
                        accept <- accept + 1
                        if ( nn > burnin ) accept_post <- accept_post + 1
                    }

                    #- increase idx:
                    idx <- idx + ni

                } # i
              
            } # biter 

        } else {

            #- get path:
            lan_path <- if ( args$use_lan && args$use_rcpp ) {
                paste0( system.file( paste0("extdata/", ddm ), package = "nsddm" ), .Platform$file.sep )
            } else {""}

            res <- ddm_chain_rcpp( rts=rts, xs=xs, info=info, I=I, K=K, 
                alpha=inits$alpha[,,1],c_MUa=inits$MUa[,1],c_SIGa=inits$SIGa[,,1],
                m=m, invM=invM, S0=S0, nu0=nu0, a_d=a_d,
                biter=biter, burnin=burnin, use_adapt_dao=args$bayes_list$use_adapt_dao,
                tau2=args$bayes_list$tau2, epsilon=args$bayes_list$epsilon, 
                pi_mix=args$bayes_list$pi_mix, pi_mix1=args$bayes_list$pi_mix1, pi_mix2=args$bayes_list$pi_mix2, 
                pi_mix3=args$bayes_list$pi_mix3, R=args$bayes_list$R,
                method=args$type_ddm, type_alpha=args$type_alpha, 
                type_proposal=args$bayes_list$type_proposal,
                type_sigma_prior=args$bayes_list$type_sigma_prior, 
                verbose=verbose, ddm=ddm, kmax=args$kmax, delta=args$delta,
                use_lan=args$use_lan, lan_path=lan_path )

        } 

    } else {

        res <- ddm_parallel_bayes( rts=rts, xs=xs, info=info, I=I, K=K,
            inits=inits, ddm=ddm, args=args )
        
    }
    
    #- ----------------
    #-   final things
    
    #- step 1: combine chains
    if ( nchain > 1 ) {

        #- combine MUa and SIGa in one array
        parms_chains <- array( 0, dim = c( biter-burnin, K + K*(K+1)/2, nchain ) )
        for ( nc in 1:nchain ) {
            parms_chains[,,nc] <- cbind( res[[nc]]$MUa[-c(1:burnin),], res[[nc]]$SIGa[-c(1:burnin),] )
        }

        #- compute convergence diagnostics
        rhat <- ddm_rhat( parms_chains )
        iact <- matrix( 0, nrow = nchain, ncol = K + K*(K+1)/2 )
        for ( nc in 1:nchain ) {
            iact[nc,]  <- apply( parms_chains[,,nc],  2, ddm_iact )
        }
        iact <- colMeans( iact )
        conv <- data.frame( iact = iact, rhat = rhat )

        #- combine chains
        parms <- do.call( rbind, lapply(1:nchain, function(nc) parms_chains[,,nc]) )

        #- alpha: Mittelwert über Ketten:
        alpha <- (1/nchain)*Reduce( "+", lapply(1:nchain, function(nc) res[[nc]]$alpha ) )
    
        #- accept:
        p_accept      <- mean( sapply(1:nchain, function(nc) res[[nc]]$accept) / (biter*I) )
        p_accept_post <- mean( sapply(1:nchain, function(nc) res[[nc]]$accept_post) / ((biter-burnin)*I) )

    } else {

        if ( args$use_rcpp ) {
            MUa         <- res$MUa
            SIGa        <- res$SIGa
            accept      <- res$accept
            accept_post <- res$accept_post
            alpha       <- res$alpha
        }
        
        #- combine matrices
        parms <- cbind( MUa, SIGa )
        parms <- parms[-c(1:burnin),]
        conv  <- data.frame( iact = apply( parms, 2, ddm_iact ) )

        #- proportion of accepted proposals:
        p_accept <- accept/(biter*I)
        p_accept_post <- accept_post / ( (biter - burnin) * I )

    }

    #- now compute posterior means, standard deviations and so on:
    parm_table$M    <- apply( parms, 2, mean )
    parm_table$Md   <- apply( parms, 2, median )
    parm_table$pSD  <- apply( parms, 2, sd )
    qu.parms        <- apply( parms, 2, quantile, probs = c( 0.025,0.975 ) )
    parm_table$lCI  <- qu.parms[1,]
    parm_table$uCI  <- qu.parms[2,]
    parm_table <- cbind( parm_table, conv )

    #- -----------------
    #- output:
    res <- list( parm_table=parm_table, p_accept=p_accept, p_accept_post=p_accept_post, 
        alpha=as.data.frame( alpha ), chains=parms )
    return( res )

}