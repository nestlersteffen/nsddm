
#- -------------------------- MAIN

ddm_parallel_bayes <- function( rts=NULL, xs=NULL, info=NULL, I=NULL, K=NULL,
    inits=NULL, ddm=NULL, args=NULL, cl=NULL ) 
{

    #- collect args:
    nchain <- args$bayes_list$nchain

    #- get path:
    lan_path <- if ( args$use_lan && args$use_rcpp ) {
        paste0( system.file( paste0("extdata/", ddm ), package = "nsddm" ), .Platform$file.sep )
    } else {""}

    #- backend setup
    if (is.null(cl)) {
        cl <- parallel::makeCluster(nchain)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
        doParallel::registerDoParallel(cl)
    }
        
    #- run chains in parallel:
    nc <- NULL
    results <- foreach::foreach( nc = seq(1, nchain, 1 ), .packages = "nsddm" ) %dopar% {
        
        nbddm::ddm_chain_rcpp( rts=rts, xs=xs, info=info, I=I, K=K, 
            alpha=inits$alpha[,,nc],c_MUa=inits$MUa[,nc],c_SIGa=inits$SIGa[,,nc],
            m=inits$m, invM=base::solve(inits$M), S0=inits$S0, nu0=inits$nu0, a_d=inits$a_d,
            biter=args$bayes_list$biter, burnin=args$bayes_list$burnin, use_adapt_dao=args$bayes_list$use_adapt_dao,
            tau2=args$bayes_list$tau2, epsilon=args$bayes_list$epsilon, 
            pi_mix=args$bayes_list$pi_mix, pi_mix1=args$bayes_list$pi_mix1, 
            pi_mix2=args$bayes_list$pi_mix2, pi_mix3=args$bayes_list$pi_mix3, R=args$bayes_list$R,
            method=args$type_ddm, type_alpha=args$type_alpha, 
            type_proposal=args$bayes_list$type_proposal, 
            verbose=FALSE, ddm=ddm, kmax=args$kmax, delta=args$delta,
            use_lan=args$use_lan, lan_path=lan_path )

    }     
    
    return( results )
}