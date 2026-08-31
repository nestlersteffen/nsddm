
#--- function to approximate log-likelihood function with Laplace

ddm7_getmodes <- function( parm = NULL, data_list = NULL, parm_list = NULL, 
    parm_table = NULL, args = NULL, check_t0 = FALSE )
{
  
    #- insert parm:
    parm_list <- ddm_include_free_parameters( parm = parm, 
        parm_list = parm_list, parm_table = parm_table )
  
    #- get BETA and SIGMA:
    parm_list <- ddm_compute_sigma( parm_list=parm_list, type_sigma=args$type_sigma )
    MU        <- parm_list[["MU"]]
    L         <- parm_list[["L"]]
    SIGMA     <- make_pd( parm_list$SIGMA, tol_factor = 1e-6 )
    invSIGMA  <- base::solve( SIGMA )
      
    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    info <- data_list[["info"]]

    #- load matrices in case of use_lan:
    if ( args$use_lan ) {
        lan_path <- paste0(system.file("extdata", package = "nbddm"), .Platform$file.sep)
        lan_load_weights_export( lan_path )
    }

    #- backend setup
    if ( is.null(args$maxlik_list$n_threads) || args$maxlik_list$n_threads <= 1 ) {
        `%loop%` <- foreach::`%do%`
    } else {
        cl <- parallel::makeCluster( args$maxlik_list$n_threads )
        doParallel::registerDoParallel( cl )
        on.exit( parallel::stopCluster(cl), add = TRUE )
        `%loop%` <- foreach::`%dopar%`
    }
  
    #- let's go:
    res <- foreach::foreach( i = seq(I), .packages = "nbddm" ) %loop% {    
      
        #- collect person data:
        idx1 <- info[i,3] - info[i,2] + 1 
        idx2 <- info[i,3]
        rti  <- rts[idx1:idx2]
        xsi  <- xs[idx1:idx2]

        #- a fix for t0:
        start_u  <- as.vector( MU )
        upper_t0 <- Inf
        if ( check_t0 ) {
            upper_t0 <- log( max( min(rti) - 1e-4, 0.05 ) )
            start_u[4] <- min( start_u[4], upper_t0 - 1e-4 )
        }

        if ( args$use_tmb_modes & !( args$use_lan ) ) {

            #- optimize:
            tdata <- list( rts = rti, xs = xsi, 
                pts = c(-0.973906528517172, -0.865063366688985, -0.679409568299024,
                        -0.433395394129247, -0.148874338981631,  0.148874338981631,
                         0.433395394129247,  0.679409568299024,  0.865063366688985, 0.973906528517172),
                wgh = c(0.066671344308688, 0.149451349150581, 0.219086362515982,
                        0.269266719309996, 0.295524224714753, 0.295524224714753,
                        0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688),
                MU = MU, SIGMA = SIGMA )
            tdata <- c( tdata, model = "ddm7_getmodes_tmb")
            tparm <- list( us = start_u )

            obj <- TMB::MakeADFun( data = tdata, parameters = tparm, 
                DLL = "nbddm_TMBExports", # DLL = "ddm4_getmodes_tmb", 
                silent = TRUE )

            opt <- suppressWarnings( tryCatch( 
                stats::nlminb( obj$par, objective = obj$fn, gradient = obj$gr, 
                    control = list( trace = 0, step.min = 1 ),
                    #upper = c(Inf, Inf, Inf, upper_t0) 
                    ),
                error   = function(e) { NULL },
                warning = function(w) { NULL }
            ) )

        } else {

            #- optimize
            opt <- suppressWarnings( tryCatch( { stats::nlminb( start_u, 
                objective = ddm7_random_nllfct_rcpp, control = list( trace = 0 ),
                #upper = c(Inf, Inf, Inf, upper_t0, Inf, Inf, Inf ), 
                rt = rti, xs = xsi, MU = MU, SIGMA = SIGMA, 
                type_alpha=args$type_alpha, type_ddm=args$type_ddm, 
                kmax=args$kmax, delta=args$delta, 
                use_lan=args$use_lan ) },
                error = function(e) { NULL } ) )

        }

        #- we check convergence etc.
        opt_ok <- !is.null( opt ) && opt$objective != 1e+06 && opt$convergence == 0

        # converged...
        if ( opt_ok ) {

            #- save modes and compute other stuff:
            MODE  <- opt$par
            if ( args$use_tmb_modes & !( args$use_lan ) ) {
                HESS <- tryCatch( obj$he(), error = function(e) NULL )
            } else {
                HESS <- tryCatch(
                    numDeriv::hessian( x = MODE, func = ddm7_random_nllfct_rcpp, 
                        rt = rti, xs = xsi, MU = MU, SIGMA = SIGMA, 
                        type_alpha=args$type_alpha, type_ddm=args$type_ddm, 
                        kmax=args$kmax, delta=args$delta, use_lan=args$use_lan,
                        method = "Richardson", method.args=list( d=0.0001 ) ),
                    error = function(e) NULL
                )
            }

            hess_fb <- is.null( HESS ) | any( is.na( HESS ) ) | any( is.infinite( HESS ) )
            if ( hess_fb ) { 

                #- save modes and compute other stuff:
                MODE <- start_u
                HESS <- invSIGMA 
            
            }
        
        } else { 

            #- save modes and compute other stuff:
            MODE <- start_u
            HESS <- invSIGMA 
            hess_fb <- TRUE

        }

        #- compute AGH stuff...
        HESS_pd   <- make_pd( HESS, tol_factor = 1e-3 )
        CHOL_HESS <- base::chol( HESS_pd )
        iHESS     <- chol2inv( CHOL_HESS )
        CHOL      <- base::chol( iHESS )
        DET       <- 1 / prod( diag( CHOL_HESS ) )^2

        #- save results:
        return( list( MODE = MODE, HESS = HESS, iHESS = iHESS, CHOL = CHOL, DET = DET, converged = opt_ok, nothess_fb = !hess_fb ) )

    }
  
    #- output:
    return( res )
}