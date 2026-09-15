
#--- function to obtain some start values for maximum likelihood estimation

ddm7_getstarts <- function( data_list=NULL, parm_table=NULL, args=NULL, 
    check_t0=FALSE, ddm=NULL, n_threads=1L )
{
  
    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    info <- data_list[["info"]]

    #- initial values for MU
    parm <- ddm_inits_mu( rts=rts, type_alpha=args$type_alpha, ddm="seven" )
    print( parm )

    #- load matrices in case of use_lan:
    if ( args$use_lan ) {
        lan_path <- paste0(system.file("extdata", package = "nsddm"), .Platform$file.sep)
        lan_load_weights_export( lan_path )
    }
  
    #- backend setup
    if ( is.null(n_threads) || n_threads <= 1 ) {
        `%loop%` <- foreach::`%do%`
    } else {
        cl <- parallel::makeCluster( n_threads )
        doParallel::registerDoParallel( cl )
        on.exit( parallel::stopCluster(cl), add = TRUE )
        `%loop%` <- foreach::`%dopar%`
    }
  
    #- let's go:
    res <- foreach::foreach( i = seq(I), .packages = "nsddm" ) %loop% {    
      
        #- collect person data:
        idx1 <- info[i,3] - info[i,2] + 1 
        idx2 <- info[i,3]
        rti  <- rts[idx1:idx2]
        xsi  <- xs[idx1:idx2]
        
        if ( args$use_tmb_modes & !( args$use_lan ) ) {

            #- optimize:
            tdata <- list( rts = rti, xs = xsi, 
                pts = c(-0.973906528517172, -0.865063366688985, -0.679409568299024,
                        -0.433395394129247, -0.148874338981631,  0.148874338981631,
                         0.433395394129247,  0.679409568299024,  0.865063366688985, 0.973906528517172),
                wgh = c(0.066671344308688, 0.149451349150581, 0.219086362515982,
                        0.269266719309996, 0.295524224714753, 0.295524224714753,
                        0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688) )
            tdata <- c( tdata, model = "ddm7_tmb")
            tparm <- list( parm = parm )

            obj   <- TMB::MakeADFun( data = tdata, parameters = tparm,  
                DLL = "nsddm_TMBExports", # DLL = "ddm4_tmb",
                silent = TRUE )
            fit_n <- suppressWarnings( tryCatch( 
                stats::nlminb( obj$par, objective = obj$fn, gradient = obj$gr, 
                control = list( trace = 0, step.min = 1 ) ),
                    error   = function(e) { NULL },
                    warning = function(w) { NULL }
            ) )

        } else {

            #- optimize
            fit_n <- suppressWarnings( tryCatch( { stats::nlminb( parm, 
                objective = ddm7_nllfct_export, control = list( trace = 0 ),
                # upper = c(Inf, Inf, Inf, upper_t0), 
                rt = rti, xs = xsi, type_alpha=args$type_alpha, type_ddm=args$type_ddm, 
                kmax=args$kmax, delta=args$delta, use_lan=args$use_lan ) },
                error = function(e) { NULL } ) )

        }

        #- error handling:
        if ( is.null( fit_n )  ) {
            us <- parm
        } else {
            
            if ( check_t0 ) {
                upper_t0 <- log( max( min( rti ) - 1e-4, 0.05 ) )
                fit_n$par[4] <- min( fit_n$par[4], upper_t0 - 1e-4)
            }
            us <- fit_n$par

        }    

        return( us )

    }

    #- make us-object:
    us <- do.call("rbind", res )

    #- we truncate very high values:
    us <- apply( us, 2, function( x ) { 
        bounds <- quantile( x, prob = c( 0.15, 0.85 ) )
        return( pmax( bounds[1], pmin( bounds[2], x ) ) ) 
    } )

    #- compute some relevant things:
    START_MU    <- colMeans( us )
    START_SIGMA <- cov( us )

    if ( args$type_sigma == "sigma" ) {
        starts <- c( START_MU, START_SIGMA[ lower.tri( START_SIGMA, diag = TRUE ) ] )
    } else if (args$type_sigma == "cholesky" ) {
        START_L <- t( chol( START_SIGMA ) )
        starts  <- c( START_MU, START_L[ lower.tri( START_L, diag = TRUE ) ] )
    }

    #- output:
    return( list( starts = starts, us = us ) )
}