
#--- function to obtain some start values for maximum likelihood estimation

ddm4_getstarts <- function( data_list=NULL, parm_table=NULL, args=NULL, 
    check_t0=FALSE, ddm=NULL )
{
  
    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    info <- data_list[["info"]]

    #- initial values for MU
    #parm <- ddm_inits_mu( rts=rts, type_alpha=args$type_alpha, ddm="four" )
    parm <- c( 0, log( 0.5 ), log(0.5), log( min(rts ) * 0.5 ) )

    #- make matrix to save values
    us <- matrix( 0, , ncol = length( parm ), nrow = I )

    #- load matrices in case of use_lan:
    if ( args$use_lan ) {
        lan_path <- paste0( system.file("extdata/four", package = "nbddm"), .Platform$file.sep )
        lan_load_weights_export( lan_path, "four" )
    }
  
    #- let's go:
    for ( i in 1:I ) {
    
        #- collect person data:
        idx1 <- info[i,3] - info[i,2] + 1 
        idx2 <- info[i,3]
        rti  <- rts[idx1:idx2]
        xsi  <- xs[idx1:idx2]
        
        if ( args$use_tmb_modes & !( args$use_lan ) ) {

            #- optimize:
            tdata <- list( rts = rti, xs = xsi, method = 0, type_alpha = 0 )
            tdata <- c( tdata, model = "ddm4_tmb")
            tparm <- list( parm = parm )

            obj   <- TMB::MakeADFun( data = tdata, parameters = tparm,  
                DLL = "nbddm_TMBExports", # DLL = "ddm4_tmb",
                silent = TRUE )
            fit_n <- suppressWarnings( tryCatch( 
                stats::nlminb( obj$par, objective = obj$fn, gradient = obj$gr, 
                control = list( trace = 0, step.min = 1 ) ),
                    error   = function(e) { NULL },
                    warning = function(w) { NULL }
            ) )

        } else if ( args$use_lan ) {

            #- make fit object
            obj <- ddm_MakeOptFun( ddm4_lanll_grad_weights_rcpp, rt=rti, xs=xsi )

            #- optimize
            fit_n <- suppressWarnings( tryCatch( { stats::nlminb( parm, 
                objective=obj$fn, gradient=obj$gr, control = list( trace = 0 ) ) },
                error = function(e) { NULL } ) )

        } else {

            #- optimize
            fit_n <- suppressWarnings( tryCatch( { stats::nlminb( parm, 
                objective = ddm4_nllfct_export, control = list( trace = 0 ),
                # upper = c(Inf, Inf, Inf, upper_t0), 
                rt = rti, xs = xsi, type_alpha=args$type_alpha, type_ddm=args$type_ddm, 
                kmax=args$kmax, delta=args$delta, use_lan=args$use_lan ) },
                error = function(e) { NULL } ) )

        }

        #- error handling:
        if ( is.null( fit_n )  ) {
            us[i,] <- parm
        } else {
            
            if ( check_t0 ) {
                upper_t0 <- log( max( min( rti ) - 1e-4, 0.05 ) )
                fit_n$par[4] <- min( fit_n$par[4], upper_t0 - 1e-4)
            }
            us[i,] <- fit_n$par

        }    

    }

    #- we truncate very high values:
    us <- apply( us, 2, function( x ) { 
        bounds <- quantile( x, prob = c( 0.05, 0.95 ) )
        return( pmax( bounds[1], pmin( bounds[2], x ) ) ) 
    } )

    #- compute some relevant things:
    START_MU    <- colMeans( us )
    START_SIGMA <- cov( us )
    
    if ( ddm == "seven" ) {
        #- deterministic plausible defaults for the s-parameters in ML
        sv  <- max( 0.8 * abs( START_MU[1] ), 0.1 )   # 80% der mittleren Drift, min 0.1
        sz  <- 0.05                                   # konservativ klein
        st0 <- 0.3 * exp( START_MU[4] )               # 30% der mittleren t0
        #- add them
        add_MU    <- c( log(sv), log(sz), log(st0) )
        START_MU  <- c( START_MU, add_MU )
        add_SIGMA <- diag(c(0,0,0,0,0.05,0.05,0.05),7)
        add_SIGMA[1:4,1:4] <- START_SIGMA
        START_SIGMA <- add_SIGMA
    }

    if ( args$type_sigma == "sigma" ) {
        starts <- c( START_MU, START_SIGMA[ lower.tri( START_SIGMA, diag = TRUE ) ] )
    } else if (args$type_sigma == "cholesky" ) {
        START_L <- t( chol( START_SIGMA ) )
        starts  <- c( START_MU, START_L[ lower.tri( START_L, diag = TRUE ) ] )
    }

    #- output:
    return( list( starts = starts, us = us ) )
}