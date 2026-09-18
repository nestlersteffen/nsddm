
#---- this function contains routines to fit the different GLMMs:

ddm_fit_ml <- function( parm_table=NULL, parm_list=NULL, data_list=NULL,
	args=NULL, ddm=NULL )
{

	#- get the start values
	starts <- parm_table$starts
	
	#- the algorithm has not converged...so far
	converged    <- FALSE
	warning_vcov <- TRUE

	#- start making points: 
	pts_list  <- ddm_makepoints( dimension = ncol( parm_list[["L"]] ), 
		nPoints = args$maxlik_list$nPoints, args=args ) 

	#- select modes-function depending on ddm:
    ddm_modes_fun <- switch( ddm,
        "four"  = ddm4_getmodes,
        "seven" = ddm7_getmodes,
        stop( paste( "The 4- or 7-parameter DDM can be estimated." ) )
    )

	#- some definitions for the while-loop:
	maxit_out  	   <- args$maxlik_list$maxit_out
	x_tol_out      <- args$maxlik_list$x_tol_out
	fx_tol_out     <- args$maxlik_list$fx_tol_out
	verbose_out    <- args$maxlik_list$verbose_out
	lambda_min     <- args$maxlik_list$lambda_min
	lambda         <- args$maxlik_list$lambda
	warning_raneff <- 0

	#- let's start the while loop:
	parm_new <- starts
	obj_new  <- 0
	diterate <- TRUE
	diter    <- 0

	while ( diterate ) {
 				
		#- old parm_vector:	 
  		parm_old <- parm_new  
   		obj_old  <- obj_new
   		diter    <- diter + 1 
   		lambda   <- max( lambda_min, lambda * (1 - diter/maxit_out) ) 
    			
   		#- Step 1: get random effect estimates given current parm:
		modes_list <- suppressWarnings( tryCatch( { 
	      	ddm_modes_fun( parm = parm_old, data_list = data_list, 
	      		parm_list = parm_list, parm_table = parm_table,
	      		args=args,check_t0=TRUE ) },
	      	error = function(e) { NULL },
	      	warning = function(w) { NULL } ) )
			    
	    # error handling:
	    if ( is.null( modes_list ) ) {
			parm_new 	   <- parm_old
			warning_raneff <- 1
			break
		} else {
			if ( args$nlminb_list$trace == 1 ) print( "Modes ok!")
		}

		#- Step 2: make pts_array:
		adapt_pts <- ddm_adaptpoints( modes_list=modes_list, pts_list=pts_list, args=args )

		#- Step 2: make the optimization object:
		obj <- ddm_MakeOptFun( ddm_agh_gradfct_parallel, data_list=data_list, parm_list=parm_list, 
		  	parm_table=parm_table, pts_array=adapt_pts$pts, wgh=adapt_pts$wgh, 
		  	args=args, ddm=ddm, both=TRUE ) 

		fit <- suppressWarnings( tryCatch( 
		    stats::nlminb( start=parm_old, 
		        objective=obj$fn, gradient=obj$gr, 
		        control=args$nlminb_list ),
		    error   = function(e) { NULL },
		    warning = function(w) { NULL }
		) )

		# fit <- suppressWarnings( tryCatch( 
		#     stats::optim( par=parm_old, fn=obj$fn, gr=obj$gr, method="BFGS",
		#         control=list(trace = 1) ),
		#     error   = function(e) { NULL },
		#     warning = function(w) { NULL }
		# ) )

		# error handling:
		if ( is.null(fit) || inherits(fit, "try-error") ) {
		    parm_new <- parm_old
		    break
		} else {
		    #parm_new <- fit$par
		    parm_new <- (1 - lambda) * fit$par + lambda * parm_old
		    obj_new  <- fit$objective
		}

	    # Step 3: check convergence
	    delta    <- crossprod( parm_old - fit$par )/crossprod( parm_old ) 
	    delta_ll <- abs( obj_new - obj_old ) / abs( obj_old )
		if ( delta < x_tol_out | delta_ll < fx_tol_out | diter > maxit_out ) {
   			converged <- TRUE
   			diterate  <- FALSE
   			# finaler Wert ist fit$par, nicht parm_new:
    		parm_new  <- fit$par
   		} 

   		#- make output
		if ( verbose_out ) {
		   	xx <- as.character()
		   	tb <- c( round( parm_new, 4 ) )
		   	for (j in 1:length(tb)) { xx <- paste0(xx, sprintf(" %g", tb[j])) }
		  	h1 <- paste0("Parms: ", xx )
			cat(h1, "\n")
			utils::flush.console()
		}
			    
	}

  	#- --------------------------------------------------
	if ( converged ) {
	
		#- compute final ll:
		ll <- -1*fit$objective
			
		#- add parm_new to parm_table:
		parm_table$est <- parm_new

		#- rescale to sigma:
        parm_table <- ddm_change_to_sigma( parm_list=parm_list, parm_table=parm_table, args=args )	

		#- compute hessian:
		tmp_args <- args 
		tmp_args$type_sigma <- "sigma"
		hessian <- tryCatch( 
			nloptr::nl.jacobian( x0=parm_table$est, fn=ddm_agh_gradfct_parallel, 
				data_list=data_list, parm_list=parm_list, parm_table=parm_table, 
				pts_array=adapt_pts$pts, wgh=adapt_pts$wgh, args=tmp_args, ddm=ddm,
				both=FALSE ),
			error = function(e) {
                NULL
    		}
    	)

    	if ( is.null( hessian ) ) {
    		
    		parm_table$est <- parm_table$se  <- NA
			parm_table$z   <- parm_table$p   <- NA
			parm_table$lCI <- parm_table$uCI <- NA
			ll <- dev <- aic <- vcov <- NULL
			converged_nlminb <- fit$convergence
	    
	    } else {

	    	hessian <- ( hessian + t( hessian ) ) / 2
		
			#- compute vcov ( with correcton in case problems occur )
			vcov <- tryCatch(
			    base::solve(hessian),
			    error = function(e) {
			        eig        <- eigen(hessian, symmetric = TRUE)
			        eig$values <- pmax(eig$values, 1e-6)
			        hessian_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
			        chol2inv(chol(hessian_pd))
			    }
			)

			#- get standard errors:
			ses <- suppressWarnings( sqrt( diag( vcov ) ) )
			ses[ is.nan(ses) | ses <= 0 ] <- NA
			warning_vcov <- any( is.na(ses) )

			#- save results:
			parm_table$se  <- ses
			parm_table$z   <- with( parm_table, est/se )   
			parm_table$p   <- 2 * pnorm( -abs( parm_table$z ) )
			parm_table$lCI <- parm_table$est - 1.96*parm_table$se
			parm_table$uCI <- parm_table$est + 1.96*parm_table$se

			#- deviance and aic
			dev <- -2*ll
			aic <- dev + 2*length( parm_new )
			converged_nlminb <- fit$convergence

			#-insert parm_new to parm_list:
			parm_list <- ddm_include_free_parameters( parm=parm_table$est, parm_list=parm_list, 
	        	parm_table=parm_table )	

	    }

	} else {
		parm_table$est <- parm_table$se  <- NA
		parm_table$z   <- parm_table$p   <- NA
		parm_table$lCI <- parm_table$uCI <- NA
		ll <- dev <- aic <- vcov <- NULL
		converged_nlminb <- 1
	}

	#- save and return: 
  	res <- list( parm = parm_new, parm_table = parm_table, parm_list = parm_list,
	    data_list = data_list, args = args, deviance = dev, ll = ll, aic = aic, 
	    vcov = vcov, converged = converged, no_iterations = diter, 
	    warning_vcov = warning_vcov, warning_raneff = warning_raneff, 
	    converged_nlminb = converged_nlminb )            
  	return( res )

}