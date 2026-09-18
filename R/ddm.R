
#' Fit a four-parameter DDM to the data of a single person
#'
#' @param rt A vector with reaction times
#' @param xs A vector with choices
#' @param estimator A string indicating use of "ML" or "Bayes", defaults to "ML"
#' @param control A list of control arguments
#' @param verbose A boolean control argument to see the optimization history, defaults to TRUE
#' @export

ddm <- function( rt=NULL, xs=NULL, estimator="ML", control=ddm_args(), verbose=TRUE )
{

	#--- step 0: make args_list:
	missCtrl <- missing( control )
	if ( !missCtrl && !inherits( control, "ddm_args" ) ) {
  		if(!is.list(control)) { stop("'control' has to be a list.") }
  		args <- do.call( ddm_args, control )
 	} else {
 		args <- control
 	}

 	#--- step 1: some simple checks:
 	if ( length( rt ) != length( xs ) ) {
 		stop("Responses rt and choices x have different length.")
 	}

 	#--- step 2: fit the model
 	if ( estimator == "ML" ) {

 		#--- make some start values
 		starts <- ddm_inits_mu( rts=rt, type_alpha=args$type_alpha, ddm="four" ) 

 		if ( !args$use_lan ) {

	 		#--- make TMB data:
	 		tmb_alpha <- if ( args$type_alpha == "dao" ) 0L else 1L
	 		tmb_ddm   <- if ( args$type_ddm == "navfuss" ) 0L else 1L
	 		tdata     <- list( rts=rt, xs=xs, method=tmb_ddm, type_alpha=tmb_alpha, model="ddm4_tmb")
	        
	        #--- make TMB function
	   		obj <- TMB::MakeADFun( data=tdata, parameters=list( parm=starts ), 
	   			DLL="nsddm_TMBExports", silent=TRUE )

	   		#--- fit data:
	   		fit_ml <- suppressWarnings( tryCatch( { 
	   			nlminb( obj$par, objective=obj$fn, gradient=obj$gr, control=list(trace=verbose) ) },
        			error = function(e) { NULL } ) )

        } else {
    
    		#--- make a closure used for optimization:    
	       	lan_optfct <- function( parm, rt=NULL, xs=NULL, both=TRUE ) {
	       		out <- ddm4_lanll_grad_weights_rcpp( parm, rt, xs )
	        	if ( both ) return( list( "objective"=out$objective, "gradient"=as.numeric(out$gradient ) ) )
	        	else return ( as.numeric(out$gradient ) )
	       	}

    		#--- load the weight matrices:
    		ddm_load_weights_rcpp( ddm="four" )

    		#--- make the function:
    		obj <- ddm_MakeOptFun( lan_optfct, rt=rt, xs=xs, both=TRUE )
   
   			#--- fit data:
    		fit_ml <- suppressWarnings( tryCatch( { 
	   			nlminb( starts, objective=obj$fn, gradient=obj$gr, control=list(trace=verbose) ) },
        			error = function(e) { NULL } ) )
    	}

    	opt_ok <- !is.null( fit_ml ) && fit_ml$objective != 1e+06 && fit_ml$convergence == 0

    	#--- error handling and standard errors:

    	if ( opt_ok ) {

        	parm <- fit_ml$par
	    	if ( !args$use_lan ) {
	    		HESS  <- tryCatch( obj$he( ), error = function(e) NULL )
        	} else {
        		HESS <- tryCatch( {
            		numDeriv::jacobian( x=parm, func=lan_optfct, rt=rt, xs=xs, both=FALSE ) },
                	error = function(e) NULL
        		)
        	}

        	#--- a check
			if ( is.null( HESS ) || any( is.na( HESS ) ) || any( is.infinite( HESS ) ) ) {
        		stop( "Hessian is not well defined.")
        	}

        	#--- get covariance matrix and standard errors
        	vcov  <- base::solve( HESS )
        	ses   <- suppressWarnings( sqrt( diag( vcov ) ) )

	 	} else {
	 		stop( "Algorithm has not converged." )
		}

 		#--- make a result table:
 		tab <- data.frame( Par=paste0("alpha",c(1:4)), Est=parm, SE=ses )
 		tab$z   <- with( tab, Est/SE )   
		tab$p   <- 2 * pnorm( -abs( tab$z ) )
		tab$lCI <- tab$Est - 1.96*tab$SE
		tab$uCI <- tab$Est + 1.96*tab$SE

		#- output:
    	result <- list( parm_table=tab, vcov=vcov )

 	} else if ( estimator == "Bayes" ) {

 		fit_bayes <- suppressWarnings( tryCatch( { 
        	ddm_bayes_singleperson( rt=rt, xs=xs, args=args, ddm="four", verbose=verbose ) },
                error = function(e) { NULL } ) )
 		if ( !is.null( fit_bayes ) ) result <- fit_bayes else stop("Something went wrong in Bayesian estimation.")
    	
 	}
 	
	#--- finally...
	class( result ) <- "nsddm"
  	return( result )

}