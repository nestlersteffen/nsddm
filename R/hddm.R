
#' Fit a hierarchical four-parameter DDM to data
#'
#' @param rt A vector with reaction times
#' @param xs A vector with choices
#' @param id A vector to identify which times and choices belongs to which level 2-unit
#' @param estimator A string indicating use of "ML" or "Bayes", defaults to "ML"
#' @param control A list of control arguments
#' @param verbose A control argument to see the optimization history, defaults to "TRUE"
#' @export

hddm <- function( rt=NULL, xs=NULL, id=NULL, estimator="ML", control=ddm_args(), verbose=TRUE )
{

	#--- step 0: make args_list:
	missCtrl <- missing( control )
	if ( !missCtrl && !inherits( control, "ddm_args" ) ) {
  		if(!is.list(control)) { stop("'control' has to be a list.") }
  		args <- do.call( ddm_args, control )
 	} else {
 		args <- control
 	}

 	#--- step 1: some simple checks
 	if ( length( rt ) != length( xs ) ) {
 		stop("Responses rt and choices x have different length.")
 	}

 	#--- step 2: make a parm_table (at present this cannot be changed)
	parm_list <- list( MU = matrix(0, nrow = 4, ncol = 1 ), L = diag( 0, 4 ), SIGMA = diag( 0, 4) )
	parm_table <- data.frame(
	  type=c("MU","MU","MU","MU","L","L","L","L","L","L","L","L","L","L"),
	  pos1=c(1,2,3,4,1,2,3,4,2,3,4,3,4,4),
	  pos2=c(1,1,1,1,1,1,1,1,2,2,2,3,3,4) )
	parm_table$index <- seq(1,nrow(parm_table),1)
	parm_table$type  <- as.character(parm_table$type)
	parm_table$ntype <- c( rep(0,4), rep(1,10) )

	#--- step 3: make data_list
    info <- as.data.frame( table( id ) )
    info[,1] <- as.numeric( levels( info[,1]) )[info[,1]]
    info <- cbind( info, base::cumsum( info[,2] ) )
    info <- as.matrix( info )
    attr( info, "dimnames" ) <- NULL
    data_list <- list( rts=rt, xs=xs, info=info, I=nrow( info ), K=4 )

 	#--- step 4: fit the model:
 	if ( estimator == "ML" ) {

		#--- obtain some start values:
		start_vals <- ddm4_getstarts( data_list=data_list, args=args, ddm="four" )
		parm_table$starts <- round( start_vals$starts, 4 )

		#--- set verbose:
		if ( verbose ) args$maxlik_list$verbose_out <- TRUE

		#--- get estimates:
		result <- ddm_fit_ml( parm_table=parm_table, parm_list=parm_list, 
			data_list=data_list, args=args, ddm="four" )

		if ( !result$converged ) warning( "Algorithm has not converged.") 

 	} else if ( estimator == "Bayes" ) {

 		fit_bayes <- suppressWarnings( tryCatch( { 
 			ddm_random_bayes( data_list=data_list, parm_table=parm_table, args=args, 
 				ddm="four",verbose=verbose ) },
                error = function(e) { NULL } ) )
 		if ( !is.null( fit_bayes ) ) result <- fit_bayes else stop("Something went wrong in Bayesian estimation.")

 	}

	#- finally...
	class( result ) <- "nsddm"
  	return( result )

}