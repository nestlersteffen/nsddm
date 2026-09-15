
#' Fit a hierarchical four-parameter DDM to data
#'
#' @param rt A vector with reaction times
#' @param xs A vector with choices
#' @param id A vector to identify which times and choices belongs to which level 2-unit
#' @param estimator A string indicating use of "ML" or "Bayes", defaults to "ML"
#' @param control A list of control arguments
#' @param debug A control argument for debugging, for internal use only
#' @export

hddm <- function( rt=NULL, xs=NULL, estimator="ML", control=ddm_args(), debug=FALSE )
{

	#--- step 0: make args_list:
	missCtrl <- missing( control )
	if ( !missCtrl && !inherits( control, "ddm_args" ) ) {
  		if(!is.list(control)) { stop("'control' has to be a list.") }
  		args_list <- do.call( ddm_args, control )
 	} else {
 		args_list <- control
 	}
 	
	#--- to be added
  
	#- finally...
	result <- "Hi"
	# class( result ) <- "nsddm"
  	return( result )

}