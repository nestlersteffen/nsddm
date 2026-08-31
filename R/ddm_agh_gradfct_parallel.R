
#------- MAIN

ddm_agh_gradfct_parallel <- function( parm=NULL, data_list=NULL, parm_list=NULL, 
    parm_table=NULL, pts_array=NULL, wgh=NULL, args=NULL, ddm="four", both=FALSE )
{
  
    #- insert parm:
    parm_list <- ddm_include_free_parameters( parm=parm, parm_list=parm_list, 
        parm_table=parm_table )
  
    #- get parameter matrices:
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

    #- get path:
    lan_path <- if ( args$use_lan && args$use_rcpp ) {
        paste0( system.file( paste0("extdata/", ddm ), package = "nbddm" ), .Platform$file.sep )
    } else {""}
  
    #- get gradient and llfct in parallel way:
    result <- ddm_agh_gradfct_allpersons_rcpp_save( 
        typevec=as.vector(parm_table$ntype), posmat=as.matrix(parm_table[,c("pos1","pos2")]), 
        rts=rts, xs=xs, info=info, I=I,
        pts_array=pts_array, wgh=wgh,  
        MU=MU, L=L, SIGMA=SIGMA, invSIGMA=invSIGMA, 
        type_ddm=args$type_ddm, type_alpha=args$type_alpha, 
        ddm=ddm, kmax=args$kmax, delta=args$delta, 
        use_lan=args$use_lan, lan_path=lan_path, 
        n_threads=args$maxlik_list$n_threads )

    #- each person-specific integral vector is divided by gi:
    grad <- lapply( c(1:I), function(i) { result$gri[i,]*(1/result$lli[i]) } )  
    grad <- Reduce("+", grad )
    
    #- return gradient
    if ( !both ) { res <- -1*grad
    } else { 
        #- log-lik value:
        llfct <- sum( log( result$lli ) )
        res <- list( "objective"=-1*llfct, "gradient"=-1*as.vector( grad ) ) 
    }
    return( res )
}

