
ddm_agh_gradfct_singleperson <- function( parm_table=NULL,  
    rti=NULL, xsi=NULL, pts=NULL, wgh=NULL, MU=NULL, L=NULL, SIGMA=NULL, 
    invSIGMA=NULL, args=NULL, lan_weights=NULL, ddm_fun=NULL ) 
{
    #- matrix to save the results:
    npts    <- nrow( pts )
    lgt_tab <- nrow( parm_table )
    out     <- matrix( 0, nrow = npts, ncol = lgt_tab + 1 )    
      
    for ( pp in seq( npts ) ) {

        #- get point:
        Ppp <- pts[pp,] 
        
        #- compute loglik of response data:
        fapp <- -1*ddm_fun( us = pts[pp,], rt=rti, xs=xsi, 
            MU=MU, SIGMA=SIGMA, args=args, weights=lan_weights )

        #- weight the resulting value 
        fi <- exp( fapp ) * wgh[ pp ]
        out[ pp, 1 ] <- fi

        #- some pre computations:
        alphainvSIGMA <- t( Ppp - MU )%*%( invSIGMA )
        
        #- now the gradient for the point:
        for ( jj in 1:lgt_tab ) {
            
            #- get information for derivative
            free_jj <- parm_table[jj,]
            type_jj <- free_jj$type
            pos <- c( free_jj$pos1, free_jj$pos2 ) 
            
            #- compute derivative  
            if ( type_jj == "L" | type_jj == "SIGMA" ) {
                tmp_out <- ddm_derivatives( MU=MU, L=L, invSIGMA=invSIGMA, pos=pos, 
                    type_jj=type_jj, alphainvSIGMA=alphainvSIGMA ) 
            } else {
                tmp_out <- alphainvSIGMA[ pos[1] ]
            }                                   
            out[ pp, jj + 1 ] <- tmp_out*fi
        } 
    }

    #- final computations:
    tmp <- colSums( out )
    gri <- tmp[-1]
    lli <- tmp[ 1] 
    fin <- list( gri = gri, lli = as.numeric( lli ) )
    return( fin )
}

ddm_agh_gradfct <- function( parm=NULL, data_list=NULL, parm_list=NULL, 
    parm_table=NULL, pts_array=NULL, wgh=NULL, args=NULL,
    lan_weights=NULL, ddm="four", both=FALSE )
{
  
    #- insert parm:
    parm_list <- ddm_include_free_parameters( parm=parm, parm_list=parm_list, 
        parm_table=parm_table )
  
    #- get parameter matrices:
    parm_list <- ddm_compute_sigma( parm_list=parm_list, type_sigma=args$type_sigma )
    MU        <- parm_list[["MU"]]
    L         <- parm_list[["L"]]
    SIGMA     <- make_pd(parm_list$SIGMA, tol_factor = 1e-6)
    invSIGMA  <- base::solve(SIGMA) 

    #- get data:
    rts  <- data_list[["rts"]]
    xs   <- data_list[["xs"]]
    I    <- data_list[["I"]]
    info <- data_list[["info"]]

    #- select function depending on ddm:
    ddm_fun <- switch( ddm,
        "four"  = ddm4_random_nllfct,
        "seven" = ddm7_random_nllfct,
        stop( paste( "The 4- or 7-parameter DDM can be estimated." ) )
    )

    #- get path:
    lan_path <- if ( args$use_lan && args$use_rcpp ) {
        paste0( system.file( paste0("extdata/", ddm ), package = "nbddm" ), .Platform$file.sep )
    } else {""}
  
    #- make vector to save ll-values
    gi <- zi <- vector("list", I )
  
    #- compute loglik for each person
    for ( i in seq( I ) ) {

        # print( i )
        
        #- collect person data:
        idx1 <- info[i,3] - info[i,2] + 1 
        idx2 <- info[i,3]
        rti  <- rts[idx1:idx2]
        xsi  <- xs[idx1:idx2]
        
        #- compute gradient parts:
        if ( args$use_rcpp ) {

            int <- ddm_agh_gradfct_singleperson_rcpp( 
                typevec=as.vector(parm_table$ntype), 
                posmat=as.matrix(parm_table[,c("pos1","pos2")]), 
                rti=rti, xsi=xsi, pts=pts_array[,,i], wgh=wgh[,i], 
                MU=MU, L=L, SIGMA=SIGMA, invSIGMA=invSIGMA, 
                type_ddm=args$type_ddm, type_alpha=args$type_alpha, 
                ddm=ddm, kmax=args$kmax, delta=args$delta, 
                use_lan=args$use_lan, lan_path=lan_path )

        } else {

            int <- ddm_agh_gradfct_singleperson(parm_table=parm_table, 
                rti=rti, xsi=xsi, pts=pts_array[,,i], wgh=wgh[,i], 
                MU=MU, L=L, SIGMA=SIGMA, invSIGMA=invSIGMA, 
                args=args, lan_weights=lan_weights, ddm_fun=ddm_fun )

        }
        
        #- save in the correct list
        gi[[ i ]] <- int[["lli"]]
        zi[[ i ]] <- int[["gri"]] 
    
    }
    
    #- each person-specific integral vector is divided by gi:
    grad <- lapply( 1:length( gi ), function(tt) { zi[[tt]]*1/gi[[tt]] })  
  
    #- compute gradient
    grad <- Reduce("+", grad )
    #- return gradient
    if ( !both ) { res <- -1*grad
    } else { 
        #- log-lik value:
        llfct <- lapply( gi, function(i) log(i) )
        llfct <- Reduce("+", llfct )
        res <- list( "objective"=-1*llfct, "gradient"=-1*as.vector( grad ) ) 
    }
    return( res )
}

