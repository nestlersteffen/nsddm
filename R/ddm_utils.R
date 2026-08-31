
#--- function to include parms into parm_table 

ddm_include_free_parameters <- function( parm = NULL, 
  parm_list = NULL , parm_table = NULL )
{
  NP <- max(parm_table$index)
  NOP <- nrow(parm_table) 
  for (nn in 1:NOP){
    free_nn  <- parm_table[nn,]
    type     <- free_nn$type
    pos      <- as.numeric( free_nn[ c("pos1","pos2") ]  )
    index_nn <- free_nn$index
    x_nn <- parm_list[[ type ]]   
    x_nn[ pos[1] , pos[2] ] <- parm[ free_nn$index ]
    if ( type %in% c("SIGMA","R") ) {
      x_nn[ pos[2] , pos[1] ] <- parm[ free_nn$index ]
    }   
    parm_list[[ type ]] <- x_nn       
  }
  return(parm_list)
}

#--- function to make GH and AGH points:

ddm_makepoints <- function( dimension = NULL, nPoints = NULL, args = NULL )  
{
  if ( args$maxlik_list$method == "AGH" ) {
    #- make a grid
    idx <- as.matrix(expand.grid(rep(list(1:nPoints),dimension)))
    #- get the quadrature points
    tmp_val <- fastGHQuad::gaussHermiteData(nPoints)
    #- save the weights
    tmp_wgh <- matrix( tmp_val$w[idx], nrow(idx), dimension)
    #- save the points
    tmp_pts <- matrix( tmp_val$x[idx], nrow(idx), dimension)
    #- transform and exp the points for agh quadrature
    wgh <- apply( ( 1/sqrt( pi ) )*tmp_wgh, 1, prod )
    pts <- sqrt( 2 )*tmp_pts
    #- things for agh:
    tmp_exp <- as.matrix( exp( 0.5*rowSums( pts * pts ) ) )
    wghexp  <- apply( cbind( wgh, tmp_exp ), 1, prod )
  } else if ( args$maxlik_list$method == "QMC" ) {
    if ( args$maxlik_list$qmcType == "Halton" ) {
      set.seed( args$maxlik_list$isSeed )
      #pts <- randtoolbox::halton( n = nPoints, dim = dimension, scrambling = 1, seed = 42 )
      pts <- qrng::ghalton(n = nPoints, d = dimension, method = "generalized") 
    } else if ( args$maxlik_list$qmcType == "Sobol" ) {
      # pts <- randtoolbox::sobol( n = nPoints, dim = dimension, scrambling = 1, seed = 42 ) 
      pts <- qrng::sobol(n = nPoints, d = dimension, randomize = "Owen", seed = 42)
    } else {
      stop("Method to obtain QMC points not available.")
    }
    pts <- as.matrix( stats::qnorm( pts, mean = 0, sd = 1 ) )
    wgh <- wghexp <- rep(1,nPoints)
  } else if ( args$maxlik_list$method == "IS" ) {
    set.seed( args$maxlik_list$isSeed )
    pts <- mvtnorm::rmvnorm( nPoints, mean = rep(0, dimension), 
      sigma = diag( 1, dimension ) ) 
    wgh <- wghexp <- rep(1,nPoints)
  } else {
    stop("Method to obtain points not available.")
  }

  #- return results:
  return( list( pts = pts, wgh = wgh, wghexp = wghexp ) )
}

# adapt the points with the modes and the hessian of the random effects

ddm_adaptpoints_old2 <- function( parm_list = NULL, pts_list = NULL,
  ranef = NULL, args=NULL )  
{
  #- get ranef infos:
  MODE  <- ranef[["MODE"]]
  CHOL  <- ranef[["CHOL"]]
  iHESS <- ranef[["iHESS"]]
  #- transform points:
  pts <- pts_list$pts%*%t(CHOL)
  pts <- sweep( pts, 2, MODE, "+")
  den <- NULL
  wgh <- pts_list[["wghexp"]] 
  return( list( pts = pts, wgh = wgh ) )
}

ddm_adaptpoints_old <- function( modes_list = NULL, pts_list = NULL, args=NULL )  
{
  #- extract relevant information:
  I  <- length( modes_list )
  K  <- ncol( modes_list[[I]]$HESS )
  #- make pts_array:
  pts <- array( 0, dim = c( nrow( pts_list$pts ), K, I ) )
  #- fill the array:
  for ( i in seq(I) ) {
    #- get ranef infos:
    MODE <- modes_list[[i]]$MODE
    CHOL <- modes_list[[i]]$CHOL
    #- transform points:
    tmp <- pts_list$pts%*%CHOL
    pts[,,i] <- sweep( tmp, 2, MODE, "+")
  }
  wgh <- pts_list[["wghexp"]] 
  return( list( pts = pts, wgh = wgh ) )
}

ddm_adaptpoints <- function( modes_list = NULL, pts_list = NULL, args=NULL )  
{
  #- extract relevant information:
  I  <- length( modes_list )
  K  <- ncol( modes_list[[I]]$HESS )
  #- make pts_array:
  pts <- array( 0, dim = c( nrow( pts_list$pts ), K, I ) )
  wgh <- matrix( 0, nrow = nrow( pts_list$pts ), ncol = I )
  #- fill the array:
  for ( i in seq(I) ) {
    #- get ranef infos:
    MODE  <- modes_list[[i]]$MODE
    CHOL  <- modes_list[[i]]$CHOL
    iHESS <- modes_list[[i]]$iHESS
    #- transform points:
    tmp <- pts_list$pts%*%CHOL
    pts[,,i] <- sweep( tmp, 2, MODE, "+")
    if ( args$maxlik_list$method == "AGH" ) {
      wgh[,i] <- pts_list$wghexp
    } else {
      wgh[,i] <- 1/mvtnorm::dmvnorm( pts[,,i], MODE, iHESS )
    }
  }
  return( list( pts = pts, wgh = wgh ) )
}

#--- sun-exp-log trick:

ddm_sumexplog <- function( x = NULL ) 
{
  #- get maximimum:
  xmax <- max( x )
  #- compute sum:
  res  <- xmax + log( sum( exp( x - xmax ) ) )
  return( res ) 
}

#--- optimize function wrapper

ddm_MakeOptFun <- function( grad_fn, ... ) 
{
  #- the cache
  last_par <- NULL
  last_val <- NULL
    
  eval_cached <- function(par) {
    if ( !identical(par, last_par) ) {
          last_par <<- par
          last_val <<- grad_fn(par, ... )
      }
      last_val
  }
    
  list(
    fn = function(par) eval_cached(par)$objective,
    gr = function(par) eval_cached(par)$gradient
  )
}

#--- eigenvalue correction

# make_pd <- function(H, tol_factor = 1e-6) {
#     eig        <- eigen(H, symmetric = TRUE)
#     # tol relativ zum größten Eigenwert:
#     tol        <- tol_factor * max(abs(eig$values))
#     eig$values <- pmax(eig$values, tol)
#     H_pd       <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
#     H_pd       <- (H_pd + t(H_pd)) / 2
#     return(H_pd)
# }

make_pd <- function(H, tol_factor = 1e-3) {
    eig        <- eigen(H, symmetric = TRUE)
    eig$values <- pmax(eig$values, tol_factor )
    H_pd       <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    return(H_pd)
}