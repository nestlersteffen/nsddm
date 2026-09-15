
#--- function that loads the weight matrices of the neural network:

ddm_load_weights <- function( ddm="four" ) 
{
    
    #- get path:
    path <- system.file( paste0("extdata/", ddm ), package = "nsddm")
    if (path == "") stop("extdata directory not found in nsddm package")
    path <- paste0( path, .Platform$file.sep )

    #- define dimensions of the matrices:
    if ( ddm=="four" ) {
        dims_W <- list(c(6,100), c(100,100), c(100,120), c(120,1))
    } else {
        dims_W <- list(c(9,100), c(100,100), c(100,120), c(120,1))
    }
    dims_b <- list(100, 100, 120, 1)
    
    #- get the matrices:
    weights <- list()
    for (i in 1:4) {
        W <- readBin(paste0(path, "W", i, ".bin"), what="double", 
                     n=prod(dims_W[[i]]))
        W <- matrix(W, nrow=dims_W[[i]][1], ncol=dims_W[[i]][2])
        
        b <- readBin(paste0(path, "b", i, ".bin"), what="double",
                     n=dims_b[[i]])
        
        weights[[i]] <- list(W=W, b=b)
    }
    
    #- return them
    return(weights)
}

ddm_load_dnn <- function( ddm="four" ) 
{
    #- get path:
    path <- system.file( paste0("extdata/", ddm ), package = "nsddm")
    if (path == "") stop("extdata directory not found in nsddm package")
    path <- paste0( path, .Platform$file.sep )

    #- load model:
    model <- keras3::load_model( paste0( path, "Trained_Model_",ddm,"_param.keras" ) )

    #- generate some fake data:
    if ( ddm=="four" ) { 
        df <- ddm4_simulate( 1, c(0,log(0.25),log(0.25),log(0.1) ), "dao" )
        invisible( ddm4_lan_llfct_dnn( model, df$rt, df$xs, 0, 0.5, 0.25, 0.1 ) )
    } else {
        df <- ddm7_simulate( 1, c(0,log(0.25),log(0.25),log(0.1), log(1), log(0.3), log(0.1) ), "dao" )
        invisible( ddm7_lan_llfct_dnn( model, df$rt, df$xs, 0, 0.5, 0.25, 0.1, 1, 0.3, 0.1 ) )
    }
    
    #- return them
    return( model )
}

ddm_load_weights_rcpp <- function( ddm="four" ) {
    path <- system.file( paste0("extdata/", ddm ), package = "nsddm")
    if (path == "") stop("extdata directory not found in nsddm package")
    path <- paste0( path, .Platform$file.sep )   # trailing slash
    lan_load_weights_export( path, ddm )
}

#- function for the forward Pass:
lan_forward <- function(input, weights) 
{
    x <- input
    #- layer 1-3: with tanh
    for (i in 1:3) {
        x <- tanh( t(weights[[i]]$W) %*% x + weights[[i]]$b) # falsch: weights[[i]]$W %*% x
    }
    #- output layer: linear
    x <- t(weights[[4]]$W) %*% x + weights[[4]]$b # falsch: weights[[4]]$W %*% x
    return( x )
}

lan_forward_and_backward <- function( input, weights, idx=7 ) 
{
    #- forward pass:
    activations <- list()
    activations[[1]] <- input
    x <- input
    for (i in 1:3) {
        x <- tanh( t(weights[[i]]$W) %*% x + weights[[i]]$b )
        activations[[i + 1]] <- x
    } 
    x <- t(weights[[4]]$W) %*% x + weights[[4]]$b   
    #- backward pass:
    d <- matrix(1.0) 
    for (i in 4:2) {
        d <- (weights[[i]]$W %*% d) * (1 - activations[[i]]^2)
    }
    grad_input <- weights[[1]]$W %*% d   # Gradient für alle 6 Inputs
    return( c( x, grad_input[1:idx] ) )
}

#- wrappers for lan log-lik that uses the weight matrices:
ddm4_lanll_weights <- function( weights, a, v, t0, z, rt, xs, K ) 
{
    tmpMat <- matrix( c(a, v, t0, z), nrow=length(rt), ncol=K, byrow=TRUE)
    tmpMat <- cbind( tmpMat, xs, rt)
    ll <- apply( tmpMat, 1, function(row) lan_forward(row, weights))
    return( sum( as.vector(ll) ) )
}

ddm7_lanll_weights <- function( weights, a, v, t0, z, sv, sz, st0, rt, xs, K ) 
{
    tmpMat <- matrix( c(a, v, t0, z, sv, sz, st0), nrow=length(rt), ncol=K, byrow=TRUE)
    tmpMat <- cbind( tmpMat, xs, rt)
    ll <- apply( tmpMat, 1, function(row) lan_forward(row, weights))
    return( sum( as.vector(ll) ) )
}

#- a wrapper for the neural network with gradient: 

ddm4_lanll_grad_weights <- function( us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL ) 
{
    #- Step 1: transform parameters
    v   <- us[1]
    z   <- exp(us[3])
    t0  <- exp(us[4])
    a   <- exp(us[2]) + z
    
    #- Step 2: make matrix
    tmpMat <- matrix( c(a, v, t0, z), nrow=length(rt), ncol=4, byrow=TRUE)
    tmpMat <- cbind( tmpMat, xs, rt)
    
    #- Step 3: compute approximate ll and gradient:
    tmpRes <- apply( tmpMat, 1, function(row) lan_forward_and_backward( row, weights, 4 ) )
     
    #- Step 4: compute gradient with Jacobi matrix:
    g <- rowSums(tmpRes[-1, ])  # g[1]=a, g[2]=v, g[3]=t0, g[4]=z
    J <- matrix(0, nrow=4, ncol=4)
    J[1,2] <- exp(us[2]); J[1,3] <- z; J[2,1] <- 1; J[3,4] <- t0; J[4,3] <- z 
    grad_us <- t(J) %*% g
    return( list( "objective"=-1*sum(tmpRes[1,]), "gradient"=-1*as.numeric( grad_us) ) ) 
}

ddm7_lanll_grad_weights <- function( us=NULL, rt=NULL, xs=NULL, args=NULL, weights=NULL ) 
{
    
    #- Step 1: transform parameters
    v   <- us[1]
    sv  <- exp(us[5])
    sz  <- exp(us[6])
    st0 <- exp(us[7])
    z   <- exp(us[3]) + 0.5*sz
    t0  <- exp(us[4]) + 0.5*st0
    a   <- exp(us[2]) + z + 0.5*sz
    
    #- Step 2: make matrix
    tmpMat <- matrix( c(a, v, t0, z, sv, sz, st0), nrow=length(rt), ncol=7, byrow=TRUE)
    tmpMat <- cbind( tmpMat, xs, rt)
    
    #- Step 3: compute approximate ll and gradient:
    tmpRes <- apply( tmpMat, 1, function(row) lan_forward_and_backward( row, weights, 7 ) )
     
    #- Step 4: compute gradient with Jacobi matrix:
    g <- rowSums(tmpRes[-1, ])  # g[1]=a, g[2]=v, g[3]=t0, g[4]=z, g[5]=sv, g[6]=sz, g[7]=st0
    J <- matrix(0, nrow=7, ncol=7)
    J[1,2] <- exp(us[2]); J[1,3] <- exp(us[3]); J[1,6] <- sz 
    J[2,1] <- 1; J[3,4] <- exp(us[4]); J[3,7] <- 0.5 * st0 
    J[4,3] <- exp(us[3]); J[4,6] <- 0.5 * sz 
    J[5,5] <- sv; J[6,6] <- sz; J[7,7] <- st0
    grad_us <- t(J) %*% g
    return( list( "objective"=-1*sum(tmpRes[1,]), "gradient"=-1*as.numeric( grad_us) ) ) 
}

#- wrappers to use lan_gradient - functions in getmodes:

ddm_lan_nllfct_gradient_wrap <- function(alpha=NULL, rt=NULL, xs=NULL, MU=NULL, SIGMA=NULL, 
    invSIGMA=NULL, ddm=NULL, both=TRUE) 
{
    out <- ddm_lan_nllfct_gradient_rcpp( alpha, rt, xs, MU, SIGMA, invSIGMA, ddm )
    if (both) return( list( "objective"=out$objective, "gradient"=as.numeric( out$gradient ) ) )
    else return ( as.numeric(out$gradient ) )
}

#- wrapper for lan-loglik-fct that uses the dnn object (not used anymore):

ddm4_lan_llfct_dnn <- function( dnn=NULL, rt=NULL, xs=NULL, v=NULL, a=NULL, z=NULL, t0=NULL )
{
    #- no. of reaction times:
    n <- length( rt)
    #- make temporary data matrix:
    tmpMat <- matrix(0, nrow = n, ncol = 6)
    tmpMat[,1] <- a; tmpMat[,2] <- v; tmpMat[, 3] <- t0; tmpMat[, 4] <- z
    tmpMat[,5] <- xs
    tmpMat[,6] <- rt
    #- compute ll values:
    ll <- dnn( tmpMat )
    return( as.array( sum( ll ) ) )
}

ddm4_lan_llfct_gradient_dnn <- function( dnn=NULL, rt=NULL, xs=NULL, v=NULL, a=NULL, z=NULL, t0=NULL )
{
    #- no. of reaction times:
    n <- length( rt )
    
    #- parameter to "watch" for:
    a0  <- tf$Variable(a,  dtype = tf$float32)
    v0  <- tf$Variable(v,  dtype = tf$float32)
    z0  <- tf$Variable(z,  dtype = tf$float32)
    t00 <- tf$Variable(t0, dtype = tf$float32)
    
    #- data are constants:
    rt_tensor <- tf$constant( rt, dtype = tf$float32)
    xs_tensor <- tf$constant( xs, dtype = tf$float32)
    ones_n    <- tf$ones( shape = list(n), dtype = tf$float32 )
    
    with( tf$GradientTape() %as% tape, {
        
        # broadcast parms to length n:
        col_a  <- a0  * ones_n
        col_v  <- v0  * ones_n
        col_t0 <- t00 * ones_n
        col_z  <- z0  * ones_n
        # make a temporary data matrix:
        tmpMat <- tf$stack( list( col_a, col_v, col_t0, col_z, xs_tensor, rt_tensor), axis = 1L)
        # compute ll value
        ll <- dnn( tmpMat )
        ll_sum <- tf$reduce_sum(ll)
    } )
    
    gr <- tape$gradient( ll_sum, list(a0, v0, z0, t00) )
    gr <- lapply( gr, function(x) as.array( x ) )
    return( list( "objective"=-1*as.array( ll_sum ), "gradient"=-1*do.call("c", gr ) ) )
}