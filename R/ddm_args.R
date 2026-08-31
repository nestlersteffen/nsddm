
ddm_args <- function( 
    # general args:  
    type_ddm = "std", # alternative: navfuss
    type_alpha = "dao", # alternative: logit
    type_sigma = "cholesky",
    use_tmb_modes = TRUE,
    use_rcpp = TRUE, 
    use_lan = FALSE,
    kmax = 5000L,
    delta = 1e-29,
    # max-lik specific agruments:
    nPoints = 3,
    method = "AGH",
    qmcType = "Sobol",
    isSeed = 123,
    maxit_out = 50,
    maxit_inn = 500,
    x_tol_out = 10e-4,
    fx_tol_out = 10e-6,
    verbose_out = TRUE, 
    verbose_inn = FALSE,
    lambda = 0.5,
    lambda_min = 0.5,
    n_threads = 1L,
    # bayes specific arguments:
    biter = 5000, 
    burnin = 2500, 
    nchain = 1, 
    tau2 = 0.1,
    tau2_A = c(0.50, 0.05, 0.05, 0.01)^2,
    tau2_B = c(0.10, 0.10, 0.05)^2, 
    use_adapt_dao = FALSE,
    pi_mix = 0.85, 
    pi_mix1= 0.65,
    pi_mix2= 0.30,
    pi_mix3= 0.05,
    epsilon = 0.10, 
    R = 20,
    type_proposal = "pmwg", # pmwg, mixture, standard
    type_inits = "none", # none, ml
    m = NULL,
    M = NULL,
    nu0 = NULL,
    S0 = NULL,
    # another list fpr nlminb specific args:
    control_nlminb = NULL ) 
{
    #- some checks:
    if ( !( type_ddm %in% c( "std", "navfuss", "dao" ) ) ) {
        stop("Method to compute DDM density not available.") }

    if ( !( type_alpha %in% c( "dao", "logit" ) ) ) {
        stop("Method to transform alpha not available.") }

    if ( !( method %in% c( "AGH", "QMC", "IS" ) ) ) {
        stop("Method to approximate integrals for ML not available.") }

    if (!is.numeric(n_threads) || n_threads < 0 || n_threads != as.integer(n_threads)) {
        stop("n_threads must be positive and a number.")
    }

    #- we make a list with all args related to ML:
    maxlik_list <- list( 
        nPoints = nPoints,
        method = method,
        qmcType = qmcType,
        isSeed = isSeed,
        lambda = lambda,
        lambda_min = lambda_min,
        maxit_out = maxit_out,
        x_tol_out = x_tol_out,
        fx_tol_out = fx_tol_out,
        verbose_out = verbose_out,
        n_threads = n_threads
    )

    #- we make a list with all args related to Bayes:
    bayes_list <- list( 
        biter = biter, 
        burnin = burnin, 
        nchain = nchain, 
        tau2 = tau2, 
        tau2_A = tau2_A,
        tau2_B = tau2_B,
        use_adapt_dao = use_adapt_dao,
        pi_mix = pi_mix, 
        pi_mix1 = pi_mix1,
        pi_mix2 = pi_mix2,
        pi_mix3 = pi_mix3,
        epsilon = epsilon, 
        R = R,
        type_proposal = type_proposal,
        type_inits = type_inits,
        m = m,
        M = M,
        nu0 = nu0,
        S0 = S0
    )

    #- we make a list for optimization in nlminb:
    default_list <- list( 
        iter.max = maxit_inn,
        eval.max = 1000,
        abs.tol = (.Machine$double.eps*10), 
        rel.tol = 1e-6,
        step.min = 1,#2.2e-10, 
        x.tol = 1.5e-8,
        trace = 0 
    )
    #- check verbose_inn:
    if ( verbose_inn ) { default_list$trace = 1 }
    #- add to final list:
    nlminb_list <- c( control_nlminb, 
        default_list[ !( names( default_list ) %in% names( control_nlminb ) ) ] ) 

    #- make args:
	args <- list( 
        type_ddm = type_ddm, 
        type_alpha = type_alpha, 
        type_sigma = type_sigma,
        use_tmb_modes = use_tmb_modes,
        use_rcpp = use_rcpp, 
        use_lan = use_lan,
        kmax = kmax,
        delta = delta,
        maxlik_list = maxlik_list,
        bayes_list = bayes_list,
        nlminb_list = nlminb_list )
	return( args )
}