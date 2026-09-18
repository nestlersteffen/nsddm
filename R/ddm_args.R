
#' Control arguments for nsddm
#'
#' @description
#'
#' The values supplied in the function call replace the default values. A list is returned
#' that contains all final values of the arguments. The list is then used as
#' the \code{args} argument to functions in \code{\link{nsddm}}. The list was inspired
#' by the well-known \code{nlme} package.
#'
#' @param type_dmm The exact density that is used. Default is \code{std} the density reported in 
#'   Tuerlinckx (2004). The alternative is the density reported in Navarro & Fuss (2009): \code{navfuss}.
#' @param type_alpha The parametrization of the parameter vector. Default is \code{dao}, the parametrization
#'   used in Dao et al. (2025). The alternative, \code{logit}, is a relative starting point parametrization.
#' @param type_sigma How the covariance matrix is parametrized during hierarchical ML estimation. Default is 
#'   \code{cholesky}. 
#' @param use_tmb_modes Controls whether the modes in the hierarchical ML estimation routine should be estimated 
#' @param use_rcpp A logical indicating whether C++ or R be used during fitting. Default is \code{TRUE}.
#' @param use_lan A logical indicating whether the log-likelihood density values are obtained with a pre-trained 
#'   neural netwok. Default is \code{FALSE}.
#' @param use_tf A logical indicating whether the pre-trained neural netwok computations are done with the
#'   TensorFlow object. Default is \code{FALSE}. Can only be used when \code{use_rcpp=FALSE}.
#' @param kmax Maximum number of terms in a partial sum for approximating the infinite series. Default is
#'   \code{5,000}. 
#' @param delta Truncation value used to finish the computation of the partial sum. Default is \code{1e-29}.
#' @param nPoints Number of points used for integral approximation. Since we use \code{AGH} per default, this
#'   Defaults to \code{5}. 
#' @param method Integral approximation method used during ML estimation. Default is \code{AGH}, currently the 
#'   only supported value.
#' @param qmcType A character string indicating the type of points used for the QMC optimization algorithm.
#'   Available are \code{"Halton"} or \code{"Sobel"}. Default is \code{"Halton"}. Currently not in use.
#' @param maxit_out Maximum number of iterations for the outer optimization. Default is \code{50}.
#' @param maxit_inn Maximum number of iterations for the inner optimization. Default is \code{500}.
#' @param x_tol_out Convergence criterion for the parameter estimates.
#'   Default is \code{10e-6}.
#' @param fx_tol_out Convergence criterion for the deviance. Default is
#'   \code{1e-6}.
#' @param verbose_out A logical indicating whether the steps in optimization
#'   should be printed to the console. Default is \code{FALSE}. This argument is used internally. Use verbose 
#'   instead!
#' @param verbose_inn A logical indicating whether the steps inside the
#'   optimization should be shown. Default is \code{FALSE}.
#' @param lambda Should there be a smoothing of the parameter vector between subsequent outer optimization runs.
#'   Can be used to stabilize the estimation. Default is \code{0}.
#' @param lambda_min When lambda is not zero use this as the minimum value. Default is \code{0}. A sensible schema 
#'   is lambda = 0.5 and lambda_min = 0.5.
#' @param n_threads Number of cores used during hierarchical ML estimation. If \code{n_threads = 1}, the default, 
#'   no parallelization is used. If \code{n_threads = 0} (!), the maximum number of cores is used.
#' @param biter Total number of MCMC iterations, including burn-in. Default is \code{5,000}.
#' @param burnin Number of initial iterations discarded as burn-in. Must be smaller than \code{biter}.
#'   Default is \code{2,500}.
#' @param nchain Number of parallel MCMC chains to run. Default is \code{1}. 
#' @param tau2 Variance of the random-walk proposal used when \code{type_proposal="standard"}.
#' Not used for \code{type_proposal="mixture"} or \code{type_proposal="pmwg"}. Default is \code{0.1}.
#' @param tau2_A: Used for the proposal distribution for the single-person four-parameter DDM.
#' @param tau2_B: Used for the additional variability parameters of the seven-parameter DDM. Currently not used.
#' @param use_adapt_dao A logical indicating whether the adaptive particle-proposal distribution described
#'   in Dao et al. (2025) is used for \code{type_proposal="pmwg"}. If \code{TRUE}, subject-specific
#'   proposal parameters are periodically re-estimated from the post-burn-in draws once adaptation becomes
#'   possible; if \code{FALSE}, the non-adaptive two-component mixture proposal (governed by \code{pi_mix} and
#'   \code{epsilon}) is used throughout. Default is \code{TRUE}.
#' @param pi_mix Mixture weight for the non-adaptive PMwG proposal (used when \code{use_adapt_dao=FALSE}, or
#'   before adaptation starts): with probability \code{pi_mix} a particle is drawn locally around the current
#'   value (scaled by \code{epsilon}), and with probability \code{1-pi_mix} from the group-level distribution
#'   \eqn{N(\mu_\alpha, \Sigma_\alpha)}. Also used for \code{type_proposal="mixture"}. Default is \code{0.85}.
#' @param pi_mix1 Mixture weight on the subject-specific adapted proposal once
#'   \code{use_adapt_dao=TRUE} adaptation is active (Dao et al., 2025). Default is \code{0.65}.
#' @param pi_mix2 Mixture weight on the random-walk proposal centered at the current value 
#'   once adaptation is active. Default is \code{0.30}.
#' @param pi_mix3 Mixture weight on the group-level prior draw once adaptation is active; a small weight here keeps #'   the importance weights bounded (Hesterberg, 1995). Default is \code{0.05}.
#' @param epsilon Scaling factor applied to the covariance of the locally-centered proposal, both in the
#'   non-adaptive PMwG mixture and in \code{type_proposal="mixture"}. Default is \code{0.10}.
#' @param R Number of particles drawn per subject in the conditional Monte Carlo step of the PMwG sampler.
#'   Only used for \code{type_proposal="pmwg"}. Default is \code{20}.
#' @param type_proposal The algorithm used to update the subject-level parameters \eqn{\alpha_i}. Default is
#'   \code{pmwg}, the Particle Metropolis within Gibbs sampler of Dao et al. (2025). Alternatives
#'   are \code{mixture}, a single-particle Metropolis-Hastings step with a two-component mixture proposal, and
#'   \code{standard}, a simple random-walk Metropolis-Hastings step (see \code{tau2}).
#' @param type_sigma_prior The prior used for the group-level covariance matrix \eqn{\Sigma_\alpha}. Default is
#'   \code{huang_wand}, the marginally noninformative prior of Huang & Wand (2013) described in Dao et al.
#'   (2025); its scale matrix is constructed internally from auxiliary variables \eqn{a_d}, so
#'   only \code{nu0} is used (\code{S0} is ignored). The alternative, \code{informative}, uses a standard
#'   Inverse-Wishart prior \eqn{IW(\code{nu0}, \code{S0})} with user-supplied hyperparameters.
#' @param type_inits The random effcets are initialized as random draws from a multivariate normal distribution
#'   \code{"none"} or with the \code{"ml"} estimates of the persons. Default is \code{"none"}.
#' @param m Prior mean vector for the group-level mean \eqn{\mu_\alpha}
#'   If \code{NULL} (default), a vector of zeros is used.
#' @param M Prior covariance matrix for the group-level mean \eqn{\mu_\alpha}. If \code{NULL} (default),
#'   the identity matrix is used.
#' @param nu0 Degrees-of-freedom hyperparameter for the prior on \eqn{\Sigma_\alpha}. Its meaning depends on
#'   \code{type_sigma_prior}: under \code{"huang_wand"} it is the Huang & Wand (2013) hyperparameter
#'   (\code{nu0=2} corresponds to the marginally noninformative specification of Dao et al. (2025)
#'   under \code{"informative"} it is the degrees of freedom of the Inverse-Wishart prior
#'   If \code{NULL} (default), an internal default is chosen depending on \code{type_sigma_prior}.
#' @param S0 Scale matrix hyperparameter for the prior on \eqn{\Sigma_\alpha}. Only used when
#'   \code{type_sigma_prior="informative"}; ignored under \code{"huang_wand"}. If \code{NULL} (default),
#'   an internal default is chosen.
#' @param control_nlminb A list of control parameters for
#'   \code{\link{nlminb}}. See the base function for details.
#'
#' @return A list with components for each of the possible arguments.
#'
#' @export

ddm_args <- function( 
    # general args:  
    type_ddm = "std", # alternative: navfuss
    type_alpha = "dao", # alternative: logit
    type_sigma = "cholesky",
    use_tmb_modes = TRUE,
    use_rcpp = TRUE, 
    use_lan = FALSE,
    use_tf = FALSE,
    kmax = 5000L,
    delta = 1e-29,
    # max-lik specific agruments:
    nPoints = 5,
    method = "AGH",
    qmcType = "Sobol",
    isSeed = 123,
    maxit_out = 50,
    maxit_inn = 500,
    x_tol_out = 10e-6,
    fx_tol_out = 10e-6,
    verbose_out = FALSE, 
    verbose_inn = FALSE,
    lambda = 0,
    lambda_min = 0,
    n_threads = 1L,
    # bayes specific arguments:
    biter = 5000, 
    burnin = 2500, 
    nchain = 1, 
    tau2 = 0.1,
    tau2_A = c(0.50, 0.05, 0.05, 0.01)^2,
    tau2_B = c(0.10, 0.10, 0.05)^2, 
    use_adapt_dao = TRUE,
    pi_mix = 0.85, 
    pi_mix1= 0.65,
    pi_mix2= 0.30,
    pi_mix3= 0.05,
    epsilon = 0.10, 
    R = 20,
    type_proposal = "pmwg", # pmwg, mixture, standard
    type_sigma_prior = "huang_wand", # informative
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

    if ( !( type_sigma_prior %in% c("huang_wand", "informative") ) ) {
        stop("Prior type for Sigma_alpha not available.") }

    if ( type_proposal %in% c("mixture","standard") ) {
        if ( !missing(type_sigma_prior) && type_sigma_prior == "huang_wand" ) {
            warning("type_sigma_prior='huang_wand' is not supported for type_proposal='",
               type_proposal, "'; we use 'informative'.")
        }
        type_sigma_prior <- "informative"
    }

    if ( type_sigma_prior == "huang_wand" && !is.null(S0) ) {
        warning("S0 is ignored when type_sigma_prior='huang_wand' ignoriert.")
    }
    
    if ( type_sigma_prior == "huang_wand" && !is.null(nu0) && nu0 != 2 ) {
        message("Dao et al. suggest to use nu0 =2, your nu0 != 2.")
    }

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
        type_sigma_prior = type_sigma_prior,
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
        use_tf = use_tf,
        kmax = kmax,
        delta = delta,
        maxlik_list = maxlik_list,
        bayes_list = bayes_list,
        nlminb_list = nlminb_list )
	return( args )
}