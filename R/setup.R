#' Set prior parameters for Truncated Normal prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param mu_p see `Mu_p`
#' @param Mu_p mean for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `mu_p`
#' @param sigmasq_p see `Sigmasq_p`
#' @param Sigmasq_p variance for the truncated normal prior on `P`, matrix
#' size K x N. Defaults to all same value `sigmasq_p`
#' @param Covar_p covariance for the truncated MVN prior on `P`, matrix
#' size K x N
#' @param mu_e see `Mu_e`
#' @param Mu_e mean for the truncated normal prior on `E`, matrix
#' size N x G. Defaults to all same value `mu_e`
#' @param sigmasq_e see `Sigmasq_e`
#' @param Sigmasq_e variance for the truncated normal prior on `P`, matrix
#' size N x G. Defaults to all same value `sigmasq_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#'
#' @return named list of prior parameters
#' @export
set_truncnorm_prior_parameters <- function(
        Theta,
        dims,
        mu_p = sqrt(100/dims$N),
        Mu_p = matrix(mu_p, nrow = dims$K, ncol = dims$N),
        sigmasq_p = mu_p, #mu_p/10,
        Sigmasq_p = matrix(sigmasq_p, nrow = dims$K, ncol = dims$N),
        Covar_p = NULL,
        mu_e = sqrt(100/dims$N), #mean(colSums(M))/(N*100),
        Mu_e = matrix(mu_e, nrow = dims$N, ncol = dims$G),
        sigmasq_e = mu_e, #mu_e/10,
        Sigmasq_e = matrix(sigmasq_e, nrow = dims$N, ncol = dims$G),
        alpha = 1, # 1
        Alpha = rep(alpha, dims$K),
        beta = 1, # 1
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    fill_list(Theta, list(
        Mu_p = Mu_p,
        Sigmasq_p = Sigmasq_p,
        Mu_e = Mu_e,
        Sigmasq_e = Sigmasq_e,
        Covar_p = Covar_p,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
}

#' Set prior parameters for Exponential prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param lambda_p see `Lambda_p`
#' @param Lambda_p rate parameter for the exponential prior on `P`, matrix
#' size K x N. Defaults to all same value `lambda_p`
#' @param lambda_e see `Lambda_e`
#' @param Lambda_e rate parameter for the exponential prior on `E`, matrix
#' size N x G. Defaults to all same value `lambda_e`
#' @param alpha see `Alpha`
#' @param Alpha shape parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `alpha`
#' @param beta see `Beta`
#' @param Beta rate parameter for the inverse-gamma prior on `sigmasq`,
#' length K. Defaults to all same value `beta`
#'
#' @return named list of prior parameters
#' @export
set_exponential_prior_parameters <- function(
        Theta,
        dims,
        lambda_p = sqrt(dims$N/100),
        Lambda_p = matrix(lambda_p, nrow = dims$K, ncol = dims$N),
        lambda_e = sqrt(dims$N/100),
        Lambda_e = matrix(lambda_e, nrow = dims$N, ncol = dims$G),
        alpha = 0.1,
        Alpha = rep(alpha, dims$K),
        beta = 2,
        Beta = rep(beta, dims$K),
        a = 0.8,
        b = 0.8
) {
    fill_list(Theta, list(
        Lambda_p = Lambda_p,
        Lambda_e = Lambda_e,
        Alpha = Alpha,
        Beta = Beta,
        a = a,
        b = b
    ))
}

#' Set prior parameters for Gamma prior
#'
#' @param Theta list of parameters
#' @param dims named list of dimensions N, K, G
#' @param alpha_p see `Alpha_p`
#' @param Alpha_p shape parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `alpha_p`
#' @param beta_p see `Beta_p`
#' @param Beta_p rate parameter for the gamma prior on `P`, matrix
#' size K x N. Defaults to all same value `beta_p`
#' @param alpha_e see `Alpha_e`
#' @param Alpha_e shape parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `alpha_e`
#' @param beta_e see `Beta_e`
#' @param Beta_e rate parameter for the gamma prior on `E`, matrix
#' size N x G. Defaults to all same value `beta_e`
#' #'
#' @return named list of prior parameters
#' @export
set_gamma_prior_parameters <- function(
        Theta,
        dims,
        alpha_p = 10,
        Alpha_p = matrix(alpha_p, nrow = dims$K, ncol = dims$N),
        beta_p = sqrt(dims$N),
        Beta_p = matrix(beta_p, nrow = dims$K, ncol = dims$N),
        alpha_e = 10,
        Alpha_e = matrix(alpha_e, nrow = dims$N, ncol = dims$G),
        beta_e = sqrt(dims$N),
        Beta_e = matrix(beta_e, nrow = dims$N, ncol = dims$G),
        a = 0.8,
        b = 0.8
) {
    fill_list(Theta, list(
        Alpha_p = Alpha_p,
        Beta_p = Beta_p,
        Alpha_e = Alpha_e,
        Beta_e = Beta_e,
        a = a,
        b = b
    ))
}

fill_list <- function(list, fill_with) {
    for (name in names(fill_with)) {
        if (!(name %in% names(list))) {
            list[[name]] = fill_with[[name]]
        }
    }
    return(list)
}

#' initialize Theta
#'
#' @param likelihood string, one of c('poisson','normal')
#' @param prior string, one of c('exponential','truncnormal','gamma')
#' @param learn_A boolean, whether A will be learned or will be fixed
#' @param dims named list of dimensions
#' @param inits list of initial values, optional
#' @param prior_parameters list of named prior parameters, optional
#'
#' @return named list of initialized unknowns
#' @noRd
initialize_Theta <- function(
        likelihood, prior,
        learn_A, dims,
        inits = NULL,
        prior_parameters = NULL
) {
    # prior parameters
    Theta = prior_parameters
    if (prior == 'truncnormal') {
        Theta = set_truncnorm_prior_parameters(Theta, dims)
    } else if (prior == 'exponential') {
        Theta = set_exponential_prior_parameters(Theta, dims)
    } else if (prior == 'gamma') {
        Theta = set_gamma_prior_parameters(Theta, dims)
    }

    # signatures P
    if (is.null(inits$P)) {
        Theta$P <- matrix(nrow = dims$K, ncol = dims$N)
        for (k in 1:dims$K) {
            for (n in 1:dims$N) {
                if (prior == 'truncnormal') {
                    # non MVN case for P
                    if (is.null(Theta$Covar_p)){
                        Theta$P[k,n] <- truncnorm::rtruncnorm(
                            1, a = 0, b = Inf,
                            mean = Theta$Mu_p[k,n],
                            sd = sqrt(Theta$Sigmasq_p[k,n])
                        )
                    }
                    # MVN for P
                    else{
                        mean_vector <- Theta$Mu_p[, n] # K x 1

                        lower <- rep(0, length(mean_vector))
                        upper <- rep(Inf, length(mean_vector))

                        # Perturb diagonal otherwise "sigma is not positive definite"
                        newsigma <- Theta$Covar_p + diag(1e-6, nrow(Theta$Covar_p))

                        newsigma <- lqmm::make.positive.definite(newsigma, tol=1e-3)

                        sample <- tmvtnorm::rtmvnorm(
                            1, mean = mean_vector, sigma = newsigma,
                            lower = lower, upper = upper
                        )

                        # sample every column of Theta$P from MVN
                        Theta$P[, n] <- sample
                    }

                } else if (prior == 'exponential') {
                    Theta$P[k,n] <- stats::rexp(1, Theta$Lambda_p[k,n])
                } else if (prior == 'gamma') {
                    Theta$P[k,n] <- stats::rgamma(1, Theta$Alpha_p[k,n], Theta$Beta_p[k,n])
                }
            }
        }
    } else {
        Theta$P <- inits$P
    }

    # exposures E
    if (is.null(inits$E)) {
        Theta$E <- matrix(nrow = dims$N, ncol = dims$G)
        for (n in 1:dims$N) {
            for (g in 1:dims$G) {
                if (prior == 'truncnormal') {
                    Theta$E[n,g] <- truncnorm::rtruncnorm(
                        1, a = 0, b = Inf,
                        mean = Theta$Mu_e[n,g],
                        sd = sqrt(Theta$Sigmasq_e[n,g])
                    )
                } else if (prior == 'exponential') {
                    Theta$E[n,g] <- stats::rexp(1, Theta$Lambda_e[n,g])
                } else if (prior == 'gamma') {
                    Theta$E[n,g] <- stats::rgamma(1, Theta$Alpha_e[n,g], Theta$Beta_e[n,g])
                }
            }
        }
    } else {
        Theta$E <- inits$E
    }

    # variance sigmasq
    if (likelihood == 'normal') {
        if (is.null(inits$sigmasq)) {
            # Theta$sigmasq <- sapply(1:dims$K, function(k) {
            #     invgamma::rinvgamma(n = 1, shape = Theta$Alpha[k], scale = Theta$Beta[k])
            #     # 1/rgamma(n = 1, shape = Theta$Alpha[k], rate = Theta$Beta[k])
            # })
            Theta$sigmasq <- armspp::arms(n_samples = dims$K, log_pdf = function(x) {-1*log(x)}, lower = 0, upper = 1000)
        } else {
            Theta$sigmasq <- inits$sigmasq
        }
    } else if (likelihood == 'poisson') {
        Theta$Z <- array(dim = c(dims$K, dims$N, dims$G))
        for (k in 1:dims$K) {
            for (g in 1:dims$G) {
                Theta$Z[k,,g] <- sample_Zkg_poisson(k, g, M, Theta, dims)
            }
        }
    }


    # signature assignment A
    if (learn_A) {
        Theta$q <- matrix(
            rbeta(dims$S * dims$N, Theta$a, Theta$b),
            nrow = dims$S, ncol = dims$N
        )
        Theta$A <- matrix(
            as.numeric(runif(dims$S * dims$N) < c(Theta$q)),
            nrow = dims$S, ncol = dims$N
        )
    } else if (!is.null(inits$A)) {
        Theta$A <- inits$A
        Theta$q <- inits$A
    } else {
        Theta$A <- matrix(1, nrow = dims$S, ncol = dims$N)
        Theta$q <- Theta$A
    }

    return(Theta)
}
