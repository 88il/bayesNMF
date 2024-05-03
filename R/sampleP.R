#' get mu and sigmasq for P[,n] in Normal-Exponential model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn_normal_exponential <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)
    sum_E_sq <- gamma * sum(Theta$E[n, ] ** 2)

    # compute mean
    mu_num_term_1 <- gamma * sweep(
        (M - Mhat_no_n), # dim KxG
        2, # multiply each row by E[n,]
        Theta$E[n, ], # length G
        "*"
    ) %>% # dim KxG
        rowSums() # length K

    mu_num_term_2 <- Theta$Lambda_p[, n] * Theta$sigmasq # length K
    mu_P <- (mu_num_term_1 - mu_num_term_2) / sum_E_sq # length K
    sigmasq_P <- Theta$sigmasq / sum_E_sq # length K

    return(list(
        mu = mu_P,
        sigmasq = sigmasq_P
    ))
}

#' get mu and sigmasq for P[,n] in Normal-TruncNormal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq or mu and covar
#' @noRd
get_mu_sigmasq_Pn_normal_truncnormal <- function(n, M, Theta, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, dims, n)

    if(is.null(Theta$Covar_p)){
        # compute mean
        mu_num_term_1 <- gamma * Theta$A[1,n] * (1/Theta$sigmasq) * (sweep(
            (M - Mhat_no_n), # dim KxG
            2, # multiply each row by E[n,]
            Theta$E[n, ], # length G
            "*"
        ) %>% # dim KxG
            rowSums()) # length K
        mu_num_term_2 <- Theta$Mu_p[, n] / Theta$Sigmasq_p[,n] # length K
        denom <- (1/Theta$Sigmasq_p[,n]) + gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2) / Theta$sigmasq


        mu_P <- (mu_num_term_1 + mu_num_term_2) / denom # length K
        sigmasq_P <- 1 / denom # length K

        return(list(
            mu = mu_P,
            sigmasq = sigmasq_P
        ))
    }

    # case when P has covariance matrix
    # NOTES: incorporate gamma and A matrix (from PAE) ?
    else{
        #' Create kxk covariance matrix for M using sigmasq values (Normal)
        #' @param Theta list of parameters
        #' @param dims list of dimensions
        #'
        #' @return matrix
        #' @noRd
        get_M_covariance_matrix <- function(Theta, dims){ # abstract this function somewhere else?
            covar_M <- matrix(0, dims$K, dims$K)
            diag(covar_M) <- Theta$sigmasq
            print(Theta$sigmasq)
            return(covar_M)
        }

        cov_p = Theta$Covar_p
        sigma_M <- get_M_covariance_matrix(Theta, dims) # KxK


        # if (kappa(sigma_M) > 1e12) {
        #     stop("Matrix is close to singular or poorly conditioned with a high kappa value.")
        # }

        # get rid of bottom 10% of vals
        cutoff = quantile(abs(cov_p), 0.1)
        cov_p[abs(cov_p) < cutoff] = 0

        # cutoff = quantile(abs(sigma_M), 0.1)
        # sigma_M[abs(sigma) < cutoff] = 0

        sigma_M_inv <- solve(sigma_M) # KxK
        covar_P_inv <- solve(cov_p) # KxK
        # sigma_M_inv <- solve(sigma_M) # KxK

        A_star <- sapply(1:dims$K, function(index) sum(Theta$E[n, ] * (M[k, ] - Mhat_no_n[k, ]))) # Kx1 vector

        A <- solve(covar_P_inv + sum(Theta$E[n, ] ** 2) * sigma_M_inv) # KxK
        B <- covar_P_inv %*% Theta$Mu_p[, n] - sigma_M_inv %*% A_star # Kx1

        A_inv <- solve(A)

        mu_P <- A_inv %*% B
        covar_P <- A_inv

        return(list(
            mu = mu_P,
            covar = covar_P
        ))
    }
}

#' sample P[,n] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_normal <- function(n, M, Theta, dims, prior = 'truncnormal', gamma = 1) {
    if (prior == 'truncnormal') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_truncnormal(n, M, Theta, dims, gamma = gamma)
    } else if (prior == 'exponential') {
        mu_sigmasq_P <- get_mu_sigmasq_Pn_normal_exponential(n, M, Theta, dims, gamma = gamma)
    }

    # non MVN case
    if (is.null(Theta$Covar_p)){
        mu_P = mu_sigmasq_P$mu
        sigmasq_P = mu_sigmasq_P$sigmasq

        # sample from truncated normal
        truncnorm::rtruncnorm(1, mean = mu_P, sd = sqrt(sigmasq_P), a = 0, b = Inf)
    }
    else{
        # mu_P = mu_sigmasq_P$mu # K x 1
        # covar_P = mu_sigmasq_P$covar # K x K
        #
        # # sample from multivariate truncated normal
        # # Lower and upper bounds set to 0 and Inf for all dimensions
        # lower <- rep(0, length(mu_P))
        # upper <- rep(Inf, length(mu_P))
        #
        # tmvtnorm::rtmvnorm(1, mean = mu_P, sigma = covar_P, lower = lower, upper = upper)

        mean_vector <- mu_sigmasq_P$mu[, n] # K x 1

        lower <- rep(0, length(mean_vector))
        upper <- rep(Inf, length(mean_vector))

        # Perturb diagonal otherwise "sigma is not positive definite"
        newsigma <- mu_sigmasq_P$covar # + diag(1e-6, nrow(Theta$Covar_p))

        sample <- tmvtnorm::rtmvnorm(
            1, mean = mean_vector, sigma = newsigma,
            lower = lower, upper = upper
        )

        # sample every column of Theta$P from MVN
        Theta$P[, n] <- sample
    }
}

#' sample P[,n] for Poisson likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_poisson <- function(n, M, Theta, dims, prior = 'gamma', gamma = 1) {
    if (prior == 'gamma') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, Theta$Alpha_p[k,n] + gamma * sum(Theta$Z[k,n,]), Theta$Beta_p[k,n] + gamma * sum(Theta$E[n,]))
        })
    } else if (prior == 'exponential') {
        sampled <- sapply(1:dims$K, function(k) {
            rgamma(1, 1 + gamma * sum(Theta$Z[k,n,]), Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,]))
        })
    }
    return(sampled)
}

#' Sample Pkn for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param k mutation type index
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pkn_poisson_exp <- function(k, n, M, Theta, gamma) {
    log_pdf <- function(Pkn) {
        gamma * sum(Theta$Z[k,n,]) * log(Pkn) -
            Pkn * (Theta$Lambda_p[k,n] + gamma * sum(Theta$E[n,]))
    }
    sample <- armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
    return(sample)
}

#' Sample Pn for Poisson-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pn_poisson_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$K, function(k) {
        sample_Pkn_poisson_exp(k, n, M, Theta, gamma)
    })
}

#' Sample Pkn for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param k mutation type index
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pkn_norm_exp <- function(k, n, M, Theta, gamma) {
    log_pdf <- function(Pkn) {
        Theta$P[k,n] <- Pkn
        Mhat <- get_Mhat(Theta)

        -Theta$Lambda_p[k,n] * Pkn -
            gamma/(2*Theta$sigmasq[k]) * sum((M[k,] - Mhat[k,]) ** 2)
    }
    sample <- armspp::arms(
        n_samples = 1,
        log_pdf = log_pdf,
        lower = 0,
        upper = 100
    )
    return(sample)
}

#' Sample Pn for Normal-Exponential Model with adaptive rejection metropolis sampling
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Pn_norm_exp <- function(n, M, Theta, dims, gamma) {
    sapply(1:dims$K, function(k) {
        sample_Pkn_norm_exp(k, n, M, Theta, gamma)
    })
}

#' sample P[,n] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','exponential')
#' @param prior string, one of c('gamma','exponential','truncnormal')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn <- function(n, M, Theta, dims, likelihood = 'normal', prior = 'truncnormal', gamma = 1) {
    if (likelihood == 'normal') {
        if (prior == 'truncnormal' | gamma > 0.5) {
            sample_Pn_normal(n, M, Theta, dims, prior, gamma)
        } else {
            sample_Pn_norm_exp(n, M, Theta, dims, gamma)
        }
    } else if (likelihood == 'poisson') {
        if (prior == 'gamma' | gamma > 0.5) {
            sample_Pn_poisson(n, M, Theta, dims, prior, gamma)
        } else {
            sample_Pn_poisson_exp(n, M, Theta, dims, gamma)
        }

    }
}