#' C-JAMP: Copula-based joint analysis of multiple phenotypes.
#'
#' Functions to perform C-JAMP: \code{\link{cjamp}} fits a joint model of two
#' phenotypes conditional on one or multiple predictors; \code{\link{cjamp_loop}}
#' uses \code{\link{cjamp}} to fit the same copula model separately for a list of
#' multiple predictors, e.g. for a genetic association study.
#'
#' Both functions \code{\link{cjamp}} and \code{\link{cjamp_loop}} fit a joint copula model
#' of two phenotypes using the Clayton copula or 2-parameter copula, conditional on none
#' or multiple predictors and covariates in the marginal models. The \code{\link[optimx]{optimx}}
#' function of the \code{optimx} package is used to maximize the log-likelihood (i.e. to minimize
#' the minus log-likelihood function \code{\link{minusloglik}}) to obtain
#' maximum likelihood coefficient estimates of all parameters.
#' For this, the BFGS optimization method is recommended. Standard error estimates
#' of the coefficient estimates can be obtained by using
#' the observed inverse information matrix. P-values from hypothesis tests of
#' the absence of effects of each predictor on each phenotype in the marginal
#' models can be obtained from large-sample Wald-type tests (see the vignette for details).
#'
#' It should be noted that \code{\link{cjamp}}, \code{\link{cjamp_loop}} and
#' \code{\link{minusloglik}}) assume quantitative predictors and use an
#' additive model, i.e. for categorical predictors with more than 2 levels,
#' dummy variables have to be created beforehand. Accordingly, if single
#' nucleotide variants (SNVs) are included as predictors, the computation is
#' based on an additive genetic model if SNVs are provided as 0-1-2 genotypes
#' and on a dominant model if SNVs are provided as 0-1 genotypes.
#'
#' The \code{\link{cjamp}} function returns point estimates of all parameters,
#' standard error estimates and p-values for all marginal parameters (i.e. all parameters
#' for \code{predictors_Y1}, \code{predictors_Y1}), the minus
#' log-likelihood value as well as information about the convergence. The
#' \code{\link{cjamp_loop}} function only returns point estimates, standard error
#' estimates, and p-values for the specified predictors \code{predictors} and not
#' the covariates \code{covariates_Y1} and \code{covariates_Y2}, in addition
#' to the minus log-likelihood value as well as convergence information.
#'
#' It is recommended that all variables are centered and scaled before the analysis,
#' which can be done through the \code{scale_var} parameter. Otherwise, if the scales
#' of the variables differ, it can lead to convergence problems of the optimization.
#'
#' @param copula String indicating whether the joint model will be computed
#'               under the Clayton (\code{"Clayton"}) or 2-parameter copula
#'               (\code{"2param"}) model.
#' @param Y1 Numeric vector containing the first phenotype.
#' @param Y2 Numeric vector containing the second phenotype.
#' @param predictors_Y1 Dataframe containing the predictors of \code{Y1}
#'                      in columns (for the \code{\link{cjamp}} function).
#' @param predictors_Y2 Dataframe containing the predictors of \code{Y2}
#'                      in columns (for the \code{\link{cjamp}} function).
#' @param scale_var Logical. Indicating whether all predictors will be centered and
#'                  scaled before the analysis or not (default: FALSE).
#' @param predictors Dataframe containing the predictors of \code{Y1}
#'                   and \code{Y2} in columns for which estimates are returned
#'                   (for the \code{\link{cjamp_loop}} function).
#' @param covariates_Y1 Dataframe containing the covariates of \code{Y1}
#'                      in columns for which estimates are not returned
#'                      (for the \code{\link{cjamp_loop}} function).
#' @param covariates_Y2 Dataframe containing the covariates of \code{Y2}
#'                      in columns for which estimates are not returned
#'                      (for the \code{\link{cjamp_loop}} function).
#' @param optim_method String passed to the \code{\link[optimx]{optimx}}
#'                     function. It specifies the optimization method in the
#'                     \code{\link[optimx]{optimx}} function that is used for
#'                     maximizing the log-likelihood function. Recommended is the
#'                     \code{"BFGS"} optimization method (default). For further
#'                     methods, see the description of the
#'                     \code{\link[optimx]{optimx}} function.
#' @param trace Integer passed to the \code{\link[optimx]{optimx}}
#'              function. It specifies the tracing information on the progress
#'              of the optimization. The default value \code{2} gives full tracing,
#'              value \code{0} blocks all details. See also the \code{\link[optimx]{optimx}}
#'              documentation.
#' @param kkt2tol Numeric. Passed to the \code{\link[optimx]{optimx}}
#'                function, default value is 1E-16. It specifies the tolerance for
#'                the eigenvalue ratio in the Karush-Kuhn-Tucker (KKT) test for a positive definite
#'                Hessian matrix. See also the \code{\link[optimx]{optimx}}
#'                documentation.
#' @param SE_est Logical indicator whether standard error estimates are
#'               computed for the parameters using the inverse of the observed
#'               information matrix (\code{TRUE}, default), or whether standard
#'               error estimates are not computed (\code{FALSE}).
#' @param pval_est Logical indicator whether p-values are computed from
#'                 hypothesis tests of the absence of effects of each predictor
#'                 on each phenotype in the marginal models (\code{TRUE}, default),
#'                 or whether p-values are not computed (\code{FALSE}). P-values
#'                 are obtained from large sample Wald-type tests based on the
#'                 maximum likelihood parameter estimates and the standard error
#'                 estimates.
#' @param n_iter_max Integer indicating the maximum number of optimization attempts
#'                   of the log-likelihood function with different starting values,
#'                   if the optimization doesn't converge (default: 10).
#' @return An object of class \code{cjamp}, for which the summary function
#'         \code{\link{summary.cjamp}} is implemented. The output is a list
#'         containing estimates of the copula parameters, marginal parameters, Kendall's
#'         tau (as well as the upper and lower tail dependence \eqn{\lambda_l,
#'         \lambda_u} if the 2-parameter copula model is fitted), the standard
#'         error estimates of all parameters, p-values of hypothesis tests
#'         of the marginal parameters (i.e. of the absence of predictor effects on the
#'         phenotypes), the convergence code of the log-likelihood maximization
#'         (from the \code{\link[optimx]{optimx}}
#'         function, where 0 indicates successful convergence), the KKT
#'         conditions 1 and 2 of the convergence, and the maximum log-likelihood value.
#'
#' @examples
#'
#' # Data generation
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 100)
#' phenodata <- generate_phenodata_2_copula(genodata = genodata$SNV1,
#'                                          MAF_cutoff = 1, prop_causal = 1,
#'                                          tau = 0.2, b1 = 0.3, b2 = 0.3)
#'
#' ## Not run. When executing, the following takes about 2 minutes running time.
#' ## Example 1: Analysis of multiple SNVs as predictors in one model
#' #predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#' #                         genodata[, 1:3])
#' #cjamp(copula = "Clayton", Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#' #      predictors_Y1 = predictors, predictors_Y2 = predictors,
#' #      optim_method = "BFGS", trace = 2, kkt2tol = 1E-16, SE_est = TRUE,
#' #      pval_est = TRUE, n_iter_max = 10)
#' #cjamp(copula = "2param", Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#' #      predictors_Y1 = predictors, predictors_Y2 = predictors,
#' #      optim_method = "BFGS", trace = 2, kkt2tol = 1E-16, SE_est = TRUE,
#' #      pval_est = TRUE, n_iter_max = 10)
#' #
#' ## Example 2: Analysis of multiple SNVs in separate models
#' #covariates <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' #predictors <- genodata
#' #cjamp_loop(copula = "Clayton", Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#' #           predictors = predictors, covariates_Y1 = covariates,
#' #           covariates_Y2 = covariates, optim_method = "BFGS", trace = 2,
#' #           kkt2tol = 1E-16, SE_est = TRUE, pval_est = TRUE,
#' #           n_iter_max = 10)
#'
#' @export
#'

cjamp <- function(copula = "Clayton", Y1 = NULL, Y2 = NULL, predictors_Y1 = NULL,
                  predictors_Y2 = NULL, scale_var = FALSE, optim_method = "BFGS",
                  trace = 2, kkt2tol = 1e-16, SE_est = TRUE, pval_est = TRUE,
                  n_iter_max = 10) {
    if (pval_est & !SE_est) {
        stop("SE_est has to be TRUE to compute p-values.")
    }
    if (!requireNamespace("optimx", quietly = TRUE)) {
        stop("Package optimx needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (is.null(Y1) | is.null(Y2)) {
        stop("Y1 and Y2 have to be supplied.")
    }
    if (stats::cor(Y1, Y2, use = "complete.obs", method = "kendall") < 0) {
        Y2 <- -Y2
        warning("Dependence between Y1, Y2 is negative but copulas require positive dependence. \n
                 Hence Y2 is transformed to -Y2 and point estimates for Y2 are in inverse direction.")
    }
    if (copula == "Clayton") {
        n_dep <- 1
        copula_param = "phi"
    }
    if (copula == "2param") {
        n_dep <- 2
        copula_param = "both"
    }
    if ((!class(predictors_Y1) == "data.frame" & !is.null(predictors_Y1)) |
        (!class(predictors_Y2) == "data.frame" & !is.null(predictors_Y2))) {
         stop("predictors_Y1 and predictors_Y2 have to be a dataframe or NULL.")
    }
    if (class(predictors_Y1) == "data.frame") {
        if (!all(sapply(predictors_Y1, is.numeric))) {
            predictors_Y1 <- as.data.frame(data.matrix(predictors_Y1))
            warning("predictors_Y1 contains non-numeric Variables,
                     which have automatically been transformed.")
        }
    }
    if (class(predictors_Y2) == "data.frame") {
        if (!all(sapply(predictors_Y2, is.numeric))) {
            predictors_Y2 <- as.data.frame(data.matrix(predictors_Y2))
            warning("predictors_Y2 contains non-numeric Variables,
                     which have automatically been transformed.")
        }
    }
    if (is.null(predictors_Y1) & is.null(predictors_Y2)) {
        idx <- (!is.na(Y1) & !is.na(Y2))
        Y1 <- Y1[idx]
        Y2 <- Y2[idx]
        n_pred <- 0
    }
    if (!is.null(predictors_Y1) & is.null(predictors_Y2)) {
        idx <- (!is.na(Y1) & !is.na(Y2) & stats::complete.cases(predictors_Y1))
        Y1 <- Y1[idx]
        Y2 <- Y2[idx]
        predictors_Y1 <- predictors_Y1[idx, , drop = FALSE]
        if (qr(predictors_Y1)$rank < dim(predictors_Y1)[2]) {
            stop("Complete cases matrix of predictors_Y1 doesn't have full rank.")
        }
        n_pred <- dim(predictors_Y1)[2]
    }
    if (is.null(predictors_Y1) & !is.null(predictors_Y2)) {
        idx <- (!is.na(Y1) & !is.na(Y2) & stats::complete.cases(predictors_Y2))
        Y1 <- Y1[idx]
        Y2 <- Y2[idx]
        predictors_Y2 <- predictors_Y2[idx, , drop = FALSE]
        if (qr(predictors_Y2)$rank < dim(predictors_Y2)[2]) {
            stop("Complete cases matrix of predictors_Y2 doesn't have full rank.")
        }
        n_pred <- dim(predictors_Y2)[2]
    }
    if (!is.null(predictors_Y1) & !is.null(predictors_Y2)) {
        idx <- (!is.na(Y1) & !is.na(Y2) & stats::complete.cases(predictors_Y1) &
                stats::complete.cases(predictors_Y2))
        Y1 <- Y1[idx]
        Y2 <- Y2[idx]
        predictors_Y1 <- predictors_Y1[idx, , drop = FALSE]
        predictors_Y2 <- predictors_Y2[idx, , drop = FALSE]
        if (qr(predictors_Y1)$rank < dim(predictors_Y1)[2]) {
            stop("Complete cases matrix of predictors_Y1 doesn't have full rank.")
        }
        if (qr(predictors_Y2)$rank < dim(predictors_Y2)[2]) {
            stop("Complete cases matrix of predictors_Y2 doesn't have full rank.")
        }
        n_pred <- dim(predictors_Y1)[2] + dim(predictors_Y2)[2]
    }
    if (scale_var) {
        predictors_Y1 <- data.frame(lapply(predictors_Y1, scale))
        predictors_Y2 <- data.frame(lapply(predictors_Y2, scale))
    }
    convcode <- 1
    kkt <- c(FALSE, FALSE)
    n_iter <- 1
    helpnum <- 0.4
    while (!(convcode == 0 & all(kkt) & !any(is.na(kkt))) & (n_iter <= n_iter_max)) {
        startvalues <- get_estimates_naive(Y1 = Y1, Y2 = Y2,
                                           predictors_Y1 = predictors_Y1,
                                           predictors_Y2 = predictors_Y2,
                                           copula_param = copula_param) - helpnum
        names_pred <- names(startvalues[(n_dep + 3):length(startvalues)])
        if (minusloglik(Y1 = Y1, Y2 = Y2, predictors_Y1 = predictors_Y1,
                        predictors_Y2 = predictors_Y2, copula = copula,
                        parameters = startvalues) == Inf |
            is.na(minusloglik(Y1 = Y1, Y2 = Y2, predictors_Y1 = predictors_Y1,
                              predictors_Y2 = predictors_Y2, copula = copula,
                              parameters = startvalues))) {
            helpnum <- (-1)^stats::rbinom(1, 1, 0.5) * stats::runif(1, min = 0.1, max = 1)
            n_iter <- n_iter + 1
            next
        }
        res <- optimx::optimx(par = startvalues, fn = minusloglik, Y1 = Y1, Y2 = Y2,
                              predictors_Y1 = predictors_Y1, predictors_Y2 = predictors_Y2,
                              copula = copula, method = optim_method, hessian = TRUE,
                              control = list(trace = trace, kkt2tol = kkt2tol))
        convcode <- res$convcode
        kkt <- c(res$kkt1, res$kkt2)
        helpnum <- (-1)^stats::rbinom(1, 1, 0.5) * stats::runif(1, min = 0.1, max = 1)
        n_iter <- n_iter + 1
    }
    names(convcode) <- "convcode"
    names(kkt) <- c("KKT1", "KKT2")
    maxloglik <- NULL
    if (!convcode == 0) {
        warning("In model with predictors \n", names_pred,
                ": \n Optimization doesn't seem to be converged.",
                "\n SE estimates and pvalues are not computed. Please check.")
    }
    if (!all(kkt) | any(is.na(kkt))) {
        warning("In model with predictors \n", names_pred,
                ": \n KKT conditions are not satistified.",
                " \n SE estimates and pvalues are not computed. Please check.")
    }
    if (!convcode == 0 | !all(kkt) | any(is.na(kkt))) {
        warning("In model with predictors \n", names_pred,
                ": \n Naive parameter estimates are computed from separate marginal models.")
        point_estimates <- get_estimates_naive(Y1 = Y1, Y2 = Y2, predictors_Y1 = predictors_Y1,
                                               predictors_Y2 = predictors_Y2,
                                               copula_param = copula_param)
        point_estimates <- as.list(point_estimates)
        if (copula == "Clayton") {
            phi <- exp(point_estimates$log_phi)
            tau <- phi/(2 + phi)
            theta <- lambda_l <- lambda_u <- NA
        }
        if (copula == "2param") {
            phi <- exp(point_estimates$log_phi)
            theta <- exp(point_estimates$log_theta_minus1) + 1
            tau <- 1 - (2/(theta * (phi + 2)))
            lambda_l <- 2^(-1/(theta * phi))
            lambda_u <- 2 - 2^(1/theta)
        }
        point_estimates <- c(point_estimates[(n_dep + 1):length(point_estimates)],
                             phi, theta, tau, lambda_l, lambda_u)
        point_estimates[[1]] <- exp(point_estimates[[1]])
        point_estimates[[2]] <- exp(point_estimates[[2]])
        names(point_estimates) <- c("Y1_sigma", "Y2_sigma", names_pred, "phi",
                                    "theta", "tau", "lambda_l", "lambda_u")
    }
    if (convcode == 0 & all(kkt) & !any(is.na(kkt))) {
        maxloglik <- res$value
        names(maxloglik) <- "maxloglik"
        point_estimates <- res[(n_dep + 1):(n_dep + n_pred + 4)]
        point_estimates[1:2] <- exp(point_estimates[1:2])
        if (copula == "Clayton") {
            phi <- exp(res$log_phi)
            tau <- phi/(2 + phi)
            theta <- lambda_l <- lambda_u <- NA
        }
        if (copula == "2param") {
            phi <- exp(res$log_phi)
            theta <- exp(res$log_theta_minus1) + 1
            tau <- 1 - (2/(theta * (phi + 2)))
            lambda_l <- 2^(-1/(theta * phi))
            lambda_u <- 2 - 2^(1/theta)
        }
        point_estimates <- c(point_estimates, phi, theta, tau, lambda_l, lambda_u)
        names(point_estimates) <- c("Y1_sigma", "Y2_sigma", names_pred,
                                    "phi", "theta", "tau", "lambda_l", "lambda_u")
    }
    SE_estimates <- NULL
    if (SE_est & convcode == 0 & all(kkt) & !any(is.na(kkt))) {
        SE_estimates <- vector(mode = "numeric", length = (n_pred + 4))
        SE_estimates[1] <- sqrt(diag(solve(attributes(res)$details[[3]]))[n_dep + 1] *
                                  (point_estimates$Y1_sigma^2))
        SE_estimates[2] <- sqrt(diag(solve(attributes(res)$details[[3]]))[n_dep + 2] *
                                  (point_estimates$Y2_sigma^2))
        SE_estimates[3:length(SE_estimates)] <-
          sqrt(diag(solve(attributes(res)$details[[3]]))[(n_dep + 3):(n_dep + n_pred + 4)])
        if (copula == "Clayton") {
            helpvar1 <- 2 * phi/((phi + 2)^2)
            phi_SE <- sqrt(diag(solve(attributes(res)$details[[3]]))[1] * (phi^2))
            tau_SE <- sqrt(diag(solve(attributes(res)$details[[3]]))[1] * (helpvar1^2))
            theta_SE <- lambda_l_SE <- lambda_u_SE <- NA
        }
        if (copula == "2param") {
            helpvar3.1 <- 2 * phi/(theta * (phi + 2)^2)
            helpvar3.2 <- 2 * (theta - 1)/(theta^2 * (phi + 2))
            helpvar3.3 <- (2^(-1/(theta * phi))) * (log(2)/(theta * phi))
            helpvar3.4 <- (2^(-1/(theta * phi))) * (log(2) * (theta - 1)/(theta^2 * phi))
            helpvar3.5 <- (2^(1/theta)) * ((theta - 1) * log(2)/theta^2)
            phi_SE <- sqrt(diag(solve(attributes(res)$details[[3]]))[1] * (phi^2))
            theta_SE <- sqrt(diag(solve(attributes(res)$details[[3]]))[2] *
                               ((theta - 1)^2))
            tau_SE <- sqrt(t(c(helpvar3.1, helpvar3.2)) %*%
                             solve(attributes(res)$details[[3]])[1:2, 1:2] %*%
                             c(helpvar3.1, helpvar3.2))
            lambda_l_SE <- sqrt(t(c(helpvar3.3, helpvar3.4)) %*%
                                  solve(attributes(res)$details[[3]])[1:2, 1:2] %*%
                                  c(helpvar3.3, helpvar3.4))
            lambda_u_SE <- sqrt(diag(solve(attributes(res)$details[[3]]))[2] *
                                  helpvar3.5^2)
        }
        SE_estimates <- c(SE_estimates, phi_SE, theta_SE, tau_SE,
                          lambda_l_SE, lambda_u_SE)
        names(SE_estimates) <- c("Y1_sigma", "Y2_sigma", names_pred, "phi",
                                 "theta", "tau", "lambda_l", "lambda_u")
    }
    if (is.null(SE_estimates)){
        SE_estimates <- rep(NA, length(point_estimates))
    }
    pval <- NULL
    point_estimates <- unlist(point_estimates)
    if (pval_est & convcode == 0 & all(kkt) & !any(is.na(kkt))) {
        pval <- vector(mode = "numeric", length = (n_pred + 2))
        pval <- 2 * stats::pnorm(as.numeric(-abs(point_estimates[3:(n_pred + 4)]/
                                            SE_estimates[3:(n_pred + 4)])))
        names(pval) <- names_pred
    }
    if (is.null(pval)){
        pval <- rep(NA, (n_pred + 2))
    }
    output <- list(point_estimates = point_estimates,
                   SE_estimates = SE_estimates,
                   pval = pval, convcode = convcode,
                   kkt = kkt, maxloglik = maxloglik)
    names(output) <- c("Parameter point estimates",
                       "Parameter standard error estimates",
                       "Parameter p-values",
                       "Convergence code of optimx function",
                       "Karush-Kuhn-Tucker conditions 1 and 2",
                       "Maximum log-likelihood")
   class(output) <- "cjamp"
   return(output)
}

#' @rdname cjamp
#' @export

cjamp_loop <- function(copula = "Clayton", Y1 = NULL, Y2 = NULL, predictors = NULL,
                       covariates_Y1 = NULL, covariates_Y2 = NULL, scale_var = FALSE,
                       optim_method = "BFGS", trace = 0, kkt2tol = 1e-16,
                       SE_est = TRUE, pval_est = TRUE, n_iter_max = 10) {
    if (is.null(Y1) | is.null(Y2)) {
        stop("Y1 and Y2 have to be supplied.")
    }
    if (is.null(predictors)) {
        stop("At least one predictor has to be supplied.")
    }
    if (!class(predictors) == "data.frame") {
        stop("Predictors has to be a dataframe.")
    }
    if (!all(sapply(predictors, is.numeric))) {
        predictors <- as.data.frame(data.matrix(predictors))
        warning("predictors contains non-numeric Variables,
                 which have automatically been transformed.")
    }
    if ((!class(covariates_Y1) == "data.frame" & !is.null(covariates_Y1)) |
        (!class(covariates_Y2) == "data.frame" & !is.null(covariates_Y2))) {
        stop("covariates_Y1 and covariates_Y2 have to be a dataframe or NULL.")
    }
    if (class(covariates_Y1) == "data.frame") {
        if (!all(sapply(covariates_Y1, is.numeric))) {
            covariates_Y1 <- as.data.frame(data.matrix(covariates_Y1))
            warning("covariates_Y1 contains non-numeric Variables,
                     which have automatically been transformed.")
        }
    }
    if (class(covariates_Y2) == "data.frame") {
        if (!all(sapply(covariates_Y2, is.numeric))) {
            covariates_Y2 <- as.data.frame(data.matrix(covariates_Y2))
            warning("covariates_Y2 contains non-numeric Variables,
                     which have automatically been transformed.")
        }
    }
    N <- dim(predictors)[2]
    Y1_point_estimates <- Y1_SE_estimates <- Y1_pval <- NULL
    Y2_point_estimates <- Y2_SE_estimates <- Y2_pval <- NULL
    convcode <- kkt <- maxloglik <- NULL
    for (i in 1:N) {
         names_pred_i <- names(predictors)[i]
         predictors_Y1 <- predictors_Y2 <- predictors[, i, drop = FALSE]
         if (!is.null(covariates_Y1)) {
             predictors_Y1 <- cbind(predictors_Y1, as.data.frame(covariates_Y1))
         }
         if (!is.null(covariates_Y2)) {
             predictors_Y2 <- cbind(predictors_Y2, as.data.frame(covariates_Y2))
         }
         n_pred_Y1 <- dim(predictors_Y1)[2]
         n_pred_Y2 <- dim(predictors_Y2)[2]
         results <- cjamp(Y1 = Y1, Y2 = Y2, predictors_Y1 = predictors_Y1,
                          predictors_Y2 = predictors_Y2, scale_var = scale_var,
                          copula = copula, optim_method = optim_method,
                          trace = trace, kkt2tol = kkt2tol, SE_est = SE_est,
                          pval_est = pval_est, n_iter_max = n_iter_max)
         Y1_point_estimates[i] <- as.numeric(results[[1]][4])
         Y2_point_estimates[i] <- as.numeric(results[[1]][5 + n_pred_Y1])
         Y1_SE_estimates[i] <- results[[2]][4]
         Y2_SE_estimates[i] <- results[[2]][5 + n_pred_Y1]
         Y1_pval[i] <- results[[3]][2]
         Y2_pval[i] <- results[[3]][3 + n_pred_Y1]
         names(Y1_point_estimates)[i] <- names(Y2_point_estimates)[i] <- names_pred_i
         names(Y1_SE_estimates)[i] <- names(Y2_SE_estimates)[i] <- names_pred_i
         names(Y1_pval)[i] <- names(Y2_pval)[i] <- names_pred_i
         convcode[i] <- results[[4]]
         kkt[i] <- all(results[[5]])
         maxloglik[i] <- results[[6]]
    }
    output <- list(Y1_point_estimates = Y1_point_estimates,
                   Y2_point_estimates = Y2_point_estimates,
                   Y1_SE_estimates = Y1_SE_estimates,
                   Y2_SE_estimates = Y2_SE_estimates,
                   Y1_pval = Y1_pval, Y2_pval = Y2_pval,
                   convcode = convcode, kkt = kkt,
                   maxloglik = maxloglik)
    names(output) <- c("Parameter point estimates for effects of predictors on Y1",
                       "Parameter point estimates for effects of predictors on Y2",
                       "Parameter standard error estimates for effects on Y1",
                       "Parameter standard error estimates for effects on Y2",
                       "Parameter p-values for effects of predictors on Y1",
                       "Parameter p-values for effects of predictors on Y2",
                       "Convergence code of optimx function",
                       "Karush-Kuhn-Tucker conditions 1 and 2",
                       "Maximum log-likelihood")
    class(output) <- "cjamp"
    return(output)
}
