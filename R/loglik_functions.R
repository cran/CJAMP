#' Minus log-likelihood of copula models.
#'
#' Function to compute the minus log-likelihood of the joint distribution
#' of \code{Y1} and \code{Y2} given the predictors \code{predictors_Y1}
#' of \code{Y1} and \code{predictors_Y2} of \code{Y2} in a copula model
#' for given parameter values \code{parameters}. Implemented are the
#' Clayton and 2-parameter copula. The function assumes quantitative
#' predictors and uses an additive model, i.e. for categorical predictors
#' with more than 2 levels, dummy variables have to be created beforehand.
#' Accordingly, if single nucleotide variants (SNVs) are included as
#' predictors, the computation is based on an additive genetic model if
#' SNVs are provided as 0-1-2 genotypes and on a dominant model if SNVs
#' are provided as 0-1 genotypes.
#'
#' The number of predictors is not fixed and also different predictors
#' can be supplied for \code{Y1} and \code{Y2}. However, for the functioning
#' of the log-likelihood function, the \code{parameters} vector has to be
#' supplied in a pre-specified form and formatting:
#'
#' The vector has to include the copula parameter(s), marginal parameters
#' for \code{Y1}, and marginal parameters for \code{Y2} in this order.
#' For example, for the 2-parameter copula with parameters \eqn{\phi, \theta}
#' and marginal models
#'
#' \deqn{Y_1 = \alpha_0 + \alpha_1 \cdot X_1 + \alpha_2 \cdot X_2 + \epsilon_1, \ \epsilon_1 \sim N(0,\sigma_1^2),}{Y1 = \alpha0 + \alpha1*X1 + \alpha2*X2 + \epsilon1, \epsilon1 ~ N(0,\sigma1^2),}
#' \deqn{Y_1 = \beta_0 + \beta_1 \cdot X_1 + \epsilon_2, \ \epsilon_2 \sim N(0,\sigma_2^2),}{Y1 = \beta0 + \beta1*X1 + \epsilon2, \epsilon2 ~ N(0,\sigma2^2),}
#'
#' the parameter vector has be
#' \deqn{(log(\phi), log(\theta-1), log(\sigma_1), log(\sigma_2), \alpha_0, \alpha_2, \alpha_1, \beta_0, \beta_1)^T.}{(log(\phi), log(\theta-1), log(\sigma1), log(\sigma2), \alpha0, \alpha2, \alpha1, \beta0, \beta1).}
#' \eqn{log(\phi)} and \eqn{log(\theta-1)} have to be provided instead of
#' \eqn{\phi, \theta} and \eqn{log(\sigma_1)}{log(\sigma1)}, \eqn{log(\sigma_2)}{log(\sigma2)} instead of
#' \eqn{\sigma_1}{\sigma1}, \eqn{\sigma_2}{\sigma2} for computational reasons when the log-likelihood
#' function is optimized in \code{\link{cjamp}}. As further necessary format,
#' the vector has to be named and the names of the copula parameter(s) has/have
#' to be \code{log_phi} and/or \code{log_theta_minus1}.
#'
#' @param copula String indicating whether the likelihood should be computed
#'               under the Clayton (\code{"Clayton"}) or 2-parameter copula
#'               (\code{"2param"}) model.
#' @param Y1 Numeric vector containing the first phenotype.
#' @param Y2 Numeric vector containing the second phenotype.
#' @param predictors_Y1 Dataframe containing the predictors of \code{Y1}
#'                      in columns.
#' @param predictors_Y2 Dataframe containing predictors of \code{Y2}
#'                      in columns.
#' @param parameters Named integer vector containing the parameter estimates
#'                   of the marginal and copula parameters, for which the
#'                   log-likelihood will be computed. For details and the
#'                   necessary format of the vector, see the details below.
#' @return  Minus log-likelihood value (integer).
#'
#' @examples
#' # Generate genetic data:
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 1000)
#'
#' # Generate phenotype data:
#' phenodata <- generate_phenodata_2_copula(genodata = genodata, MAF_cutoff = 1,
#'                                          prop_causal = 0.5, tau = 0.2,
#'                                          b1 = 0.3, b2 = 0.3)
#'
#' # Example 1: Log-likelihood of null model without covariates & genetic effects
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = NULL, predictors_Y2 = NULL,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = NULL,
#'             predictors_Y2 = NULL, parameters = estimates, copula = "2param")
#'
#' # Example 2: Log-likelihood of null model with covariates, without genetic effects
#' predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors,
#'                                  predictors_Y2 = predictors,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors,
#'             predictors_Y2 = predictors, parameters = estimates, copula = "2param")
#'
#' # Example 3: Log-likelihood of model with covariates & genetic effect on Y1 only
#' predictors_Y1 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#'                             SNV = genodata$SNV1)
#' predictors_Y2 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors_Y1,
#'                                  predictors_Y2 = predictors_Y2,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors_Y1,
#'             predictors_Y2 = predictors_Y2, parameters = estimates,
#'             copula = "2param")
#'
#' # Example 4: Log-likelihood of model with covariates & genetic effect on Y2 only
#' predictors_Y1 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' predictors_Y2 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#'                             SNV = genodata$SNV1)
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors_Y1,
#'                                  predictors_Y2 = predictors_Y2,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors_Y1,
#'             predictors_Y2 = predictors_Y2, parameters = estimates,
#'             copula = "2param")
#'
#' # Example 5: Log-likelihood of model without covariates, with genetic effects
#' predictors <- data.frame(SNV = genodata$SNV1)
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors,
#'                                  predictors_Y2 = predictors,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors,
#'             predictors_Y2 = predictors, parameters = estimates, copula = "2param")
#'
#' # Example 6: Log-likelihood of model with covariates & genetic effects
#' predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2, SNV = genodata$SNV1)
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors,
#'                                  predictors_Y2 = predictors,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors,
#'             predictors_Y2 = predictors, parameters = estimates, copula = "2param")
#'
#' # Example 7: Log-likelihood of model with covariates & multiple genetic effects
#' predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2, genodata[, 1:5])
#' estimates <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                  predictors_Y1 = predictors,
#'                                  predictors_Y2 = predictors,
#'                                  copula_param = "both")
#' minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2, predictors_Y1 = predictors,
#'             predictors_Y2 = predictors, parameters = estimates, copula = "2param")
#'
#' @export
#'

minusloglik <- function(copula = "Clayton", Y1 = NULL, Y2 = NULL,
                        predictors_Y1 = NULL, predictors_Y2 = NULL,
                        parameters = NULL) {
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
    if (!is.null(predictors_Y2)) {
        n_Y2_pred <- dim(predictors_Y2)[2]
    } else {
        n_Y2_pred <- 0
    }
    if (!is.null(predictors_Y1)) {
        n_Y1_pred <- dim(predictors_Y1)[2]
    } else {
        n_Y1_pred <- 0
    }
    if (!is.null(predictors_Y2)) {
        n_Y2_pred <- dim(predictors_Y2)[2]
    } else {
        n_Y2_pred <- 0
    }
    if (copula == "Clayton") {
        param_Y1_Y2 <- parameters[1]
        n_dep <- 1
        if (!names(param_Y1_Y2) == "log_phi") {
            stop("Estimate for phi is not provided.")
        }
    }
    if (copula == "2param") {
        param_Y1_Y2 <- parameters[1:2]
        n_dep <- 2
        if (!identical(names(param_Y1_Y2), c("log_phi", "log_theta_minus1"))) {
            stop("Estimates for phi/theta are not provided.")
        }
    }
    param_Y1_sigma <- parameters[n_dep + 1]
    param_Y2_sigma <- parameters[n_dep + 2]
    no_p <- length(parameters)
    param_Y1_predictors <- parameters[(n_dep + 3):(n_dep + 3 + n_Y1_pred)]
    param_Y2_predictors <- parameters[(no_p - n_Y2_pred):no_p]
    if (is.null(Y1) | is.null(Y2)) {
        stop("Both Y1 and Y2 have to be supplied.")
    }
    n_ind <- length(Y1)
    if (!is.null(predictors_Y1)) {
        predictors_Y1 <- cbind(data.frame(One = rep(1, n_ind)), predictors_Y1)
    } else {
        predictors_Y1 <- data.frame(One = rep(1, n_ind))
    }
    if (!is.null(predictors_Y2)) {
        predictors_Y2 <- cbind(data.frame(One = rep(1, n_ind)), predictors_Y2)
    } else {
        predictors_Y2 <- data.frame(One = rep(1, n_ind))
    }
    if (!n_ind == length(Y2)) {
        stop("Variables must have same length.")
    }
    if (!n_ind == dim(predictors_Y1)[1]) {
        stop("Variables must have same length.")
    }
    if (!n_ind == dim(predictors_Y2)[1]) {
        stop("Variables must have same length.")
    }
    if (is.null(param_Y1_Y2) | any(is.na(param_Y1_Y2))) {
        stop("Estimate(s) of dependence betwee Y1 and Y2 is (are) not provided.")
    }
    if (any(is.na(param_Y1_sigma)) | any(is.na(param_Y2_sigma))) {
        stop("One or both marginal variances are not provided")
    }
    if (any(is.na(param_Y1_predictors))) {
        warning(paste("Effect estimate of ",
                      names(predictors_Y1)[which(is.na(param_Y1_predictors))],
                      " on Y1 is missing. Variable ",
                      names(predictors_Y1)[which(is.na(param_Y1_predictors))],
                      " is excluded from the analysis.", sep = ""))
        predictors_Y1 <- predictors_Y1[, !is.na(param_Y1_predictors)]
        param_Y1_predictors <- param_Y1_predictors[!is.na(param_Y1_predictors)]
    }
    if (any(is.na(param_Y2_predictors))) {
        warning(paste("Effect estimates of ",
                      names(predictors_Y2)[which(is.na(param_Y2_predictors))],
                      " on Y2 is missing. Variable ",
                      names(predictors_Y2)[which(is.na(param_Y1_predictors))],
                      " is excluded from the analysis.", sep = ""))
        predictors_Y2 <- predictors_Y2[, !is.na(param_Y2_predictors)]
        param_Y2_predictors <- param_Y2_predictors[!is.na(param_Y2_predictors)]
    }
    if (copula == "Clayton") {
        phi <- exp(param_Y1_Y2[1])
    } else if (copula == "2param") {
        phi <- exp(param_Y1_Y2[1])
        theta <- exp(param_Y1_Y2[2]) + 1
    } else {
        stop("copula has to be specified")
    }
    minusloglik <- 0
    predictors_Y1 <- as.matrix(predictors_Y1)
    predictors_Y2 <- as.matrix(predictors_Y2)
    for (i in 1:n_ind) {
        lik.pt1.1 <- stats::pnorm(Y1[i],
                           mean = as.numeric(param_Y1_predictors) %*% predictors_Y1[i, ],
                           sd = exp(as.numeric(param_Y1_sigma)))
        lik.pt1.2 <- stats::pnorm(Y2[i],
                           mean = as.numeric(param_Y2_predictors) %*% predictors_Y2[i, ],
                           sd = exp(as.numeric(param_Y2_sigma)))
        lik.pt2 <- stats::dnorm(Y1[i],
                         mean = as.numeric(param_Y1_predictors) %*% predictors_Y1[i, ],
                         sd = exp(as.numeric(param_Y1_sigma)))
        lik.pt3 <- stats::dnorm(Y2[i],
                         mean = as.numeric(param_Y2_predictors) %*% predictors_Y2[i, ],
                         sd = exp(as.numeric(param_Y2_sigma)))
        if (copula == "Clayton") {
            fun <- ((lik.pt1.1)^(-phi) - 1) + ((lik.pt1.2)^(-phi) - 1)
            joints <- (fun + 1)^(-1/phi)
            d1.joints <- (fun + 1)^(-1/phi - 1) * (lik.pt1.1)^(-phi - 1) * (-lik.pt2)
            d2.joints <- (fun + 1)^(-1/phi - 1) * (lik.pt1.2)^(-phi - 1) * (-lik.pt3)
            dd.joints <- d1.joints * d2.joints * fun^(-1) * (fun + 1)^(1/phi + 1)
            lik <- dd.joints * (fun * (1 + phi) * (fun + 1)^(-1))
        }
        if (copula == "2param") {
            fun <- ((lik.pt1.1)^(-phi) - 1)^theta + ((lik.pt1.2)^(-phi) - 1)^theta
            joints <- (fun^(1/theta) + 1)^(-1/phi)
            d1.joints <- (fun^(1/theta) + 1)^(-1/phi - 1) * fun^(1/theta - 1) *
                         ((lik.pt1.1)^(-phi) - 1)^(theta - 1) *
                         (lik.pt1.1)^(-phi - 1) * (-lik.pt2)
            d2.joints <- (fun^(1/theta) + 1)^(-1/phi - 1) * fun^(1/theta - 1) *
                         ((lik.pt1.2)^(-phi) - 1)^(theta - 1) *
                         (lik.pt1.2)^(-phi - 1) * (-lik.pt3)
            dd.joints <- d1.joints * d2.joints * fun^(-1/theta) *
                         (fun^(1/theta) + 1)^(1/phi + 1)
            lik <- dd.joints * ((theta - 1) * phi + fun^(1/theta) *
                   (1 + phi) * (fun^(1/theta) + 1)^(-1))
        }
        minusloglik <- minusloglik - log(lik)
    }
    return(as.numeric(minusloglik))
}
