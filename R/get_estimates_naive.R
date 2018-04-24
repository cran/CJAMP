#' Naive estimates of the copula and marginal parameters.
#'
#' Function to compute naive estimates of the copula parameter(s)
#' and maximum likelihood (ML) estimates of the marginal parameters in a joint
#' copula model of \code{Y1} and \code{Y2} given the predictors of \code{Y1}
#' and \code{Y2}. The main use of the function is to provide parameter
#' starting values for the optimization of the log-likelihood function of the
#' joint copula model in \code{\link{cjamp}} in order to obtain maximum
#' likelihood estimates in the copula model.
#'
#' The estimates of the copula parameter(s) include estimates of \eqn{\phi}
#' (if \code{copula_param == "phi"}), \eqn{\theta} (if
#' \code{copula_param == "theta"}) or both (if \code{copula_param == "both"}).
#' They are obtained by computing Kendall's tau between \code{Y1} and \code{Y2}
#' and using the relationship \eqn{\tau = \phi/(\phi+2)} of the Clayton
#' copula to obtain an estimate of \eqn{\phi} and \eqn{\tau = (\theta-1)/\theta}
#' of the Gumbel copula to obtain an estimate of \eqn{\theta}.
#'
#' The ML estimates of the marginal parameters include estimates of the log standard
#' deviations of \code{Y1}, \code{Y2} given their predictors (\eqn{log(\sigma1), log(\sigma2)})
#' and of the effects of \code{predictors_Y1} on \code{Y1} and
#' \code{predictors_Y2} on \code{Y2}. The estimates of the marginal effects are
#' obtained from linear regression models of \code{Y1} given \code{predictors_Y1}
#' and \code{Y2} given \code{predictors_Y2}, respectively. If single nucleotide
#' variants (SNVs) are included as predictors, the genetic effect estimates
#' are obtained from an underlying additive genetic model if SNVs are provided
#' as 0-1-2 genotypes and from an underlying dominant model if SNVs are provided
#' as 0-1 genotypes.
#'
#' @param Y1 Numeric vector containing the first phenotype.
#' @param Y2 Numeric vector containing the second phenotype.
#' @param predictors_Y1 Dataframe containing the predictors of \code{Y1}
#'                      in columns.
#' @param predictors_Y2 Dataframe containing the predictors of \code{Y2}
#'                      in columns.
#' @param copula_param String indicating whether estimates should be computed
#'                     for \eqn{\phi} (\code{"phi"}), for \eqn{\theta}
#'                     (\code{"theta"}), or both (\code{"both"}).
#' @return Vector of the numeric estimates of the copula parameters
#'         \eqn{log(\phi)} and/or \eqn{log(\theta-1)}, of the marginal
#'         parameters (\eqn{log(\sigma1), log(\sigma2)}, and estimates
#'         of the effects of the predictors \code{predictors_Y1}
#'         on \code{Y1} and \code{predictors_Y2} on \code{Y2}).
#'
#' @examples
#' # Generate genetic data:
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 1000)
#'
#' # Generate phenotype data:
#' phenodata <- generate_phenodata_2_copula(genodata = genodata, MAF_cutoff = 1,
#'                                          prop_causal = 0.5, tau = 0.2,
#'                                          b1 = 0.3, b2 = 0.3)
#' predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#'                          SNV = genodata$SNV1)
#'
#' get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                     predictors_Y1 = predictors, predictors_Y2 = predictors,
#'                     copula_param = "both")
#'
#' @export
#'

get_estimates_naive <- function(Y1 = NULL, Y2 = NULL, predictors_Y1 = NULL,
                                predictors_Y2 = NULL, copula_param = "both") {
    if (is.null(Y1) | is.null(Y2)) {
        stop("Both Y1 and Y2 have to be supplied.")
    }
    dataframe_Y1 <- data.frame(Y1 = Y1)
    dataframe_Y2 <- data.frame(Y2 = Y2)
    n_ind <- dim(dataframe_Y1)[1]
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
    if (!is.null(predictors_Y1)) {
        predictors_Y1 <- data.frame(predictors_Y1)
    }
    if (!is.null(predictors_Y2)) {
        predictors_Y2 <- data.frame(predictors_Y2)
    }
    if (!n_ind == dim(dataframe_Y2)[1]) {
        stop("Variables must have same length.")
    }
    if (!is.null(predictors_Y1)) {
        if (!n_ind == dim(predictors_Y1)[1]) {
            stop("Variables must have same length.")
        }
    }
    if (!is.null(predictors_Y2)) {
        if (!n_ind == dim(predictors_Y2)[1]) {
            stop("Variables must have same length.")
        }
    }
    if (!is.null(predictors_Y1)) {
        dataframe_Y1 <- cbind(dataframe_Y1, predictors_Y1)
    }
    if (!is.null(predictors_Y2)) {
        dataframe_Y2 <- cbind(dataframe_Y2, predictors_Y2)
    }
    res_Y1 <- stats::lm(Y1 ~ ., data = dataframe_Y1)
    res_Y2 <- stats::lm(Y2 ~ ., data = dataframe_Y2)
    param_Y1_sigma <- log(stats::sd(res_Y1$residuals, na.rm = T))
    names(param_Y1_sigma) <- "Y1_log_sigma"
    param_Y2_sigma <- log(stats::sd(res_Y2$residuals, na.rm = T))
    names(param_Y2_sigma) <- "Y2_log_sigma"
    tau <- stats::cor(Y1, Y2, use = "complete.obs", method = c("kendall"))
    if (tau == 0) {
        tau <- 1e-04
    } # to avoid technical breakdown of function
    log_phi <- log(2 * tau/(1 - tau))
    log_theta <- log((1/(1 - tau)) - 1)
    param_Y1_Y2 <- c(log_phi, log_theta)
    names(param_Y1_Y2) <- c("log_phi", "log_theta_minus1")
    if (copula_param == "phi") {
        param_Y1_Y2 <- param_Y1_Y2[1]
    }
    if (copula_param == "theta") {
        param_Y1_Y2 <- param_Y1_Y2[2]
    }
    param_Y1_predictors <- res_Y1$coefficients
    names(param_Y1_predictors) <- paste("Y1_", names(param_Y1_predictors), sep = "")
    param_Y2_predictors <- res_Y2$coefficients
    names(param_Y2_predictors) <- paste("Y2_", names(param_Y2_predictors), sep = "")
    if (any(is.na(param_Y1_Y2))) {
        warning("One or both dependence parameters could not be estimated.")
    }
    if (any(is.na(param_Y1_sigma)) | any(is.na(param_Y2_sigma))) {
        warning("One or both marginal variances could not be estimated.")
    }
    if (any(is.na(param_Y1_predictors)) | any(is.na(param_Y2_predictors))) {
        warning("One or more marginal parameters could not be estimated.")
    }
    estimates <- c(param_Y1_Y2, param_Y1_sigma, param_Y2_sigma,
                   param_Y1_predictors, param_Y2_predictors)
    return(estimates)
}
