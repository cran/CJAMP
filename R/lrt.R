#' Compute likelihood ratio tests.
#'
#' Functions to compute likelihood ratio tests of two nested copula models
#' with the same marginal models, and of marginal parameters in the same copula
#' model but with different nested marginal models where one or more parameters
#' are set to 0 (Yilmaz & Lawless, 2011).
#'
#' References:
#'
#' Yilmaz YE, Lawless JF (2011). Likelihood ratio procedures and tests of fit in parametric and semiparametric copula models with censored data. Lifetime Data Analysis, 17(3): 386-408.
#'
#' @param minlogl_null Minus log-likelihood of the null model, which is nested
#'                     within the larger alternative model \code{minlogl_altern}.
#' @param minlogl_altern Minus log-likelihood of the alternative model, which
#'                       contains the null model \code{minlogl_null}.
#' @param df Degrees of freedom in the likelihood ratio test of marginal
#'              parameters in the same copula model, i.e. the number of
#'              parameters that are in the alternative model but that are not
#'              contained in the smaller nested null model.
#' @return List including the \eqn{\chi^2} value and p-value of the likelihood
#'         ratio test.
#'
#' @examples
#'
#' # Example 1: Test whether 2-parameter copula model has a better
#' #            model fit compared to Clayton copula (no).
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 100)
#' phenodata <- generate_phenodata_2_copula(genodata = genodata$SNV1,
#'                                          MAF_cutoff = 1, prop_causal = 1,
#'                                          tau = 0.2, b1 = 0.3, b2 = 0.3)
#' predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#'                          SNV = genodata$SNV1)
#' estimates_c <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                    predictors_Y1 = predictors,
#'                                    predictors_Y2 = predictors,
#'                                    copula_param = "phi")
#' minusloglik_Clayton <- minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                    predictors_Y1 = predictors,
#'                                    predictors_Y2 = predictors,
#'                                    parameters = estimates_c, copula = "Clayton")
#' estimates_2p <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                     predictors_Y1 = predictors,
#'                                     predictors_Y2 = predictors,
#'                                     copula_param = "both")
#' minusloglik_2param <- minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                   predictors_Y1 = predictors,
#'                                   predictors_Y2 = predictors,
#'                                   parameters = estimates_2p, copula = "2param")
#' lrt_copula(minusloglik_Clayton, minusloglik_2param)
#'
#' # Example 2: Test marginal parameters (alternative model has better fit).
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 100)
#' phenodata <- generate_phenodata_2_copula(genodata = genodata$SNV1,
#'                                          MAF_cutoff = 1, prop_causal = 1,
#'                                          tau = 0.2, b1 = 2, b2 = 2)
#' predictors_1 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' estimates_1 <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                    predictors_Y1 = predictors_1,
#'                                    predictors_Y2 = predictors_1,
#'                                    copula_param = "phi")
#' minusloglik_1 <- minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                              predictors_Y1 = predictors_1,
#'                              predictors_Y2 = predictors_1,
#'                              parameters = estimates_1, copula = "Clayton")
#' predictors_2 <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#'                            SNV = genodata$SNV1)
#' estimates_2 <- get_estimates_naive(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                                    predictors_Y1 = predictors_2,
#'                                    predictors_Y2 = predictors_2,
#'                                    copula_param = "phi")
#' minusloglik_2 <- minusloglik(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#'                              predictors_Y1 = predictors_2,
#'                              predictors_Y2 = predictors_2,
#'                              parameters = estimates_2, copula = "Clayton")
#' lrt_param(minusloglik_1, minusloglik_2, df=2)
#'

#' @name lrt
NULL

#' @rdname lrt
#' @export

lrt_copula <- function(minlogl_null = NULL, minlogl_altern = NULL) {
    if (is.null(minlogl_null) | is.null(minlogl_altern)) {
        stop("minlogl_null and minlogl_altern have to be supplied.")
    }
    chisq <- 2 * minlogl_null - 2 * minlogl_altern
    pval <- 1 - (0.5 + 0.5 * stats::pchisq(chisq, df = 1))
    return(list(chisq = chisq, pval = pval))
}

#' @rdname lrt
#' @export

lrt_param <- function(minlogl_null = NULL, minlogl_altern = NULL, df = NULL) {
    if (is.null(minlogl_null) | is.null(minlogl_altern) | is.null(df)) {
        stop("minlogl_null, minlogl_altern and df have to be supplied.")
    }
    chisq <- 2 * minlogl_null - 2 * minlogl_altern
    pval <- 1 - stats::pchisq(chisq, df = df)
    return(list(chisq = chisq, pval = pval))
}

