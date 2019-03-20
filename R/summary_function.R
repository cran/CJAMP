#' Summary function.
#'
#' Summary function for the \code{\link{cjamp}} and \code{\link{cjamp_loop}}
#' functions.
#'
#' @param object \code{cjamp} object (output of the \code{\link{cjamp}}
#'                or \code{\link{cjamp_loop}} function).
#' @param ... Additional arguments affecting the summary produced.
#'
#' @return Formatted data frame of the results.
#'
#' @examples
#' ## Not run. When executing, the following takes about 2 minutes running time.
#' ## Summary of regular cjamp function
#' #set.seed(10)
#' #genodata <- generate_genodata(n_SNV = 20, n_ind = 100)
#' #phenodata <- generate_phenodata_2_copula(genodata = genodata$SNV1,
#' #                                         MAF_cutoff = 1, prop_causal = 1,
#' #                                         tau = 0.2, b1 = 0.3, b2 = 0.3)
#' #predictors <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2,
#' #                        genodata[, 1:3])
#' #results <- cjamp(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#' #                 predictors_Y1 = predictors, predictors_Y2 = predictors,
#' #                 copula = "2param", optim_method = "BFGS", trace = 0,
#' #                 kkt2tol = 1E-16, SE_est = TRUE, pval_est = TRUE,
#' #                 n_iter_max = 10)
#' #summary(results)
#' #
#' ## Summary of looped cjamp function
#' #covariates <- data.frame(X1 = phenodata$X1, X2 = phenodata$X2)
#' #predictors <- genodata
#' #results <- cjamp_loop(Y1 = phenodata$Y1, Y2 = phenodata$Y2,
#' #                      covariates_Y1 = covariates,
#' #                      covariates_Y2 = covariates,
#' #                      predictors = predictors, copula = "Clayton",
#' #                      optim_method = "BFGS", trace = 0, kkt2tol = 1E-16,
#' #                      SE_est = TRUE, pval_est = TRUE, n_iter_max = 10)
#' #summary(results)
#'
#' @export
#'

summary.cjamp <- function(object = NULL, ...) {
    if (is.null(object) | !class(object) == "cjamp") {
        stop("cjamp object has to be supplied.")
    }
    # for input from cjamp function:
    if ("Parameter point estimates" %in% names(object)) {
        length_res <- length(object$"Parameter p-values")
        res_out <- data.frame(point_estimates = object$"Parameter point estimates"[3:(length_res+2)],
                              SE_estimates = object$"Parameter standard error estimates"[3:(length_res+2)],
                              pvalues = object$"Parameter p-values")
        print(paste("C-JAMP estimates of marginal parameters."))
        print(res_out)
    }
    # for input from cjamp_loop function:
    if ("Parameter p-values for effects of predictors on Y1" %in% names(object)) {
        res_out_1 <- data.frame(point_estimates = object$"Parameter point estimates for effects of predictors on Y1",
                                SE_estimates = object$"Parameter standard error estimates for effects on Y1",
                                pvalues = object$"Parameter p-values for effects of predictors on Y1")
        rownames(res_out_1) <- paste("Y1", rownames(res_out_1), sep = "_")
        print(paste("C-JAMP estimates of marginal parameters on Y1."))
        print(res_out_1)
        res_out_2 <- data.frame(point_estimates = object$"Parameter point estimates for effects of predictors on Y2",
                                SE_estimates = object$"Parameter standard error estimates for effects on Y2",
                                pvalues = object$"Parameter p-values for effects of predictors on Y2")
        rownames(res_out_2) <- paste("Y2", rownames(res_out_2), sep = "_")
        print(paste("C-JAMP estimates of marginal parameters on Y2."))
        print(res_out_2)
        res_out <- data.frame(rbind(res_out_1, res_out_2))
    }
    invisible(res_out)
}
