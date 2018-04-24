#'  Phenotypic variance explained by genetic variants.
#'
#' Function to estimate the percentage of the variance of a phenotype
#' that can be explained by given single nucleotide variants
#' (SNVs).
#'
#' Four different approaches are available to estimate the percentage of
#' explained phenotypic variance (Laird & Lange, 2011):
#'
#' (1) \code{"Rsquared_unadj"}: Unadjusted \eqn{R^2}{R2} from a linear regression of
#'     the phenotype conditional on all provided SNVs.
#'
#' (2) \code{"Rsquared_adj"}: Adjusted \eqn{R^2}{R2} from a linear regression of
#'     the phenotype conditional on all provided SNVs.
#'
#' (3) \code{"MAF_based"}: Expected explained phenotypic variance computed based on the
#'     MAF and effect size of the provided causal SNVs.
#'
#' (4) \code{"MAF_based_Y_adjusted"}: Expected explained phenotypic variance computed
#'     based on the MAF and effect size of the causal SNVs, with respect to the empirical
#'     phenotypic variance, which is the broad-sense heritability relative to the
#'     empirical phenotypic variance.
#'
#' References:
#'
#' Laird NM, Lange C (2011). The fundamentals of modern statistical genetics. New York: Springer.
#'
#' @param genodata Numeric vector or dataframe containing the genetic variant(s) in
#'                 columns. Must be in allelic coding 0, 1, 2.
#' @param phenodata Numeric vector or dataframe of the phenotype.
#' @param type String (vector) specifying the estimation approach(es) that are computed.
#'             Available are the methods \code{"Rsquared_unadj"},
#'             \code{"Rsquared_adj"}, \code{"MAF_based"}, and
#'             \code{"MAF_based_Y_adjusted"}. See below for more details.
#' @param causal_idx Vector with entries \code{TRUE}, \code{FALSE} specifying which
#'                   SNVs are causal. Has to be supplied for the approaches
#'                   \code{MAF_based} and \code{MAF_based_Y_adjusted}.
#' @param effect_causal Numeric vector containing the effect sizes of the causal SNVs.
#'                      Has to be supplied for the approaches \code{MAF_based} and
#'                      \code{MAF_based_Y_adjusted}.
#' @return A list containing the estimated percentage of explained phenotypic variance.
#'
#' @examples
#'
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 1000)
#' phenodata <- generate_phenodata_1_simple(genodata = genodata[,1],
#'                                          type = "quantitative", b = 0)
#' compute_expl_var(genodata = genodata, phenodata = phenodata$Y,
#'                  type = c("Rsquared_unadj", "Rsquared_adj"),
#'                  causal_idx = NULL, effect_causal = NULL)
#'
#' @export
#'

compute_expl_var <- function(genodata = NULL, phenodata = NULL,
                             type = "Rsquared_unadj", causal_idx = NULL,
                             effect_causal = NULL) {
    ExplVar <- NULL
    if (is.null(phenodata) | is.null(genodata)) {
        stop("Genodata and phenodata have to be supplied.")
    }
    if (!dim(genodata)[1] == length(phenodata)) {
        stop("Number of observations in genodata and phenodata has to be the same.")
    }
    if (class(genodata) == "data.frame") {
        if (!all(sapply(genodata, is.numeric))) {
            genodata <- as.data.frame(data.matrix(genodata))
        }
    }
    if (is.null(dim(genodata))) {
        if (!is.numeric(genodata)) {
            stop("genodata has to be numeric.")
        }
    }
    if (class(phenodata) == "data.frame") {
        if (!all(sapply(phenodata, is.numeric))) {
            phenodata <- as.data.frame(data.matrix(phenodata))
            warning("phenodata contains non-numeric Variables,
                     which have automatically been transformed.")
        }
    }
    if (is.null(dim(phenodata))) {
        if (!is.numeric(phenodata)) {
            stop("phenodata has to be numeric.")
        }
    }
    if (is.null(causal_idx)) {
        warning("A vector indicating the causal SNVs is not supplied \n
                so the estimate of explained variance is based on all SNVs.")
    }
    if (any(type %in% c("MAF_based", "MAF_based_Y_adjusted")) & is.null(effect_causal)) {
        stop("A vector indicating the genetic effect sizes of the causal SNVs has \n
             to be supplied for the MAF_based and MAF_based_Y_adjusted approach.")
    }
    if (!is.null(causal_idx)) {
        genodata <- as.data.frame(genodata[, causal_idx])
    }
    MAF_causal <- compute_MAF(genodata)
    phenodata <- as.data.frame(phenodata)
    dat <- as.data.frame(append(phenodata, genodata))
    dat <- dat[stats::complete.cases(dat),]
    names(dat)[1] <- "Y"
    ExplVar <- list()
    if ("Rsquared_unadj" %in% type) {
        ExplVar$Rsquared_unadj <- summary(stats::lm(Y ~ ., data = dat))$r.squared
    }
    if ("Rsquared_adj" %in% type) {
        ExplVar$Rsquared_adj <- summary(stats::lm(Y ~ ., data = dat))$adj.r.squared
    }
    if ("MAF_based" %in% type) {
        ExplVar$MAF_based <- sum(2 * MAF_causal * (1 - MAF_causal) * effect_causal^2)
    }
    if ("MAF_based_Y_adjusted" %in% type) {
        ExplVar$MAF_based_Y_adjusted <- sum(2 * MAF_causal *
                                            (1 - MAF_causal) *
                                            effect_causal^2)/stats::var(dat$Y)
    }
    return(ExplVar)
}
