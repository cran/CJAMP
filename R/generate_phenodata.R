#' Functions to generate phenotype data.
#'
#' Functions to generate standard normal or binary phenotypes based on provided genetic
#' data, for specified effect sizes.
#' The functions \code{\link{generate_phenodata_1_simple}} and
#' \code{\link{generate_phenodata_1}} generate one phenotype Y conditional on
#' single nucleotide variants (SNVs) and two covariates.
#' \code{\link{generate_phenodata_2_bvn}} as well as \code{\link{generate_phenodata_2_copula}}
#' generate two phenotypes \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2} with dependence Kendall's tau conditional on
#' the provided SNVs and two covariates.
#'
#' In more detail, the function \code{\link{generate_phenodata_1_simple}}
#' generates a quantitative or binary phenotype Y with n observations,
#' conditional on the specified SNVs with given effect sizes and conditional
#' on one binary and one standard normally-distributed covariate with
#' specified effect sizes. n is given through the provided SNVs.
#'
#' \code{\link{generate_phenodata_1}} provides an extension of
#' \code{\link{generate_phenodata_1_simple}} and allows to further select
#' the percentage of causal SNVs, a minor allele frequency cutoff on the
#' causal SNVs, and varying effect directions. n is given through the
#' provided SNVs.
#'
#' The function \code{\link{generate_phenodata_2_bvn}} generates
#' two quantitative phenotypes \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2} conditional on one binary and one
#' standard normally-distributed covariate \eqn{X_1}{X1}, \eqn{X_2}{X2} from the bivariate
#' normal distribution so that they have have dependence \eqn{\tau} given
#' by Kendall's \code{tau}.
#'
#' The function \code{\link{generate_phenodata_2_copula}} generates
#' two quantitative phenotypes \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2} conditional on one binary and one
#' standard normally-distributed covariate \eqn{X_1}{X1}, \eqn{X_2}{X2} from the Clayton copula
#' so that \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2} are marginally normally distributed and have dependence
#' Kendall's tau specified by \code{tau} or \code{phi}, using the function
#' \code{\link{generate_clayton_copula}}.
#'
#' The genetic effect sizes are the specified numeric values \code{b} and
#' \code{b1, b2}, respectively, in the functions \code{\link{generate_phenodata_1_simple}}
#' and \code{\link{generate_phenodata_2_bvn}}. In
#' \code{\link{generate_phenodata_1}} and \code{\link{generate_phenodata_2_copula}},
#' the genetic effect sizes are computed by multiplying \code{b} or \code{b1, b2},
#' respectively, with the absolute value of the log10-transformed
#' minor allele frequencies, so that rarer variants have larger effect sizes.
#'
#' @param genodata Numeric input vector or dataframe containing the genetic
#'                 variant(s) in columns. Must be in allelic coding 0, 1, 2.
#' @param type String with value \code{"quantitative"} or \code{"binary"}
#'             specifying whether normally-distributed or binary phenotypes
#'             are generated.
#' @param a Numeric vector specifying the effect sizes of the covariates \eqn{X_1}{X1}, \eqn{X_2}{X2}
#'          in the data generation.
#' @param a1 Numeric vector specifying the effect sizes of the covariates \eqn{X_1}{X1}, \eqn{X_2}{X2}
#'           on the first phenotype in the data generation.
#' @param a2 Numeric vector specifying the effect sizes of the covariates \eqn{X_1}{X1}, \eqn{X_2}{X2}
#'           on the second phenotype in the data generation.
#' @param b Integer or vector specifying the genetic effect size(s) of
#'          the provided SNVs (\code{genodata}) in the data generation.
#' @param b1 Integer or vector specifying the genetic effect size(s) of
#'           the provided SNVs (\code{genodata}) on the first phenotype
#'           in the data generation.
#' @param b2 Integer or vector specifying the genetic effect size(s) of
#'           the provided SNVs (\code{genodata}) on the second phenotype
#'           in the data generation.
#' @param phi Integer specifying the parameter \eqn{\phi} for
#'            the dependence between the two generated phenotypes.
#' @param tau Integer specifying Kendall's tau, which determines the
#'            dependence between the two generated phenotypes.
#' @param MAF_cutoff Integer specifying a minor allele frequency cutoff to
#'                   determine among which SNVs the causal SNVs are
#'                   sampled for the phenotype generation.
#' @param prop_causal Integer specifying the desired percentage of causal SNVs
#'                    among all SNVs.
#' @param direction String with value \code{"a"}, \code{"b"}, or \code{"c"}
#'                  specifying whether all causal SNVs have a positive effect on
#'                  the phenotypes (\code{"a"}), 20\% of the causal SNVs have a
#'                  negative effect and 80\% a positive effect on the phenotypes
#'                  (\code{"b"}), or 50\% of the causal SNVs have a negative
#'                  effect and 50\% a positive effect on the phenotypes (\code{"c"}).
#' @return A dataframe containing n observations of the phenotype Y or phenotypes
#'         \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2} and of the covariates \eqn{X_1}{X1}, \eqn{X_2}{X2}.
#'
#' @examples
#'
#' # Generate genetic data:
#' set.seed(10)
#' genodata <- generate_genodata(n_SNV = 20, n_ind = 1000)
#' compute_MAF(genodata)
#'
#' # Generate different phenotype data:
#' phenodata1 <- generate_phenodata_1_simple(genodata = genodata[,1],
#'                                           type = "quantitative", b = 0)
#' phenodata2 <- generate_phenodata_1_simple(genodata = genodata[,1],
#'                                           type = "quantitative", b = 2)
#' phenodata3 <- generate_phenodata_1_simple(genodata = genodata,
#'                                           type = "quantitative", b = 2)
#' phenodata4 <- generate_phenodata_1_simple(genodata = genodata,
#'                                           type = "quantitative",
#'                                           b = seq(0.1, 2, 0.1))
#' phenodata5 <- generate_phenodata_1_simple(genodata = genodata[,1],
#'                                           type = "binary", b = 0)
#' phenodata6 <- generate_phenodata_1(genodata = genodata[,1],
#'                                    type = "quantitative", b = 0,
#'                                    MAF_cutoff = 1, prop_causal = 0.1,
#'                                    direction = "a")
#' phenodata7 <- generate_phenodata_1(genodata = genodata,
#'                                    type = "quantitative", b = 0.6,
#'                                    MAF_cutoff = 0.1, prop_causal = 0.05,
#'                                    direction = "a")
#' phenodata8 <- generate_phenodata_1(genodata = genodata,
#'                                    type = "quantitative",
#'                                    b = seq(0.1, 2, 0.1),
#'                                    MAF_cutoff = 0.1, prop_causal = 0.05,
#'                                    direction = "a")
#' phenodata9 <- generate_phenodata_2_bvn(genodata = genodata[,1],
#'                                        tau = 0.5, b1 = 0, b2 = 0)
#' phenodata10 <- generate_phenodata_2_bvn(genodata = genodata,
#'                                         tau = 0.5, b1 = 0, b2 = 0)
#' phenodata11 <- generate_phenodata_2_bvn(genodata = genodata,
#'                                         tau = 0.5, b1 = 1,
#'                                         b2 = seq(0.1,2,0.1))
#' phenodata12 <- generate_phenodata_2_bvn(genodata = genodata,
#'                                         tau = 0.5, b1 = 1, b2 = 2)
#' par(mfrow = c(3, 1))
#' hist(phenodata12$Y1)
#' hist(phenodata12$Y2)
#' plot(phenodata12$Y1, phenodata12$Y2)
#'
#' phenodata13 <- generate_phenodata_2_copula(genodata = genodata[,1],
#'                                            MAF_cutoff = 1, prop_causal = 1,
#'                                            tau = 0.5, b1 = 0, b2 = 0)
#' phenodata14 <- generate_phenodata_2_copula(genodata = genodata,
#'                                            MAF_cutoff = 1, prop_causal = 0.5,
#'                                            tau = 0.5, b1 = 0, b2 = 0)
#' phenodata15 <- generate_phenodata_2_copula(genodata = genodata,
#'                                            MAF_cutoff = 1, prop_causal = 0.5,
#'                                            tau = 0.5, b1 = 0, b2 = 0)
#' phenodata16 <- generate_phenodata_2_copula(genodata = genodata,
#'                                            MAF_cutoff = 1, prop_causal = 0.5,
#'                                            tau = 0.2, b1 = 0.3,
#'                                            b2 = seq(0.1, 2, 0.1))
#' phenodata17 <- generate_phenodata_2_copula(genodata = genodata,
#'                                            MAF_cutoff = 1, prop_causal = 0.5,
#'                                            tau = 0.2, b1 = 0.3, b2 = 0.3)
#' par(mfrow = c(3, 1))
#' hist(phenodata17$Y1)
#' hist(phenodata17$Y2)
#' plot(phenodata17$Y1, phenodata17$Y2)
#'
#'
#'
#'
#'
#'

#' @name generate_phenodata
NULL

#' @rdname generate_phenodata
#' @export

generate_phenodata_1_simple <- function(genodata = NULL, type = "quantitative",
                                        b = 0, a = c(0, 0.5, 0.5)) {
    if (is.null(genodata)) {
        stop("Genotype data has to be supplied.")
    }
    genodata <- as.matrix(genodata)
    genodata <- genodata[stats::complete.cases(genodata),]
    genodata <- as.matrix(genodata)
    n_ind <- dim(genodata)[1]
    n_SNV <- dim(genodata)[2]
    if ((!n_SNV == length(b)) & (!n_SNV == length(b)) & (!length(b) == 1)) {
        stop("Genetic effects b has to be of same length as \n
             the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b)) & (length(b) == 1)) {
        b <- rep(b, n_SNV)
    }
    X1 <- stats::rnorm(n_ind, mean = 0, sd = 1)
    X2 <- stats::rbinom(n_ind, size = 1, prob = 0.5)
    if (type == "quantitative") {
        Y <- a[1] + a[2] * X1 + a[3] * X2 + genodata %*% b + stats::rnorm(n_ind, 0, 1)
    }
    if (type == "binary") {
        P <- 1/(1 + exp(-(a[1] + a[2] * X1 + a[3] * X2 + b * genodata)))
        Y <- stats::rbinom(n_ind, 1, P)
    }
    phenodata <- data.frame(Y = Y, X1 = X1, X2 = X2)
    return(phenodata)
}

#' @rdname generate_phenodata
#' @export

generate_phenodata_1 <- function(genodata = NULL, type = "quantitative", b = 0.6,
                                 a = c(0, 0.5, 0.5), MAF_cutoff = 1,
                                 prop_causal = 0.1, direction = "a") {
    if (is.null(genodata)) {
        stop("Genotype data has to be supplied")
    }
    genodata <- as.data.frame(genodata)
    genodata <- genodata[stats::complete.cases(genodata),]
    genodata <- as.data.frame(genodata)
    n_ind <- dim(genodata)[1]
    n_SNV <- dim(genodata)[2]
    if ((!n_SNV == length(b)) & (!length(b) == 1)) {
        stop("Genetic effects b has to be of same length as \n
             the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b)) & (length(b) == 1)) {
        b <- rep(b, n_SNV)
    }
    X1 <- stats::rnorm(n_ind, mean = 0, sd = 1)
    X2 <- stats::rbinom(n_ind, size = 1, prob = 0.5)
    sigma1 <- 1
    causal_idx <- FALSE
    help_causal_idx_counter <- 0
    while (!any(causal_idx)) {
      # if no causal variants are selected (MAF <= MAF_cutoff), do it again up to 10 times
      help_causal_idx_counter <- help_causal_idx_counter + 1
      if(help_causal_idx_counter == 10){ stop("MAF cutoff too low, no causal SNVs") }
      help_idx <- sample(1:dim(genodata)[2],
                           ceiling(prop_causal * dim(genodata)[2]))
        causal_idx <- ((1:dim(genodata)[2] %in% help_idx) &
                       (compute_MAF(genodata) < MAF_cutoff))
    }
    geno_causal <- as.matrix(genodata[, causal_idx])
    b <- b[causal_idx]
    alpha_g <- b * abs(log10(compute_MAF(geno_causal)))
    # if direction b or c is specified, change direction of effect for some variants
    if (direction == "b") {
        help_idx_2 <- sample(1:dim(geno_causal)[2],
                             round(0.2 * dim(geno_causal)[2]))
        alpha_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in%help_idx_2)) + 1) * alpha_g
    }
    if (direction == "c") {
        help_idx_3 <- sample(1:dim(geno_causal)[2],
                             round(0.5 * dim(geno_causal)[2]))
        alpha_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in% help_idx_3)) + 1) * alpha_g
    }
    epsilon1 <- stats::rnorm(n_ind, mean = 0, sd = sigma1)
    if (type == "quantitative") {
        Y <- a[1] + a[2] * X1 + a[3] * X2 + geno_causal %*% alpha_g + epsilon1
        Y <- as.numeric(Y)
    }
    if (type == "binary") {
        P <- 1/(1 + exp(-(a[1] + a[2] * scale(X1, center = T, scale = F) +
                          a[3] * scale(X2, center = T, scale = F) +
                          geno_causal %*% alpha_g)))
        Y <- stats::rbinom(n_ind, 1, P)
        Y <- as.numeric(Y)
    }
    phenodata <- data.frame(Y = Y, X1 = X1, X2 = X2)
    return(phenodata)
}

#' @rdname generate_phenodata
#' @export

generate_phenodata_2_bvn <- function(genodata = NULL, tau = NULL, b1 = 0,
                                     b2 = 0, a1 = c(0, 0.5, 0.5),
                                     a2 = c(0, 0.5, 0.5)) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if (is.null(genodata)) {
        stop("Genotype data has to be supplied.")
    }
    genodata <- as.matrix(genodata)
    genodata <- genodata[stats::complete.cases(genodata),]
    genodata <- as.matrix(genodata)
    n_ind <- dim(genodata)[1]
    n_SNV <- dim(genodata)[2]
    if ((!n_SNV == length(b1)) & (!n_SNV == length(b1)) & (!length(b1) == 1)) {
        stop("Genetic effects b1 has to be of same length \n
             as the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b2)) & (!n_SNV == length(b2)) & (!length(b2) == 1)) {
        stop("Genetic effects b2 has to be of same length \n
             as the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b1)) & (length(b1) == 1)) {
        b1 <- rep(b1, n_SNV)
    }
    if ((!n_SNV == length(b2)) & (length(b2) == 1)) {
        b2 <- rep(b2, n_SNV)
    }
    # generate the two covariates
    X1 <- stats::rnorm(n_ind, mean = 0, sd = 1)
    X2 <- stats::rbinom(n_ind, size = 1, prob = 0.5)
    # set the weight of the nongenetic covariates set the residual variances
    sigma1 <- 1
    sigma2 <- 1
    if (is.null(tau)) {
        stop("Tau has to be provided.")
    }
    rho <- sin(tau * pi/2)
    mu <- c(0, 0)
    sigma <- matrix(c(sigma1, rho, rho, sigma2), 2, 2)
    epsilon <- MASS::mvrnorm(n = n_ind, mu = mu, Sigma = sigma)
    epsilon1 <- epsilon[, 1]
    epsilon2 <- epsilon[, 2]
    Y1 <- a1[1] + a1[2] * X1 + a1[3] * X2 + genodata %*% b1 + epsilon1
    Y2 <- a2[1] + a2[2] * X1 + a2[3] * X2 + genodata %*% b2 + epsilon2
    phenodata <- data.frame(Y1 = Y1, Y2 = Y2, X1 = X1, X2 = X2)
    return(phenodata)
}

#' @rdname generate_phenodata
#' @export

generate_phenodata_2_copula <- function(genodata = NULL, phi = NULL, tau = 0.5,
                                        b1 = 0.6, b2 = 0.6, a1 = c(0, 0.5, 0.5),
                                        a2 = c(0, 0.5, 0.5), MAF_cutoff = 1,
                                        prop_causal = 0.1, direction = "a") {
    if (is.null(genodata)) {
        stop("Genotype data has to be supplied.")
    }
    genodata <- as.matrix(genodata)
    genodata <- genodata[stats::complete.cases(genodata),]
    genodata <- as.matrix(genodata)
    n_ind <- dim(genodata)[1]
    n_SNV <- dim(genodata)[2]
    if ((!n_SNV == length(b1)) & (!n_SNV == length(b1)) & (!length(b1) == 1)) {
        stop("Genetic effects b1 has to be of same length \n
             as the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b2)) & (!n_SNV == length(b2)) & (!length(b2) == 1)) {
        stop("Genetic effects b2 has to be of same length \n
             as the number of provided SNVs or of length 1.")
    }
    if ((!n_SNV == length(b1)) & (length(b1) == 1)) {
        b1 <- rep(b1, n_SNV)
    }
    if ((!n_SNV == length(b2)) & (length(b2) == 1)) {
        b2 <- rep(b2, n_SNV)
    }
    X1 <- stats::rnorm(n_ind, mean = 0, sd = 1)
    X2 <- stats::rbinom(n_ind, size = 1, prob = 0.5)
    sigma1 <- 1
    sigma2 <- 1
    causal_idx <- FALSE
    help_causal_idx_counter <- 0
    while (!any(causal_idx)) {
        # if no causal variants are selected (MAF <= MAF_cutoff), do it again up to 10 times
        help_causal_idx_counter <- help_causal_idx_counter + 1
        if(help_causal_idx_counter == 10){ stop("MAF cutoff too low, no causal SNVs") }
        help_idx <- sample(1:dim(genodata)[2],
                           ceiling(prop_causal * dim(genodata)[2]))
        causal_idx <- ((1:dim(genodata)[2] %in% help_idx) &
                       (compute_MAF(genodata) < MAF_cutoff))
    }
    geno_causal <- as.matrix(genodata[, causal_idx])
    b1 <- b1[causal_idx]
    b2 <- b2[causal_idx]
    alpha_g <- b1 * abs(log10(compute_MAF(geno_causal)))
    beta_g <- b2 * abs(log10(compute_MAF(geno_causal)))
    # if b or c is specified, change effect direction for some variants
    if (direction == "b") {
        help_idx_2 <- sample(1:dim(geno_causal)[2],
                             round(0.2 * dim(geno_causal)[2]))
        alpha_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in% help_idx_2)) + 1) * alpha_g
        beta_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in% help_idx_2)) + 1) * beta_g
    }
    if (direction == "c") {
        help_idx_3 <- sample(1:dim(geno_causal)[2],
                             round(0.5 * dim(geno_causal)[2]))
        alpha_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in% help_idx_3)) + 1) * alpha_g
        beta_g <- (-1)^(as.numeric(!(1:dim(geno_causal)[2] %in% help_idx_3)) + 1) * beta_g
    }
    # compute copula parameter between the two phenotypes from tau if tau is given
    if (!is.null(tau)) {
        phi <- (2 * tau)/(1 - tau)
    }
    res_copula <- generate_clayton_copula(n = n_ind, phi = phi)
    epsilon1 <- res_copula$Y1
    epsilon2 <- res_copula$Y2
    if (is.null(tau)) {
        tau <- phi/(phi + 2)
    }
    Y1 <- a1[1] + a1[2] * X1 + a1[3] * X2 + geno_causal %*% alpha_g + epsilon1
    Y1 <- as.numeric(Y1)
    Y2 <- a2[1] + a2[2] * X1 + a2[3] * X2 + geno_causal %*% beta_g + epsilon2
    Y2 <- as.numeric(Y2)
    phenodata <- data.frame(Y1 = Y1, Y2 = Y2, X1 = X1, X2 = X2)
    print(paste("Kendall's tau between Y1, Y2 = ", tau, sep=""))
    return(phenodata)
}
