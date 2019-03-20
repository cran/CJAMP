#' Functions to generate genetic data.
#'
#' Functions to generate genetic data in the form of single
#' nucleotide variants (SNVs). The function
#' \code{\link{generate_singleton_data}} generates singletons (i.e. SNVs with
#' one observed minor allele); \code{\link{generate_doubleton_data}}
#' generates doubletons (i.e. SNVs with two observed minor alleles), and the
#' function \code{\link{generate_genodata}} generates \code{n_ind}
#' observations of \code{n_SNV} SNVs with random minor allele frequencies.
#'
#' @param n_ind Integer specifying the number of observations that are generated.
#' @param n_SNV Integer specifying the number of SNVs that are generated.
#' @return A dataframe containing \code{n_ind} observations of \code{n_SNV} SNVs.
#'
#' @examples
#' set.seed(10)
#' genodata1 <- generate_singleton_data()
#' compute_MAF(genodata1)
#'
#' genodata2 <- generate_doubleton_data()
#' compute_MAF(genodata2)
#'
#' genodata3 <- generate_genodata()
#' compute_MAF(genodata3)
#'
#' @export
#'

generate_genodata <- function(n_SNV = 100, n_ind = 1000) {
    genodata <- data.frame(matrix(nrow = n_ind, ncol = n_SNV))
    for (i in 1:n_SNV) {
        MAFhelp <- round(stats::runif(1, 0, 0.5), 3)
        n0 <- round((1 - MAFhelp)^2 * n_ind)
        n2 <- round(MAFhelp^2 * n_ind)
        n1 <- n_ind - n0 - n2
        helpvec <- c(rep(0, n0), rep(1, n1), rep(2, n2))
        genodata[, i] <- sample(helpvec, n_ind, replace = F)
        names(genodata)[i] <- paste("SNV", i, sep = "")
    }
    return(genodata)
}

#' @rdname generate_genodata
#' @export

generate_singleton_data <- function(n_SNV = 100, n_ind = 1000) {
    genodata <- data.frame(matrix(nrow = n_ind, ncol = n_SNV))
    helpvec <- c(1, rep(0, n_ind - 1))
    for (i in 1:n_SNV) {
        genodata[, i] <- sample(helpvec, n_ind, replace = F)
        names(genodata)[i] <- paste("SNV", i, sep = "")
    }
    return(genodata)
}

#' @rdname generate_genodata
#' @export

generate_doubleton_data <- function(n_SNV = 100, n_ind = 1000) {
    genodata <- data.frame(matrix(nrow = n_ind, ncol = n_SNV))
    helpvec <- c(1, 1, rep(0, n_ind - 2))
    for (i in 1:n_SNV) {
        genodata[, i] <- sample(helpvec, n_ind, replace = F)
        names(genodata)[i] <- paste("SNV", i, sep = "")
    }
    return(genodata)
}
