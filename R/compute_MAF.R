#' Compute minor allele frequency of genetic variants.
#'
#' Function to compute the minor allele frequency (MAF) of
#' one or more genetic variants.
#'
#' @param genodata Numeric vector or dataframe containing the genetic variants
#'                 in columns. Must be in allelic coding 0, 1, 2.
#' @return A vector containing the minor allele frequencies of the variants.
#'
#' @examples
#' # Example of a single variant
#' genodata <- stats::rbinom(2000, 2, 0.3)
#' compute_MAF(genodata)
#'
#' # Example of a set of variants
#' genodata <- generate_genodata()
#' compute_MAF(genodata)
#'
#' @export
#'

compute_MAF <- function(genodata) {
    MAF <- NULL
    if (is.null(dim(genodata))) {
        dat <- genodata[!is.na(genodata)]
        if (!all(unique(dat) %in% c(0, 1, 2))) {
            stop("SNP has to be supplied as genotypes 0,1,2.")
        }
        if (is.na(table(dat)["1"])) {
            genocount1 <- 0
        }
        if (!is.na(table(dat)["1"])) {
            genocount1 <- table(dat)["1"][[1]]
        }
        if (is.na(table(dat)["2"])) {
            genocount2 <- 0
        }
        if (!is.na(table(dat)["2"])) {
            genocount2 <- table(dat)["2"][[1]]
        }
        nobs <- length(dat)
        MAF <- (genocount1/nobs + 2 * genocount2/nobs)/2
        if (MAF > 0.5) {
            warning("MAF is larger than 0.5, maybe recode alleles?")
        }
    }
    if (!is.null(dim(genodata))) {
        for (i in 1:dim(genodata)[2]) {
            dat <- genodata[, i][!is.na(genodata[, i])]
            if (!all(unique(dat) %in% c(0, 1, 2))) {
                stop("SNP has to be supplied as genotypes 0,1,2.")
            }
            if (is.na(table(dat)["1"])) {
                genocount1 <- 0
            }
            if (!is.na(table(dat)["1"])) {
                genocount1 <- table(dat)["1"][[1]]
            }
            if (is.na(table(dat)["2"])) {
                genocount2 <- 0
            }
            if (!is.na(table(dat)["2"])) {
                genocount2 <- table(dat)["2"][[1]]
            }
            nobs <- length(dat)
            MAF[i] <- (genocount1/nobs + 2 * genocount2/nobs)/2
            if (MAF[i] > 0.5) {
                warning(paste("MAF of SNV ",i ," is larger than 0.5, maybe recode alleles?",sep=""))
            }
        }
    }
    names(MAF) <- names(genodata)
    return(MAF)
}
