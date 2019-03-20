#' Generate data from the Clayton copula.
#'
#' Function to generate two quantitative phenotypes \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2},
#' from the bivariate Clayton copula with standard normal
#' marginal distributions.
#'
#' @param phi Integer specifying the value of the copula parameter \eqn{\phi} for
#'            the dependence between the two generated phenotypes.
#' @param n Sample size.
#' @return A dataframe containing \code{n} observations of \eqn{Y_1}{Y1}, \eqn{Y_2}{Y2}.
#'
#' @examples
#'
#' set.seed(10)
#' dat1a <- generate_clayton_copula(n = 1000, phi = 0.5)
#' dat1b <- generate_clayton_copula(n = 1000, phi = 2)
#' dat1c <- generate_clayton_copula(n = 1000, phi = 8)
#' par(mfrow = c(3, 1))
#' plot(dat1a$Y1, dat1a$Y2, main="Clayton copula, tau = 0.2")
#' plot(dat1b$Y1, dat1b$Y2, main="Clayton copula, tau = 0.5")
#' plot(dat1c$Y1, dat1c$Y2, main="Clayton copula, tau = 0.8")
#'
#' @export
#'

generate_clayton_copula <- function(n = NULL, phi = NULL) {
    U1 <- stats::runif(n, min = 0, max = 1)
    U2 <- stats::runif(n, min = 0, max = 1)
    Y1 <- stats::qnorm(U1, mean = 0, sd = 1)
    Y2.help <- (1 - U1^(-phi) * (1 - U2^(-phi/(1 + phi))))^(-1/phi)
    tau <- phi/(phi + 2)
    Y2 <- stats::qnorm(Y2.help, mean = 0, sd = 1)
    phenodata <- data.frame(Y1 = Y1, Y2 = Y2)
    print(paste("Kendall's tau between Y1, Y2 = ", tau, sep=""))
    return(phenodata)
}
