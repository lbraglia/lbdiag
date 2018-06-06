#' Confidence interval for the difference between two
#' paired proportions in a 2x2 table
#'
#' This function computes the Wald confidence interval with
#' Bonett-Price Laplace adjustment confidence interval for the
#' difference between two paired proportions. It gives the estimate as
#' well.
#' 
#' @param n12 corresponding cell frequency in 2x2 setting
#' @param n21 corresponding cell frequency in 2x2 setting
#' @param N table total frequency
#' @param alpha type I error
#'
#' @references
#'
#' Fagerland MW, Lydersen S, Laake P. Recommended tests and confidence
#'     intervals for paired binomial proportions.  Stat Med. 2014 Jul
#'     20;33(16):2850-75. doi: 10.1002/sim.6148. Epub 2014 Mar 20.
#' 
#' @examples
#' 
#' # table V fagerland etal, recommended tests and confindence intervals
#' # for paired binomial proportions  ,stat in med 2014
#' 
#' bonett_price(n12 = 1, n21 = 7, N =21)
#' @export
bonett_price <- function(n12, n21, N, alpha = 0.05){
    p12 <- (n12 + 1)/(N+2)
    p21 <- (n21 + 1)/(N+2)
    pdiff <- p12 - p21
    z <- qnorm(1 - alpha/2)
    se <- sqrt( (p12 + p21 - (pdiff)^2) / (N+2))
    est <- (n12 - n21) / N
    lower <- max(pdiff - z * se, -1L)
    upper <- min(pdiff + z * se, +1L)
    c('estimate' = est, 'lower' = lower, 'upper' = upper)
}


#' compare diagnostic accuracy measures for 2 binary tests evaluated
#' in a paired setting (each test evaluated each patients) having a
#' reference tandard
#'
#' @param test1 oldest test available
#' @param test2 new test
#' @param refstd reference standard diagnosis
#' @param alpha Type I error

da_compare <- function(test1 = NULL, test2 = NULL, refstd = NULL,
                       alpha = 0.05,
                       test1_lab  = "Test 1",
                       test2_lab  = "Test 2",
                       refstd_lab = "Refence standard")
{

    test1  <- binary_preproc(test1, "test1")
    test2  <- binary_preproc(test2, "test2")
    refstd <- binary_preproc(test2, "test2")

    db <- data.frame(test1, test2, refstd)


}
