
#' Generate truncated random deviates
#'
#' @importFrom Rdpack reprompt
#' @description This function generates n random deviates that are drawn from the specified truncated distribution via Devroye sampling method for log-concave distributions. It is meant to improve on truncdist package over a subset a distributions, only exponential family models for now.
#' @param n a positive integer for the number of random deviates generated
#' @param spec a character value that specifies the underlying probability distribution
#' \itemize{
#'   \item{"binom" - stats' binomial distribution}
#'   \item{"pois" - stats' Poisson distribution}
#'   \item{"nbinom" - stats' negative binomial distribution}
#'   \item{"geom" - stats' geometric distribution}
#'   \item{"gamma" - stats' gamma distribution}
#'   \item{"norm" - stats' normal distribution}
#'   \item{"invgauss" - actuar's inverse-Gaussian distribution}
#' }
#' @param a a numeric value for the lower bound of the random variable
#' @param b a numeric value for the upper bound of the random variable
#' @param ... other arguments are are passed to the corresponding distribution-related functions
#' @param its boolean, whether classical quantile-based sampling à-la truncdist should be performed. Default is FALSE, so Devroye sampling is used instead.
#' @return A vector with one or more random deviates. An attribute reports on acceptance rates.
#' @details
#' This function returns a vector of size n, irrespective of parameters supplied.
#' \itemize{
#'  \item{Log-concave distributions are such that their logarithmic transformation is concave. A rejection sampling algorithm ensuring acceptance rate > 1/4 exists for log-concave densities \insertCite{Devroye_1986}{truncLCdist}. Its counterpart for log-concave discrete distributions ensures acceptance rate > 1/5 \insertCite{Devroye_1987}{truncLCdist}.}
#'  \item{Gamma distribution with shape parameter < 1 sampled is simulated by transforming samples from an exponential power distribution, obtained à-la Devroye. Exponential power distribution is implemented only for internal use for now, because it is not exponential class and we should provide ad hoc checks, distinct from those we meant to focus in this package.}
#' }
#' @author Euloge Clovis Kenne Pagui, \email{kenneeg@gmail.com}
#' @author Michele Lambardi di San Miniato, \email{michele.lambardidisanminiato@gmail.com}
#' @references
#'   \insertAllCited{}
#' @examples
#' rtruncLC(10000, "pois", a=15, lambda=1)
#' rtruncLC(10000, "binom", a=50, size=100, prob=0.15)
rtruncLC <- function(n, spec, a=-Inf, b=Inf, ..., its=F) {
  if (any(a > b)) {
    stop("argument a must be less-than-or-equal-to b")
  }
  if (its) {
    return(inverseTransform(n, spec, a, b, ...))
  }
  if (!is.null(e <- truncLCdistoptions(paste0(spec, ".exception.raise")))) {
    if (e(a=a, b=b, ...)) {
      e <- truncLCdistoptions(paste0(spec, ".exception.handle"))
      return(e(n=n, a=a, b=b, ...))
    }
  }
  cont <- truncLCdistoptions(paste0(spec, ".continuous"))
  if (is.null(cont)) {
    stop(paste0(spec, " not found. Start by setting up ", spec, ".continuous=F or T in options"))
  } else {
    if (cont) { # continuous
      rtruncLCcontinuous(n, spec, a, b, ...)
    } else { # discrete
      rtruncLCdiscrete(n, spec, floor(a), floor(b), ...)
    }
  }
}
