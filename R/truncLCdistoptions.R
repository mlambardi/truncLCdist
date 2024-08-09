
#' Add/modify distributions
#'
#' @importFrom Rdpack reprompt
#' @description Sets the package behavior with each distribution, which can be added or modified. Supplies rtruncLC with the relevant information on how each distribution shall be sampled, based on mode (.m), probability mass/density function (.d), cumulative distribution function(.p), whether it is continuous or not (.continuous). Derived from futile.options package by \insertCite{futileopts}{truncLCdist}.
#' @param ... a character string or named argument with value, to read/write the behavior of the package
#' @param simplify a boolean, telling whether multiple options should be presented as vectors in read mode
#' @param update allows specific values to be updated dynamically rather than via named key=value pairs, can ignore it, see futile.options for more details
#' @references
#'   \insertAllCited{}
#' @examples
#' truncLCdistoptions("binom.d")
#' truncLCdistoptions(binom.d=stats::dbinom)
truncLCdistoptions <- futile.options::OptionsManager(
  'truncLCdistoptions',
  default=list(
    # log-probability of no early stopping under non-problematic settings
    early.stopping.prob=1e-5,

    # binomial
    binom.continuous=F,
    binom.m=function(size, prob) floor(prob*(size + 1)),
    binom.d=Rmpfr::dbinom,
    binom.p=pbinom,
    binom.q=qbinom,
    binom.dld=function(x, size, prob) -digamma(x + 1) + digamma(size - x + 1) + log(prob) - log(1 - prob),

    # Poisson
    pois.continuous=F,
    pois.m=function(lambda) floor(lambda),
    pois.d=Rmpfr::dpois,
    pois.p=ppois,
    pois.q=qpois,
    pois.dld=function(x, lambda) log(lambda) - digamma(x + 1),

    # negative binomial
    nbinom.continuous=F,
    nbinom.m=function(size, prob) pmax(floor((size - 1)*(1 - prob)/prob), 0),
    nbinom.d=Rmpfr::dnbinom,
    nbinom.p=pnbinom,
    nbinom.q=qnbinom,
    nbinom.dld=function(x, size, prob) digamma(x + size) - digamma(x + 1) + log(1 - p),

    # geometric
    geom.continuous=F,
    geom.m=function(prob) 0,
    geom.d=dgeom,
    geom.p=pgeom,
    geom.q=qgeom,
    geom.dld=function(x, prob) log(1-prob),

    # gamma
    gamma.continuous=T,
    gamma.exception.raise=function(shape, rate = 1, scale = 1/rate, a=-Inf, b=+Inf) {
      any(shape < 1)
    },
    gamma.exception.handle=function(n, shape, rate = 1, scale = 1/rate, a=-Inf, b=+Inf, maxit=NULL) {
      a <- pmax(a, 0)
      out <- rtruncLC(n=n, spec="epd", beta=1/shape, a=(a/scale)^shape, b=(b/scale)^shape, maxit=maxit)
      aux <- attributes(out)
      out <- scale*out^(1/shape)
      attributes(out) <- aux
      attr(out, "note") <- "generated gamma w shape < 1 via epd"
      return(out)
    },
    gamma.m=function(shape, rate = 1, scale = 1/rate) {
      if (any(shape < 1)) {
        stop("shape < 1 unsupported")
      } else {
        (shape - 1)/rate
      }
    },
    gamma.d=Rmpfr::dgamma,
    gamma.p=pgamma,
    gamma.q=qgamma,
    gamma.dld=function(x, shape, rate = 1, scale = 1/rate) (shape - 1)/x - rate,

    # exponential
    exp.continuous=T,
    exp.m=function(rate = 1) 0,
    exp.d=dexp,
    exp.p=pexp,
    exp.q=qexp,
    exp.dld=function(x, rate=1) -rate,

    # normal
    norm.continuous=T,
    norm.m=function(mean = 0, sd = 1) mean,
    norm.d=Rmpfr::dnorm,
    norm.p=Rmpfr::pnorm,
    norm.q=qnorm,
    norm.dld=function(x, mean = 0, sd = 1) (mean - x)/sd^2,

    # inverse-Gaussian
    invgauss.continuous=T,
    invgauss.m=function(mean, shape = 1, dispersion = 1/shape) mean*((1 + 9/4*(mean*dispersion)^2)^0.5 - 3/2*mean*dispersion),
    invgauss.d=actuar::dinvgauss,
    invgauss.p=actuar::pinvgauss,
    invgauss.q=actuar::qinvgauss,
    invgauss.dld=function(x, mean, shape = 1, dispersion = 1/shape) -1.5/x - 0.5/mean^2/dispersion + 0.5/dispersion/x^2,

    # exponential power distribution
    epd.continuous=T,
    epd.m=function(beta=1, mu=0, sigma=1, thr=0, betathr=Inf) {
      if (any(beta < 1)) {
        stop("beta < 1 unsupported")
      } else {
        mu
      }
    },
    epd.d=function(x, beta=1, mu=0, sigma=1, thr=0, betathr=Inf, z=(x-mu)/sigma, log=F) {
      aux <- -abs(z)^beta - lgamma(1 + 1/beta) - log(2*sigma)
      if (log) aux else exp(aux)
    },
    epd.p=function(q, beta=1, mu=0, sigma=1, thr=1-exp(-0.2*(log(beta) - 4.8)^2), betathr=140, z=(q-mu)/sigma, log.p=F) {
      f <- function(z, beta) {
        aux <- pgamma(abs(z)^beta, shape=1/beta, log.p = T)
        ifelse(z <= 0, Rmpfr::log1mexp(-aux), Rmpfr::log1pexp(aux)) - log(2)
      }
      logy <- f(z, beta)
      if (length(i <- abs(z) < thr & beta > betathr) > 0) {
        fa <- exp(f(-thr, beta))
        fb <- exp(f(+thr, beta))
        logy[i] <- log(fa + (fb - fa)*(z[i] + thr)/(2*thr))
      }
      if (log.p) logy else exp(logy)
    },
    epd.q=function(x, beta=1, mu=0, sigma=1, thr=0, betathr=Inf) {
      stop("stub, not yet implemented")
    },
    epd.dld=function(x, beta=1, mu=0, sigma=1, thr=0, betathr=Inf, z=(x-mu)/sigma) {
      -sign(z)*abs(z)^(beta-1)/sigma
    }
  )
)
