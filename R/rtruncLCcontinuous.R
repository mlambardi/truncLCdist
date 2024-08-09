
rtruncLCcontinuous <- function(n, spec, a = -Inf, b = +Inf, maxit=NULL, ...) {
  # gammashape <- function(shape, rate = 1, scale = 1/rate) shape
  # gammascale<- function(shape, rate = 1, scale = 1/rate) scale
  auxmode <- truncLCdistoptions(paste0(spec, ".m"))
  if (is.null(auxmode)) stop("mode unspecified")
  ld <- truncLCdistoptions(paste0(spec, ".d"))
  if (is.null(ld)) stop("probability density unspecified")
  lcdf <- truncLCdistoptions(paste0(spec, ".p"))
  if (is.null(lcdf)) stop("cumulative distribution function unspecified")
  if (b < a) stop("must be a <= b")
  m <- auxmode(...) # base distribution's mode
  m <- pmax(pmin(m, b), a) # truncated mode, log-concavity applies
  ldm <- ld(m, ..., log=T) # log-probability of the mode
  # log-normalization constant under truncation:
  lcdfa <- lcdf(a,..., log.p=T)
  lcdfb <- lcdf(b,..., log.p=T)
  # reads like: lk = log(cdf(b) - cdf(a)) = lcdfb + log(1 - exp(-(lcdfb - lcdfa)))
  lk <- lcdfb + Rmpfr::log1mexp(abs(lcdfb - lcdfa))
  if (!is.finite(lk)) {
    # can overestimate lk, which is ok, the algorithm is still exact, acceptance rate is lower
    ldm1 <- abs(truncLCdistoptions(paste0(spec, ".dld"))(m, ...))
    if ((m > a & m < b) | (ldm1 == 0)) {
      lk <- ldm + log(b-a)
    } else {
      lk <- ldm - log(ldm1) + Rmpfr::log1mexp(ldm1*(b-a))
    }
  }
  # Devroye 1986, Ch. 7, Section 2.3, "Rejection method for log-concave densities. Exponential version"
  x <- rep(NA, n)
  i <- seq(n)
  acc <- 0
  j <- 0
  # 25% is lower bound on acceptance rate in continuous case
  if (is.null(maxit)) maxit <- Rmpfr::log1mexp(truncLCdistoptions("early.stopping.prob")/n)/log(1 - 0.25)
  while ((N <- length(i)) & j < maxit) {
    if (!is.finite(lk)) break
    acc <- acc + N
    U <- 2*runif(N)
    E <- rexp(N)
    Estar <- rexp(N)
    S <- sign(runif(N) - 0.5)
    x[i] <- m + S*ifelse(U > 1, 1 + Estar, U)*exp(lk - ldm)
    Z <- -ifelse(U > 1, E + Estar, E)
    rej <- x[i] <= a | x[i] > b | Z > ld(x[i], ..., log=T) - ldm
    i <- i[ifelse(is.finite(rej), rej, T)] # otherwise, "i" may contain NAs
    j <- j + 1
  }
  if (N) x[i] <- m
  attr(x, "iterations") <- j
  attr(x, "imputingmode") <- N/n
  attr(x, "generated") <- acc
  attr(x, "acceptance.rate") <- n/acc
  attr(x, "location") <- m
  attr(x, "scale") <- exp(lk - ldm)
  return(x)
}
