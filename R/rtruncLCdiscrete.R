
rtruncLCdiscrete <- function(n, spec, a = -Inf, b = +Inf, ...) {
  auxmode <- truncLCdistoptions(paste0(spec, ".m"))
  if (is.null(auxmode)) stop("mode unspecified")
  ld <- truncLCdistoptions(paste0(spec, ".d"))
  if (is.null(ld)) stop("probability mass function unspecified")
  lcdf <- truncLCdistoptions(paste0(spec, ".p"))
  if (is.null(lcdf)) stop("cumulative distribution function unspecified")
  m <- auxmode(...)
  m <- pmax(pmin(m, b), a + 1) # truncated mode, due to concavity
  ldm <- ld(m, ..., log=T) # log-density at the mode
  # log-normalization constant under truncation:
  lcdfb <- lcdf(b, ..., log.p=T)
  lcdfa <- lcdf(a, ..., log.p=T)
  # reads like: lk = log(cdf(b) - cdf(a)) = lcdfb + log(1 - exp(-(lcdfb - lcdfa)))
  lk <- lcdfb + Rmpfr::log1mexp(abs(lcdfb - lcdfa))
  # see Devroye 1987, algorithm n. 2
  w <- 1 + exp(ldm - lk)/2
  thr <- w/(1 + w)
  # following is shared with the continuous case
  x <- rep(NA, n)
  i <- 1:n
  acc <- 0
  j <- 0
  # 20% is lower bound on acceptance rate in discrete case
  maxit <- Rmpfr::log1mexp(truncLCdistoptions("early.stopping.prob")/n)/log(1 - 0.2)
  while ((N <- length(i)) & j < maxit) {
    acc <- acc + N
    U <- runif(N)
    V <- runif(N)
    W <- runif(N)
    E <- rexp(N)
    S <- 2*rbinom(N,1,0.5)-1
    y <- ifelse(U > thr, w + E, w*V)
    x[i] <- ifelse(is.finite(y), m + S*round(y*exp(lk - ldm)), Inf)
    # with log-concave distributions, infinite is almost surely rejected.
    rej <- x[i] <= a | x[i] > b | log(W) + pmin(w - y, 0) > ld(x[i], ..., log=T) - ldm
    rej <- ifelse(is.finite(rej), rej, T) # otherwise, "i" may contain NAs
    i <- i[rej]
    j <- j + 1
  }
  if (N > 0) {
    x[i] <- rep(m, length.out=N)
  }
  attr(x, "iterations") <- j
  attr(x, "imputingmode") <- N/n
  attr(x, "generated") <- acc
  attr(x, "acceptance.rate") <- n/acc
  return(x)
}
