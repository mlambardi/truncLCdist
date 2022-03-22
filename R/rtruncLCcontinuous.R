
rtruncLCcontinuous <- function(n, spec, a = -Inf, b = +Inf, ...) {
  # gammashape <- function(shape, rate = 1, scale = 1/rate) shape
  # gammascale<- function(shape, rate = 1, scale = 1/rate) scale
  auxmode <- truncLCdistoptions(paste0(spec, ".m"))
  if (is.null(auxmode)) stop("mode unspecified")
  ld <- truncLCdistoptions(paste0(spec, ".d"))
  if (is.null(ld)) stop("probability density unspecified")
  lcdf <- truncLCdistoptions(paste0(spec, ".p"))
  if (is.null(lcdf)) stop("cumulative distribution function unspecified")
  m <- auxmode(...) # base distribution's mode
  m <- pmax(pmin(m, b), a) # truncated mode, log-concavity applies
  ldm <- ld(m, ..., log=T) # log-probability of the mode
  # log-normalization constant under truncation:
  lcdfa <- lcdf(a,..., log.p=T)
  lcdfb <- lcdf(b,..., log.p=T)
  # reads like: lk = log(cdf(b) - cdf(a)) = lcdfb + log(1 - exp(-(lcdfb - lcdfa)))
  lk <- lcdfb + Rmpfr::log1mexp(abs(lcdfb - lcdfa))
  # see Devroye 1986, section 2.3, "Rejection method for log-concave densities. Exponential version"
  x <- rep(NA, n)
  i <- 1:n
  acc <- 0
  j <- 0
  # if (spec=="gamma" & any(m < 0)) {
  #   aux <- which(i & m < 0)
  #   shp <- gammashape(...)[m < 0]
  #   scl <- gammascale(...)[m < 0]
  #   a <- pmax(a, 0)
  #   b <- pmax(b, 0)
  #   y <- rtruncLCcontinuous(
  #     n=length(aux),
  #     spec="epd",
  #     a = (a/scl)^shp,
  #     b = (b/scl)^shp,
  #     beta=1/shp
  #   )
  #   acc <- attr(y, "generated")
  #   x[aux] <- scl*y^(1/shp)
  #   i <- i[-aux]
  # }
  # 25% is lower bound on acceptance rate in continuous case
  maxit <- Rmpfr::log1mexp(truncLCdistoptions("early.stopping.prob")/n)/log(1 - 0.25)
  while ((N <- length(i)) & j < maxit) {
    acc <- acc + N
    U <- 2*runif(N)
    E <- rexp(N)
    Estar <- rexp(N)
    S <- sign(runif(N) - 0.5)
    x[i] <- ifelse(U > 1, 1 + Estar, U)
    Z <- -ifelse(U > 1, E + Estar, E)
    x[i] <- m + S*x[i]*exp(lk - ldm)
    rej <- x[i] <= a | x[i] > b | Z > ld(x[i], ..., log=T) - ldm
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
