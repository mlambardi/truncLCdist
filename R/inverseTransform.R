
inverseTransform <- function(n, spec, a, b, ...) {
  qfun <- truncLCdistoptions(paste0(spec, ".q"))
  pfun <- truncLCdistoptions(paste0(spec, ".p"))
  pa <- pfun(a, ...)
  pb <- pfun(b, ...)
  qfun(pmax(pmin(runif(n)*(pb - pa) + pa, 1), 0), ...)
}
