

library(normalp)

#exponential power distribution
#https://cran.r-project.org/web/packages/normalp/normalp.pdf
dnormp2 <- function(x, mu, sigmap, shape) dnormp(x = x, mu = mu, sigmap, p = shape)
pnormp2 <- function(q, mu, sigmap, shape) pnormp(q = q, mu = mu, sigmap, p = shape)
qnormp2 <- function(p, mu, sigmap, shape) qnormp(pr = p, mu = mu, sigmap, p = shape)




