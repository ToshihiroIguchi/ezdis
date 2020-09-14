#Laplace

dLaplace <- function(x, mu, sigma) extraDistr::dlaplace(x, mu = mu, sigma = sigma)
pLaplace <- function(q, mu, sigma) extraDistr::plaplace(q, mu = mu, sigma = sigma)
qLaplace <- function(p, mu, sigma) extraDistr::qlaplace(p, mu = mu, sigma = sigma)


