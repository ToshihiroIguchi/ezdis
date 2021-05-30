#アーラン分布

derlang <- function(x, k, mu) dgamma(x, shape = (k), rate = 1/mu)
perlang <- function(q, k, mu) pgamma(q, shape = (k), rate = 1/mu)
qerlang <- function(p, k, mu) qgamma(p, shape = (k), rate = 1/mu)









