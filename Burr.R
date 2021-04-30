#https://en.wikipedia.org/wiki/Burr_distribution
#https://www.wolframalpha.com/input/?i=solve%5BF+%3D+1-%281%2B%28x%2Fl%29%5Ec%29%5E%28-k%29%2C+x%5D&lang=ja

dBurr <- function(x, c, k, lambda){
  ((c * k / lambda) * (x / lambda)^(c - 1))* ((1 + (x / lambda)^c)^(-k - 1))
}

pBurr <- function(q, c, k, lambda){
  1 - (1 + (q / lambda)^c)^(-k)
}

qBurr <- function(p, c, k, lambda){
  lambda * ((1 - p)^(-1 / k) - 1)^(1 / c)
}



