#https://en.wikipedia.org/wiki/Lomax_distribution
#https://www.wolframalpha.com/input/?i=Solve%5BF+%3D+1+-%281%2Bx%2Fr%29%5E%28-a%29%2C+x%5D&lang=ja


dLomax <- function(x, alpha, lambda){
  (alpha/lambda)*(1 + x/lambda)^(-(alpha + 1))
}

pLomax <- function(q, alpha, lambda){
  1 - (1 + q/lambda)^(-alpha)
}

qLomax <- function(p, alpha, lambda){
  lambda*(-((1 - p)^(1/alpha) - 1))*((1 - p)^(-1/alpha))
}

