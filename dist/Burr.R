#https://en.wikipedia.org/wiki/Burr_distribution
#https://www.wolframalpha.com/input/?i=Solve%5BF+%3D+1-%281%2Bx%5Ec%29%5E%28-k%29%2C+x%5D&lang=ja



dBurr <- function(x, k, c){
  c*k*(x ^ (c - 1))/((1 + x^c)^(k + 1))
}

pBurr <- function(q, k, c){
 1 - (1 + q^c)^(-k) 
}

qBurr <- function(p, k, c){
  ((1 - p)^(-1/k) - 1)^(1/c)
}

