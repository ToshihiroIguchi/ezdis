#Hyperbolic Secant Distribution
#https://rpubs.com/FJRubio/HSD

dhs <- function(x, mu, sigma){
  1/(2 * sigma) * (1/cosh(pi * (x - mu) / (2 * sigma)))
}

phs <- function(q, mu, sigma){
  (2 / pi)*atan(exp((pi * (q - mu))/ (2 * sigma)))
}

qhs <- function(p, mu, sigma){
  #https://www.wolframalpha.com/input/?i=Solve%5B%282+%2F+pi%29*atan%28exp%28%28pi+*+%28x+-+mu%29%29%2F+%282+*+sigma%29%29%29+%3DF%2C+x%5D&lang=ja
  (2 * sigma) * log(tan((pi * p)/2))/pi + mu
}

