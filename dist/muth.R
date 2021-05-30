#Muth分布
#https://www.researchgate.net/publication/279164529_On_the_Muth_Distribution

#確率密度関数
dmuth <- function(x, alpha){
  (exp(alpha * x) - alpha) * exp(alpha * x - 1/alpha * (exp(alpha * x) - 1))
}

#累積分布関数
pmuth <- function(q, alpha){
  1 - exp(alpha * q - 1/alpha*(exp(alpha * q) - 1))
}

#累積分布関数の逆関数
qmuth <- function(p, alpha){
  #https://www.wolframalpha.com/input/?i=Solve%5B+1+-+exp%28alpha+*+x+-+1%2Falpha*%28exp%28alpha+*+x%29+-+1%29%29+%3D+F%2C+x%5D&lang=ja
  (alpha * log(1-p) - alpha * lambert_W0((exp(-1/alpha) * (p - 1))/alpha) - 1)/(alpha^2)
}






