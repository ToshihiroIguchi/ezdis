#一般化パレート分布
#ライブラリを使用しない


dGPD <- function(x, loc, scale, shape) {
  
  #戻り値のベクトルを作っておく。すべてNAを入れる。
  ret <- rep(NA, length(x))
  
  #値が計算できるはずの位置
  good.pos <- which(x > loc)
  
  #計算できる値のみを入れる
  ret[good.pos] <- (1/scale)*(1 + shape*(x[good.pos] - loc)/scale)^(-1/(shape) - 1)
  
  #戻り値
  return(ret)
  
  
}

pGPD <- function(q, loc, scale, shape) {
  
  #戻り値のベクトルを作っておく。すべてNAを入れる。
  ret <- rep(NA, length(q))
  
  #値が計算できるはずの位置
  good.pos <- which(q > loc)
  
  #計算できる値のみを入れる
  ret[good.pos] <- 1 - (1 + shape*(q[good.pos] - loc)/scale)^(-1/shape)
  
  #戻り値
  return(ret)
  
}

qGPD <- function(p, loc, scale, shape){
  ret <- (scale/shape)*((1 - p)^(-shape) - 1) + loc
  return(ret)
}
