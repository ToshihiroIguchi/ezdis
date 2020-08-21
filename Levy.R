#Levy分布

#ライブラリ読み込み
library(rmutil)

#確率密度
dLevy <- function(x, m, s){
  
  #戻り値のベクトルを作っておく。すべてNAを入れる。
  ret <- rep(0, length(x))
  
  #値が計算できるはずの位置
  good.pos <- which(x > m)
  
  #計算できる値のみを入れる
  ret[good.pos] <- rmutil::dlevy(y = x[good.pos], m = m, s = s)
  
  #戻り値
  return(ret)
}

#累積確率
pLevy <- function(q, m, s){
  
  #戻り値のベクトルを作っておく。すべてNAを入れる。
  ret <- rep(0, length(q))
  
  #値が計算できるはずの位置
  good.pos <- which(q > m)
  
  #計算できる値のみを入れる
  ret[good.pos] <- rmutil::plevy(q = q[good.pos], m = m, s = s)
  
  #戻り値
  return(ret)
  
  
  
  
}

#確率点
qLevy <- function(p, m, s) rmutil::qlevy(p = p, m = m, s = s)




