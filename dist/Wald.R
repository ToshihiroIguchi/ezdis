
#Wald分布
#https://ja.wikipedia.org/wiki/%E9%80%86%E3%82%AC%E3%82%A6%E3%82%B9%E5%88%86%E5%B8%83

#関数にあるWは大文字
#確率密度関数と累積分布関数でxとqにゼロ以下の値を入れるとエラーになる。
#しかし、その範囲はWald関数の下限以下で、必ずゼロになるはずなのでゼロを入れておく。
dWald <- function(x, mu, lambda) {

  if(min(x) > 0){
    #xが0よりおおきければ普通に計算
    ret <- rmutil::dinvgauss(y = x, m = mu, s = lambda)
  }else{
    
    #まず、戻り値を定義しておく
    ret <- rep(NA, length(x))
    
    #xが0以下はありえないのでゼロを入れておく
    ret[x <= 0] <- 0
    
    #xが0より大きい場合だけ計算
    ret[x > 0] <- rmutil::dinvgauss(y = x[x > 0], m = mu, s = lambda)
  }
  
  #戻り値
  return(ret)
}

pWald <- function(q, mu, lambda) {
  
  if(min(q) > 0){
    #qが0よりおおきければ普通に計算
    ret <- rmutil::pinvgauss(q = q, m = mu, s = lambda)
  }else{
    
    #まず、戻り値を定義しておく
    ret <- rep(NA, length(q))
    
    #qが0以下はありえないのでゼロを入れておく
    ret[q <= 0] <- 0
    
    #qが0より大きい場合だけ計算
    ret[q > 0] <- rmutil::pinvgauss(q = q[q > 0], m = mu, s = lambda)
  }
  
  #戻り値
  return(ret)
}

qWald <- function(p, mu, lambda) rmutil::qinvgauss(p = p, m = mu, s = lambda)



