#boxcox分布
#https://en.wikipedia.org/wiki/Box%E2%80%93Cox_distribution


#確率密度
dboxcox <- function(x, m, s = 1, f = 1){
  #I(f<0)を計算
  If0 <- if(f < 0){1}else{0}
  
  #box-cox分布の確率密度
  ret <- 1/((1 - If0 - sign(f) * pnorm(0, m, sqrt(s))) * sqrt(2 * pi * s^2)) * exp(-(1/(2 * s^2)) * ((x^f)/f - m)^2)
  
  #戻り値
  return(ret)
}

#累積分布
pboxcox <- function(q, m, s = 1, f = 1){
  
  norm <- sign(f) * pnorm(0, m, sqrt(s))
  ret <- (pnorm((q^f)/f, m, sqrt(s)) - (f > 0) * norm)/(1 - (f < 0) - norm)
  
  return(ret)
}

#確率点
qboxcox <- function(p, m, s = 1, f = 1){
  #エラーチェック
  if(is.null(p) || is.null(m)){return(NULL)}
  
  #retを作る
  ret <- vector()
  
  #初期値
  start.par <- m
  
  for(i in 1:length(p)){
    
    #最適化する関数
    q.opt <- function(x){
      ret <- (pboxcox(q = x, m = m, s = s, f = f) - p[i])^2 
      return(ret)
    }
    
    
    #最適化
    opt.res <- optim(par = start.par, fn = q.opt, method = "CG")
    
    #求めたい値
    ret[i] <- opt.res$par
    
    #初期値を変更
    #start.par <- opt.res$par
    
  }
  
  
  #戻り値
  return(ret)
  
  
  
}

