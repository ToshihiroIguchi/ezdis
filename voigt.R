#Voigt関数
library(RcppFaddeeva)

#確率密度
dvoigt <- function(x, x0, sigma, gamma){
  
  ret <- Voigt(x = x, x0 = x0, sigma = sigma, gamma = gamma)
  return(ret)
  
  
}

#累積
pvoigt <- function(q, x0, sigma, gamma){
  
  fn <- function(x){
    dvoigt(x = x, x0 = x0, sigma = sigma, gamma = gamma)
  }
  
  #変数を作る
  ret <- rep(NA, length(q))
  
  #各pで計算
  for(i in 1:length(q)){

    #積分するとエラーがでることがあるのでtryでつつむ 
    ret.i <- try(integrate(fn, lower = -Inf, upper = q[i])$value, silent =TRUE)
    
    #エラー処理
    if(class(ret.i) == "try-error"){
      
      #積分に失敗した場合はNA
      ret[i] <- NA
      
    }else{
      
      #積分に成功したら
      ret[i] <- ret.i
    }
    
  }
  
  
  return(ret)
}

#確率点
qvoigt <- function(p, x0, sigma, gamma){

  #戻り値
  ret <- rep(NA, length(p))
  
  #初期値を決める
  start.par <- x0
  
  #各pの値
  for(i in 1:length(p)){
    
    #最小化する関数
    q.opt <- function(x){
      
      #累積確率を計算
      res.p <- pvoigt(q = x, x0 = x0, sigma = sigma, gamma = gamma)
      
      #累積確率の計算（積分）に失敗するとNAが帰る。
      if(is.na(res.p)){
        #NAなら適当なでかい数を返す。
        ret <- 1e30
      }else{
        #NAでなければ誤差の二乗和を返す
        ret <- (res.p - p[i])^2
        
      }
      
      
      #対数変換 ret=0 だと-Infになってしまうので、小さな数を足す
      ret <- log(ret + 1e-10)
      
      return(ret)
    }
    
    #最適化(L-BFGS-B)
    opt.res <- optim(par = start.par, q.opt, method = "L-BFGS-B")
    
    #BFGSは変に収束
    #CGは0.156239 secs
    #L-BFGS-Bは0.1366379 secs
    #SANNは0.789161 secs
    
    #求めたい値
    ret[i] <- opt.res$par
    
    #初期値を変更
    start.par <- opt.res$par
  }
  
  #戻り値
  return(ret)
  
  
  
}
