#多重モードワイブル分布

#確率密度
dmultiweibull <- function(x, shape1, scale1, shape2, scale2, rate){
  ret <- rate*dweibull(x, shape = shape1, scale = scale1) + 
    (1 - rate)*dweibull(x, shape = shape2, scale = scale2)
  return(ret)
}

#累積分布
pmultiweibull <- function(q, shape1, scale1, shape2, scale2, rate){
  ret <- rate*pweibull(q, shape = shape1, scale = scale1) + 
    (1 - rate)*pweibull(q, shape = shape2, scale = scale2)
  return(ret)
}

#確率点
qmultiweibull <- function(p, shape1, scale1, shape2, scale2, rate){
  
  #戻り値
  ret <- rep(NA, length(p))
  
  #各pの値
  for(i in 1:length(p)){
    
    #最小化する関数
    q.opt <- function(x){
      (pmultiweibull(q = x, shape1 = shape1, scale1 = scale1, 
                     shape2 = shape2, scale2 = scale2, rate = rate) - p[i])^2
      }
    
    #最適化
    opt.res <- optim(par = mean(c(scale1, scale2)), q.opt, method = "L-BFGS-B")
    
    #求めたい値
    ret[i] <- opt.res$par
    
  }
  
  #戻り値
  return(ret)
}


