#多重モードワイブル分布

#確率密度
dmultiweibull <- function(x, shape1, scale1, shape2, scale2){
  ret <- (exp(-((x/scale1)^shape1) - ((x/scale2)^shape2))) *
    ((shape1 * (x / scale1)^shape1) + (shape2 * (x / scale2)^shape2))/x
  return(ret)
}

#累積分布
pmultiweibull <- function(q, shape1, scale1, shape2, scale2){
  ret <- 1 - (exp(-(q/scale1)^shape1)) * (exp(-(q/scale2)^shape2)) 
  
   return(ret)
}

#確率点
qmultiweibull <- function(p, shape1, scale1, shape2, scale2){
  
  #戻り値
  ret <- rep(NA, length(p))
  
  #各pの値
  for(i in 1:length(p)){
    
    #最小化する関数
    #xは負の値を取らないのでexpで変換
    q.opt <- function(x){
      ret <- try(
        log((pmultiweibull(q = exp(x), 
                      shape1 = shape1, scale1 = scale1, 
                      shape2 = shape2, scale2 = scale2) - p[i])^2), 
        
        silent = TRUE)
      
      if(class(ret)[1] == "try-error"){ret <- 1e10}
      if(is.na(ret)){ret <- 1e10}
      
      return(ret)
      }
    

    
    #最適化
    #L-BFGS-B, CG
    opt.res <- optim(par = log(mean(scale1, scale2))-3 , q.opt, method = "L-BFGS-B")
    
    #求めたい値
    ret[i] <- exp(opt.res$par)
    
  }
  
  #戻り値
  return(ret)
}

