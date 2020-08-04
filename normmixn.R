#多変量混合正規分布

#エラーチェック
normmixn.chk <- function(mean, sd, rate){
  #ベクトルか
  if(!is.vector(mean) || !is.vector(sd) || !is.vector(rate)){
    stop("mean, sd, and rate must be vectors.")
  }
  
  #長さが同じか
  if(length(mean) != length(sd) || length(mean) != length(rate)){
    stop("mean, sd, and rate vectors must have the same length.")
  }
  
  #rateの最小値が0以上か
  if(min(rate) < 0){
    stop("Minimum value of rate is less than 0.")}
  
  #rateの合計が0より大きいか
  if(sum(rate) <= 0){
    stop("Sum of rates is greater than 0.")
  }
  
}

#比率を正規化
norm.rate <- function(rate){
  ret <- rate/sum(rate)
}

#確率密度
dnormmixn <- function(x, mean, sd, rate){
  #エラーチェック
  if(!is.numeric(x)){stop("x is not numeric.")}
  
  #エラーチェック
  normmixn.chk(mean, sd, rate)

  #比率の正規化
  rate <- norm.rate(rate)
  
  #初期値
  ret <- rep(0, length(x))
  
  #足し合わせ
  for(i in 1:length(mean)){
    ret <- ret + rate[i]*dnorm(x, mean[i], sd[i])
  }
  
  #戻り値
  return(ret)
}

#累積分布
pnormmixn <- function(q, mean, sd, rate){
  #エラーチェック
  if(!is.numeric(q)){stop("x is not numeric.")}
  
  #エラーチェック
  normmixn.chk(mean, sd, rate)
  
  #比率の正規化
  rate <- norm.rate(rate)
  
  #初期値
  ret <- rep(0, length(q))
  
  #足し合わせ
  for(i in 1:length(mean)){
    ret <- ret + rate[i]*pnorm(q, mean[i], sd[i])
  }
  
  #戻り値
  return(ret)
}

#確率点
qnormmixn <- function(p, mean, sd, rate){
  
  #戻り値
  ret <- rep(NA, length(p))
  
  #初期値を決める
  start.par <- min(mean)

  #各pの値
  for(i in 1:length(p)){

    #最小化する関数
    q.opt <- function(x){
      
      ret <- (pnormmixn(q = x, mean = mean, sd = sd, rate = rate) - p[i])^2
      
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


#2変量混合正規分布の確率密度
dnormmix2 <- function(x, mean1, sd1, rate1, mean2, sd2, rate2){
  dnormmixn(x, mean = c(mean1, mean2), 
            sd = c(sd1, sd2), rate = c(rate1, rate2))
  
}

#3変量混合正規分布の確率密度
pnormmix2 <- function(p, mean1, sd1, rate1, mean2, sd2, rate2){
  pnormmixn(p, mean = c(mean1, mean2), 
            sd = c(sd1, sd2), rate = c(rate1, rate2))
  
}

#3変量混合正規分布の確率密度
qnormmix2 <- function(q, mean1, sd1, rate1, mean2, sd2, rate2){
  qnormmixn(q, mean = c(mean1, mean2), 
            sd = c(sd1, sd2), rate = c(rate1, rate2))
  
}



#3変量混合正規分布の確率密度
dnormmix3 <- function(x, mean1, sd1, rate1, mean2, sd2, rate2, mean3, sd3, rate3){
  dnormmixn(x, mean = c(mean1, mean2, mean3), 
                   sd = c(sd1, sd2, sd3), rate = c(rate1, rate2, rate3))

}

#3変量混合正規分布の確率密度
pnormmix3 <- function(p, mean1, sd1, rate1, mean2, sd2, rate2, mean3, sd3, rate3){
  pnormmixn(p, mean = c(mean1, mean2, mean3), 
            sd = c(sd1, sd2, sd3), rate = c(rate1, rate2, rate3))
  
}

#3変量混合正規分布の確率密度
qnormmix3 <- function(q, mean1, sd1, rate1, mean2, sd2, rate2, mean3, sd3, rate3){
  qnormmixn(q, mean = c(mean1, mean2, mean3), 
            sd = c(sd1, sd2, sd3), rate = c(rate1, rate2, rate3))
  
}


#4変量混合正規分布の確率密度
dnormmix4 <- function(x, mean1, sd1, rate1, mean2, sd2, rate2, 
                      mean3, sd3, rate3, mean4, sd4, rate4){
  dnormmixn(x, mean = c(mean1, mean2, mean3, mean4), 
            sd = c(sd1, sd2, sd3, sd4), rate = c(rate1, rate2, rate3, rate4))
  
}

#4変量混合正規分布の確率密度
pnormmix4 <- function(p, mean1, sd1, rate1, mean2, sd2, rate2, 
                      mean3, sd3, rate3, mean4, sd4, rate4){
  pnormmixn(p, mean = c(mean1, mean2, mean3, mean4), 
            sd = c(sd1, sd2, sd3, sd4), rate = c(rate1, rate2, rate3, rate4))
  
}

#4変量混合正規分布の確率密度
qnormmix4 <- function(q, mean1, sd1, rate1, mean2, sd2, rate2, 
                      mean3, sd3, rate3, mean4, sd4, rate4){
  qnormmixn(q, mean = c(mean1, mean2, mean3, mean4), 
            sd = c(sd1, sd2, sd3, sd4), rate = c(rate1, rate2, rate3, rate4))
}









