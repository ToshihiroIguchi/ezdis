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

#誤差関数
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)

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
pnormmixn <- function(x, mean, sd, rate){
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
    ret <- ret + rate[i]*pnorm(x, mean[i], sd[i])
  }
  
  #戻り値
  return(ret)
}

#確率点
qnormmixn <- function(p, mean, sd, rate){
  
  
  
}




dnormmixn(0.5, c(1,5), c(1,2),c(0.5,0.5))





qnormmixn <- function(){
  
  
  #https://www.medi-08-data-06.work/entry/2018/12/18/232204
}