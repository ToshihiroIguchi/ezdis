#ライブラリ読み込み
library(ggplot2)
library(tibble)
library(fitdistrplus)
library(ismev)

library(FAdist) #actuar,evd, EnvStarsより先に読み込ませる
library(extraDistr)
library(evd) #evdはactuarより先に読み込ませて、evd::dgumbelなどをマスクする

library(actuar)

library(EnvStats)
library(mixtools)

library(goftest) #CVMのomega2からp-valueを計算

#Gumbel関数の読み込み
source("gumbel.R")

#ベクトルに強制変換
as.vec <- function(x){
  as.matrix(x) %>% as.vector()
}

#NULLとNAであればNULLを返す
is.null.na.null <- function(x){

  if(is.null(x)){return(NULL)}
  if(is.na(x)){return(NULL)}
  if(x == "NA"){return(NULL)}
  return(TRUE)
  
}

#ヒストグラムを描く
gg.hist <- function(vec, bw = NULL){
  
  #https://qiita.com/hoxo_b/items/13d034ab0ed60b4dca88
  
  #エラーチェック
  if(is.null(vec)){return(NULL)}
  
  #エラーチェック2
  if(!is.vector(vec)){return(NULL)}
  
  
  vec <- na.omit(vec)
  
  #確率密度計算
  dens <- density(vec)
  
  #バンド幅が指定されていないければ設定
  if(is.null(bw)){
    bw <- diff(range(vec))/20
  }
  
  #ggplot本体
  ret <- ggplot(data.frame(x = vec), aes(x = x)) +
    geom_histogram(binwidth = bw, fill = "white", color = "black") +
    geom_density(eval(bquote(aes(y=..count..*.(bw)))), fill='black', alpha=0.3) +
    xlim(range(dens$x))
  
  #戻り値
  return(ret)
}

#NULLをNAに変換
null.na <- function(x){
  if(is.null(x)){return(NA)}else{return(x)}
}

#KSのD値からp値を計算
kspval <- function(n, D, k.max = 100){
  
  #エラーチェック
  if(!is.numeric(n) || !is.numeric(D)){return(NA)}
  
  #https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/ks.test.R
  
  pkstwo.fn <- function(x, k.max = 10000){
    

    ret <- 1
    for(i in c(1:k.max)){
      ret <- ret + 2*((-1)^i)*exp(-2*(i^2)*(x^2))
    }
    
 
    #戻り値
    return(ret)
  }
  
  #p値の計算
  ret <- 1 - pkstwo.fn(sqrt(n) * D, k.max = k.max)
  
  #戻り値
  return(ret)
  
}

#CVMのomega2からp-valueを計算
cvmpval <- function(n, omega2){
  #https://github.com/cran/goftest/blob/master/R/cvmtest.R
  
  #エラーチェック
  if(!is.numeric(n) || !is.numeric(omega2)){return(NA)}
  
  #cvm test
  ret <- try(pCvM(q = omega2, n = n, lower.tail = FALSE), silent = FALSE)
  
  #エラーだったらNAを返す
  if(class(ret) == "try-error"){return(NA)}
  
  #戻り値
  return(ret)
}

#ADのAnからp-valueを計算
adpval <- function(n, An){
  
  #エラーチェック
  if(!is.numeric(n) || !is.numeric(An)){return(NA)}
  
  #cvm test
  ret <- try(pAD(q = An, n = n, lower.tail = FALSE), silent = FALSE)
  
  #エラーだったらNAを返す
  if(class(ret) == "try-error"){return(NA)}
  
  #戻り値
  return(ret)
  
  
}

#分布関数にフィッティング
fit.dist <- function(data, distr = "norm", method = "mle"){
  
  #エラーチェック
  if(is.null(data)){return(NULL)}
  
  #初期時間
  t0 <- Sys.time()
  
  #初期値
  fitdist.start <- NULL
  fix.arg <- NULL
  fitdist.lower <- -Inf
  fitdist.upper <- Inf
  optim.method <- "Nelder-Mead"
  
  #エラーの場合の戻り値を定義
  error.ret <- function(st = Sys.time()){
    ret <- list()
    class(ret) <- "fitdist.error"
    ret$CalculationTime <- st - t0
    return(ret)
  }
  
  #各分布関数の初期値を指定
  
  #対数正規分布の場合の初期値
  if(distr == "lnorm"){
    
    #最小値がゼロ以下で対数になりえない場合は空の結果を返す
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(meanlog = mean(log(data)), sdlog = sd(log(data)))
    fitdist.lower <- c(-Inf, 0)
  }
  
  #ガンベル分布の場合の初期値
  if(distr == "Gumbel"){
    gum.res <- gum.fit(data)
    fitdist.start <- list(alpha = gum.res$mle[1], scale = gum.res$mle[2])
    fitdist.lower <- c(-Inf, 0)
  }
  
  #逆ワイブル分布の場合の初期値
  if(distr == "invweibull"){
    fitdist.start <- list(shape = 1, scale = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #3母数ワイブル分布の初期値
  if(distr == "weibull3"){
    
    fitdist.start <- list(shape = 1, scale = 1, thres = min(data) - 1)
    fitdist.lower <- c(0, 0, -Inf)
    
  }
  
  #レイリー分布の初期値
  if(distr == "rayleigh"){
    
    #最小値がゼロ以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    
    fitdist.start <- list(sigma = 1)
    fitdist.lower <- c(0)
    
  }
  
  #一般化極値分布の場合の初期値
  if(distr == "gev"){
    
    #evdパッケージのgev関数
    gev.res <- try(gev.fit(data), silent = FALSE)
    
    #結果がエラーなら空の結果を返す
    if(class(gev.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }

    fitdist.start <- list(loc = gev.res$mle[1], scale = gev.res$mle[2], shape = gev.res$mle[3])
    
    fitdist.lower <- c(-Inf, 0, -Inf)
    
  }

  #一般化パレート分布の場合の初期値
  if(distr == "gpd"){
    gen.pareto.res <- try(gev.fit(data), silent  = FALSE)
    
    #結果がエラーなら空の結果を返す
    if(class(gen.pareto.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(
      loc = gen.pareto.res$mle[1], 
      scale = gen.pareto.res$mle[2], 
      shape = gen.pareto.res$mle[3]
    )
    
    fitdist.lower <- c(-Inf, 0, -Inf)
  }
  
  #指数分布の場合の初期値
  if(distr == "exp"){
    
    fitdist.lower <- c(0)
  
  }
  
  #パレート分布の場合の初期値(EnvStatsを想定)
  if(distr == "pareto"){
    fitdist.start <- list(shape = 1, location = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #タイプ2パレート分布の場合の初期値
  if(distr == "pareto2"){
    fitdist.start <- list(shape = 1, scale = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #ベータ分布の初期値
  if(distr == "beta"){
    fitdist.start <- list(shape1 = 1, shape2 = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #カイ二乗分布
  if(distr == "chisq"){
    fitdist.start <- list(df = 1)
  }
  
  #t分布の場合の初期値
  if(distr == "t"){
    fitdist.start <- list(df = 1, ncp = 1)
  }
  
  #F分布の場合の初期値
  if(distr == "f"){
    fitdist.start <- list(df1 = 1, df2 = 1, ncp = 1)
  }
  
  #2変量混合正規分布の場合の初期値
  if(distr == "normMix"){
    
    
    normalmixEM.res <- normalmixEM(data)
    
    fitdist.start <- list(
      mean1 = normalmixEM.res$mu[1], 
      sd1 = normalmixEM.res$sigma[1],
      mean2 = normalmixEM.res$mu[2],
      sd2 = normalmixEM.res$sigma[2],
      p.mix = 1-normalmixEM.res$lambda[1])
    
    
    fitdist.lower <- c(-Inf, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, Inf, Inf, 1)
    

  }
  
  #計算中の分布関数を表示
  print(distr)
  
  #フィッティングを関数に
  fitdist.fn <- function(method = method){
    
    #フィッティング
    ret <- suppressWarnings(
      try(fitdist(data = data, distr = distr, method = method, 
                  start = fitdist.start, fix.arg = fix.arg,
                  lower = fitdist.lower, upper = fitdist.upper,
                  optim.method = optim.method), silent = FALSE)
    )

    #戻り値
    return(ret)
    
  }
 
  #フィッティング
  ret <- fitdist.fn(method = method)

  #エラーならmgeで計算
  if(class(ret)[1] == "try-error"){
    ret <- fitdist.fn(method = "mge")
  }

  #エラーなら空のリストを返す
  if(class(ret)[1] == "try-error"){
    return(error.ret(Sys.time()))
  }

  #時間を書き込み
  ret$CalculationTime <- Sys.time() - t0
  
  #戻り値
  return(ret)
}

#分布関数の結果を一覧表示
summary.fit.dist <- function(data){
  
  #戻り値をデータフレームに
  ret <- NULL

  #各結果をデータフレームに変換
  for(i in 1:length(data)){

    #適合度の統計量計算
    gofstat.res <- try(gofstat(data[[i]]), silent = TRUE)
    
    #エラーでなければデータフレーム作成
    if(class(gofstat.res)[1] != "try-error"){
      
      #データフレーム作成
      df0 <- tibble(
        distr = names(data)[i][1],
        name = data[[i]]$name[1],
        AIC = null.na(data[[i]]$aic)[1],
        BIC = null.na(data[[i]]$bic)[1],
        "log likelihood" = null.na(data[[i]]$loglik)[1],
        
        #連続分布の場合
        "Kolmogorov-Smirnov statistic(D)" = null.na(gofstat.res$ks)[1],
        "Kolmogorov-Smirnov test p-value" = kspval(data[[i]]$n, gofstat.res$ks)[1],
        "Cramer-von Mises statistic(omega2)" = null.na(gofstat.res$cvm)[1],
        "Cramer-von Mises test p-value" = cvmpval(data[[i]]$n, gofstat.res$cvm)[1],
        "Anderson-Darling statistic(An)" = null.na(gofstat.res$ad)[1],
        "Anderson-Darling test p-value" = adpval(data[[i]]$n, gofstat.res$ad)[1],
        
        
        #離散分布の場合
        "Chi-squared p-value" = null.na(gofstat.res$chisqpvalue)[1],
        
        
        "Calculation time" = data[[i]]$CalculationTime[1]
      )

    }

    #結合
    if(!is.null(ret)){
      ret <- dplyr::bind_rows(ret, df0)
    }else{
      ret <- df0
    }

  }
  
  #AICの小さな順に並び替え
  ret <- ret %>% dplyr::arrange(AIC)
  
  #重複を削除
  ret <- dplyr::distinct(ret)

  #戻り値
  return(ret)
}



