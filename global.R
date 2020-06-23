#ライブラリ読み込み
library(ggplot2)
library(fitdistrplus)
library(ismev)

library(FAdist) #actuar,evd, EnvStarsより先に読み込ませる
library(extraDistr)
library(evd) #evdはactuarより先に読み込ませて、evd::dgumbelなどをマスクする

library(actuar)

library(EnvStats)
library(mixtools)


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
  
  #各分布関数の初期値を指定
  
  #ガンベル分布の場合の初期値
  if(distr == "gumbel"){
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
    fitdist.start <- list(sigma = 1)
    fitdist.lower <- c(0)
    
  }
  
  
  
  #一般化極値分布の場合の初期値
  if(distr == "gev"){
    gev.res <- gev.fit(data)
    fitdist.start <- list(loc = gev.res$mle[1], scale = gev.res$mle[2], shape = gev.res$mle[3])
  }

  #一般化パレート分布の場合の初期値
  if(distr == "gpd"){
    gen.pareto.res <- gev.fit(data)
    
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

  #フィッティング
  ret <- try(fitdist(data = data, distr = distr, method = method, 
                     start = fitdist.start, fix.arg = fix.arg,
                     lower = fitdist.lower, upper = fitdist.upper,
                     optim.method = optim.method), silent = FALSE)
  
  #エラーならmgeで計算
  if(class(ret)[1] == "try-error"){
    ret <- try(fitdist(data = data, distr = distr, method = "mge", 
                       start = fitdist.start, fix.arg = fix.arg,
                       lower = fitdist.lower, upper = fitdist.upper,
                       optim.method = optim.method
                       ), silent = FALSE)
  }
  

  #エラーならNAを返す
  if(class(ret)[1] == "try-error"){
    ret <- list()
    class(ret) <- "fitdist.error"
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
      
      df0 <- data.frame(
        distr = names(data)[i],
        name = data[[i]]$name,
        AIC = data[[i]]$aic,
        BIC = data[[i]]$bic,
        log.likelihood = data[[i]]$loglik,
        Kolmogorov.Smirnov.statistic = gofstat.res$ks,
        Cramer.von.Mises.statistic = gofstat.res$cvm,
        Anderson.Darling.statistic = gofstat.res$ad,
        CalculationTime = data[[i]]$CalculationTime
      )
      
      #結合
      if(!is.null(ret)){
        ret <- ret %>% rbind(df0)
      }else{
        ret <- df0
      }
    }

  }
  
  #AICの小さな順に並び替え
  ret <- ret %>% dplyr::arrange(AIC)
  
  #戻り値
  return(ret)
}



