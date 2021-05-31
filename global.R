#ライブラリ読み込み
library(ggplot2)
library(tibble)
library(readr)
library(readxl)
library(R.utils)
library(fitdistrplus)
library(ismev)
library(FAdist) #actuar,evd, EnvStarsより先に読み込ませる
library(extraDistr)
library(evd) #evdはactuarより先に読み込ませて、evd::dgumbelなどをマスクする
library(actuar)
library(EnvStats)
library(mixtools)
library(RcppFaddeeva)
library(goftest) #CVMのomega2からp-valueを計算
library(rmutil)
library(PearsonDS)
library(gsl)
library(ExtDist)
library(lmomco)
library(hydroApps)

library(normalp)

library(triangle)

#指定したディレクトリのファイルを一括で読み込む
source.dir <- function(dir){
  
  #指定されたディレクトリのファイル一覧を取得
  files.vec <- list.files(dir)
  
  for(filename in files.vec){
    
    source(paste0(dir, "/", filename))
  }
  
}

#統計分布のファイルをdistディレクトリから一括で呼び出し
source.dir("dist")

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
  
  #欠損値を除く
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

#整数か判断(型チェックではない)
is_integer <- function(vec){
  
  #そもそも数値でなかったらFALSE
  if(!is.numeric(vec)){return(FALSE)}
    
  #誤差二乗和
  e2 <- sum((vec - round(vec, 0))^2)
  
  #戻り値
  if(e2 == 0){return(TRUE)}else{return(FALSE)}
  
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
fit.dist <- function(data, distr = "norm", method = "mle", timeout = 10){
  
  #エラーチェック
  if(is.null(data)){return(NULL)}
  
  #尤度計算でおかしいと判断する対数尤度の値
  loglik.th <- 1e30
  
  
  #初期時間
  t0 <- Sys.time()
  
  #初期値
  fitdist.start <- NULL
  fix.arg <- NULL
  fitdist.lower <- -Inf
  fitdist.upper <- Inf
  optim.method <- "Nelder-Mead"
  
  data.length.th <- 10000
  
  
  #エラーの場合の戻り値を定義
  error.ret <- function(st = Sys.time()){
    ret <- list()
    class(ret) <- "fitdist.error"
    ret$CalculationTime <- st - t0
    return(ret)
  }
  
  #初期値推定用のデータサンプリング
  if(length(data) > data.length.th){
    set.seed(108)
    sample.data <- sample(data, size = data.length.th)
  }else{
    sample.data <- data
  }

  
  #各分布関数の初期値を指定
  
  #対数正規分布の場合の初期値
  if(distr == "lnorm"){
    
    #最小値がゼロ以下で対数になりえない場合は空の結果を返す
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    #fitdist.start <- list(meanlog = mean(log(data)), sdlog = sd(log(data)))
    #fitdist.lower <- c(-Inf, 0)
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
  
  #ワイブル分布の場合の初期値
  if(distr == "weibull"){
    
    
    #最小値がゼロ以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    
    fitdist.lower <- c(0, 0)
  }
  
  #3母数ワイブル分布の初期値
  if(distr == "weibull3"){
    
    fitdist.start <- list(shape = 1, scale = 1, thres = min(data) - 1)
    fitdist.lower <- c(0, 0, -Inf)
    
  }
  
  #逆ワイブル分布の場合の初期値
  if(distr == "rweibull"){
    
    fitdist.start <- list(loc = max(data) + 1, scale = 1, shape = 1)

    fitdist.lower <- c(-Inf, 0, 0)
  }
  
  #多重モードワイブル分布の初期値
  if(distr == "multiweibull"){
    
    #最小値がゼロ以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    #小さい値と大きな値にわける
    data.1 <- sort(data)[c(1:round(length(data)/2, 0))]
    data.2 <- sort(data)[c(round(length(data)/2, 0):length(data))]
    
    #小さい値と大きな値でWeibull分布にフィッティングしてパラメータ推定
    wp1 <- fitdistrplus::fitdist(data.1, "weibull")$estimate
    wp2 <- fitdistrplus::fitdist(data.2, "weibull")$estimate

    #推定したパラメータで二つのワイブルの初期値をとする
    fitdist.start <- list(shape1 = data.1[1], scale1 = data.1[2], 
                          shape2 = data.2[1], scale2 = data.2[2])
    
    #上限と下限設定
    fitdist.lower <- c(0, 0, 0, 0)
    fitdist.upper <- c(Inf, Inf, Inf, Inf)
    
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
    gev.res <- try(gev.fit(data, show = FALSE), silent = FALSE)
    
    #結果がエラーなら空の結果を返す
    if(class(gev.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }

    fitdist.start <- list(loc = gev.res$mle[1], scale = gev.res$mle[2], shape = gev.res$mle[3])
    
    fitdist.lower <- c(-Inf, 0, -Inf)
    
  }

  #一般化パレート分布の場合の初期値
  if(distr == "GPD"){
    gen.pareto.res <- try(gpd.fit(data, show = FALSE, threshold = 2), silent  = TRUE)
    
    #結果がエラーなら空の結果を返す
    if(class(gen.pareto.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    #対数尤度の計算
    dgpd.loglikelihood  <- function(x, loc, scale, shape){
      ret <- sum(log(dGPD(x, loc, scale, shape)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
      
      
    }
    
    #パラメータから対数尤度を求める関数
    dgpd.opt <- function(x){
      dgpd.loglikelihood(x = data, x[1], x[2], x[3])
    }
    
    #対数尤度を最大化
    gpd.opt <- optim(par = c(min(0, data) - 0.1 , 1, 2), 
                     fn = dgpd.opt, control = list(fnscale  = -1))
    
    #初期値を定義
    fitdist.start <- list(loc = gpd.opt$par[1], 
                          scale = gpd.opt$par[2], 
                          shape = gpd.opt$par[3])
    
    
    #最小値と最大値
    fitdist.lower <- c(-Inf, 0, -Inf)
    fitdist.upper <- c(min(data), Inf, Inf)
    
  }
  
  #指数分布の場合の初期値
  if(distr == "exp"){
    
    
    #最小値がゼロより小さい場合は空の結果を返す
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }
    
    fitdist.lower <- c(0)
  
  }
  
  #exponential power distribution
  if(distr == "normp2"){
    fitdist.start <- list(mu = 0, sigmap = 1, shape = 2)
    fitdist.lower <- c(-Inf, 1e-10, 1)
  }
  
  #Wald分布
  if(distr == "Wald"){
    
    #最小値がゼロより小さい場合は空の結果を返す
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(mu = 1, lambda = 1)
    fitdist.lower <- c(1e-10, 1e-10)
  }
  
  #シングルパラメータパレート分布の場合の初期値
  if(distr == "pareto1"){
    
    #最小値が0以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    #対数尤度の計算
    dpareto1.ll <- function(x, shape, min){
      if(min(shape, min) <= 0){return(-Inf)}
      ret <- sum(log(dpareto1(x = x, shape = shape, min = min)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpareto1.opt <- function(x){
      ret <- dpareto1.ll(x = data, shape = x[1], min = x[2])
      return(ret)
    }
    
    #対数尤度を最大化
    gpd.opt <- optim(par = c(1, max(min(data), 1e-10)), 
                     fn = dpareto1.opt, control = list(fnscale  = -1))
    
    fitdist.start <- list(shape = gpd.opt$par[1], min = gpd.opt$par[2])
    
    fitdist.lower <- c(0, -Inf)
    fitdist.upper <- c(Inf, Inf)
  }
  
  #パレート分布の場合の初期値(actuarを想定)
  if(distr == "pareto_ac"){
    
    #最小値が0以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    #対数尤度の計算
    dpareto.ll <- function(x, shape, scale){
      if(min(shape, scale) <= 0){return(-Inf)}
      ret <- sum(log(dpareto_ac(x = x, shape = shape, scale = scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpareto.opt <- function(x){
      ret <- dpareto.ll(x = data, shape = x[1], scale = x[2])
      return(ret)
    }

    #対数尤度を最大化
    pareto.opt <- optim(par = c(1, 1), 
                     fn = dpareto.opt, control = list(fnscale  = -1))


    fitdist.start <- list(shape = pareto.opt$par[1], scale = pareto.opt$par[2])
    fitdist.lower <- c(1e-10, 1e-10)
    fitdist.upper <- c(Inf, Inf)

  }
  
  #タイプ2パレート分布の場合の初期値
  if(distr == "pareto2"){

    #対数尤度の計算
    dpareto2.ll <- function(x, min, shape, scale){
      ret <- sum(log(dpareto2(x = x, min = min, shape = shape, scale = scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpareto2.opt <- function(x){
      ret <- dpareto2.ll(x = data, min = x[1], shape = x[2], scale = x[3])
      return(ret)
    }
    
    #対数尤度を最大化
    gpd.opt <- optim(par = c(min(data) - 1, 2, 2), 
                     fn = dpareto2.opt, control = list(fnscale  = -1))

    fitdist.start <- list(min = gpd.opt$par[1], 
                          shape = gpd.opt$par[2], scale = gpd.opt$par[3])
    fitdist.lower <- c(-Inf, 0, 0)
  }
  
  #タイプ3パレート分布の場合の初期値
  if(distr == "pareto3"){

    #対数尤度の計算
    dpareto3.ll <- function(x, min, shape, scale){
      ret <- sum(log(dpareto3(x = x, min = min, shape = shape, scale = scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpareto3.opt <- function(x){
      ret <- dpareto3.ll(x = data, min = x[1], shape = x[2], scale = x[3])
      return(ret)
    }
    
    #対数尤度を最大化
    pa3.opt <- optim(par = c(min(data) - 1, 2, 2), 
                     fn = dpareto3.opt, control = list(fnscale  = -1))
    
    fitdist.start <- list(min = pa3.opt$par[1], 
                          shape = pa3.opt$par[2], scale = pa3.opt$par[3])
    
    fitdist.lower <- c(-Inf, 0, 0)
  }
  
  #タイプ4パレート分布の場合の初期値
  if(distr == "pareto4"){
    
    #対数尤度の計算
    dpareto4.ll <- function(x, min, shape1, shape2, scale){
      ret <- sum(log(dpareto4(x = x, min = min, shape1 = shape1, shape2 = shape2, scale = scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpareto4.opt <- function(x){
      ret <- dpareto4.ll(x = data, min = x[1], shape1 = x[2], shape2 = x[3], scale = x[4])
      return(ret)
    }
    
    #対数尤度を最大化
    pa4.opt <- optim(par = c(min(data), 1, 1, 1), 
                     fn = dpareto4.opt, control = list(fnscale  = -1))

    fitdist.start <- list(min = pa4.opt$par[1], 
                          shape1 = pa4.opt$par[2], shape2 = pa4.opt$par[3], 
                          scale = pa4.opt$par[4])
    
    fitdist.lower <- c(-Inf, 0, 0, 0)
  }
  
  
  #Lomax分布の場合の初期値
  if(distr == "Lomax"){
    
    #最小値がゼロ未満だとエラー
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }
    
    #対数尤度の計算
    dlomax.ll <- function(x, alpha, lambda){
      ret <- sum(log(dLomax(x = x, alpha = abs(alpha), lambda = abs(lambda))))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dplomax.opt <- function(x){
      ret <- dlomax.ll(x = data, alpha = x[1], lambda = x[2])
      return(ret)
    }
    
    #対数尤度を最大化
    lomax.opt <- optim(par = c(1, 1), 
                     fn = dplomax.opt, control = list(fnscale  = -1))
    
    fitdist.start <- list(alpha = abs(lomax.opt$par[1]), 
                          lambda = abs(lomax.opt$par[2]))
    
    print(fitdist.start)
    
    
    fitdist.lower <- c(0, 0)
  }
  
  #ピアソン　タイプI分布
  if(distr == "pearson1"){
    
    #対数尤度の計算
    dPearson1.ll <- function(x, a, b, location, scale){
      if(min(a, b) <= 0){return(-Inf)}
      ret <- sum(log(dPearson6(x = x, a, b, location, scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dPearson1.opt <- function(x){
      ret <- dPearson1.ll(x = data, a = x[1], b = x[2], location = x[3], scale = x[4])
      return(ret)
    }
    
    #対数尤度を最大化
    Pearson1.opt <- optim(par = c(1, 1, min(data) - 3, 1), 
                          fn = dPearson1.opt, control = list(fnscale  = -1))
    
    fitdist.start <- list(a = Pearson1.opt$par[1], b = Pearson1.opt$par[2], 
                          location = Pearson1.opt$par[3], scale = Pearson1.opt$par[4])

    
    #fitdist.start <- list(a = 5, b = 1, location = min(data) + 3, scale = 1)
    fitdist.lower <- c(0, 0, -Inf, -Inf)
  }
  
  #ピアソン　タイプII分布
  if(distr == "pearson2"){
    fitdist.start <- list(a = 5, location = 1, scale = 1)
    fitdist.lower <- c(0, -Inf, -Inf)
  }
  
  #ピアソン　タイプIII分布
  if(distr == "pearson3"){
    
    #shape, location, scale
    #s<>0, a>0 and (x-lambda)/s>=0., a = shape, lambda = location, s = scale

    #対数尤度の計算
    dpearson3.ll <- function(x, shape, location, scale){
      
      if(scale == 0){return(-Inf)}
      if(shape <= 0){return(-Inf)}
      if((min(data) - location)/scale < 0){return(-Inf)}

      ret <- sum(log(dpearson3(x = x, shape = shape, location = location, scale = scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dpearson3.opt <- function(x){
      ret <- dpearson3.ll(x = data, shape = x[1], location = x[2], scale = x[3])
      return(ret)
    }
    
    #対数尤度を最大化
    pearson3.opt <- optim(par = c(3, min(data) - 3, 3), 
                        fn = dpearson3.opt, control = list(fnscale  = -1))
    
    
    fitdist.start <- list(
      shape = pearson3.opt$par[1], 
      location = pearson3.opt$par[2], 
      scale = pearson3.opt$par[3])

    #fitdist.start <- list(shape = 1, location = mean(data), scale = 1)
    fitdist.lower <- c(1e-10, -Inf, -Inf)
  }

  #ピアソン　タイプIV(4)分布
  if(distr == "pearson4"){
    fitdist.start <- list(m = 5, nu = 1, location = mean(data), scale = 1)
    fitdist.lower <- c(1+1e-10, -Inf, -Inf, -Inf) #m > 1
  }
  
  #ピアソン　タイプV分布
  if(distr == "pearson5"){
    fitdist.start <- list(shape = 1, location = mean(data), scale = 1)
    fitdist.lower <- c(1e-10, -Inf, -Inf)
  }
  
  #ピアソン　タイプVI(6)分布
  if(distr == "Pearson6"){
    #Pearsonが大文字で始まることに注意。actuarパッケージのpearson6と重複するため。
    
    #対数尤度の計算
    dPearson6.ll <- function(x, a, b, location, scale){
      if(min(a, b) <= 0){return(-Inf)}
      ret <- sum(log(dPearson6(x = x, a, b, location, scale)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dPearson6.opt <- function(x){
      ret <- dPearson6.ll(x = data, a = x[1], b = x[2], location = x[3], scale = x[4])
      return(ret)
    }
    
    #対数尤度を最大化
    Pearson6.opt <- optim(par = c(1, 1, min(data) - 3, 1), 
                        fn = dPearson6.opt, control = list(fnscale  = -1))

    fitdist.start <- list(a = Pearson6.opt$par[1], b = Pearson6.opt$par[2], 
                          location = Pearson6.opt$par[3], scale = Pearson6.opt$par[4])
    
    fitdist.lower <- c(1e-10, 0.1, -Inf, -Inf)
  }
  
  #ピアソン　タイプVII分布
  if(distr == "pearson7"){
    fitdist.start <- list(df = 1, location = mean(data), scale = 1)
    fitdist.lower <- c(0.1, -Inf, -Inf)
  }
  
  #Burr分布
  if(distr == "Burr"){
    
    #最小値がゼロ以下だとエラー
    if(min(data) <= 0){
      return(error.ret(Sys.time()))
    }
    
    #初期値
    fitdist.start <- list(k = 1, c = 1)
    
    #最小値（ゼロより少し大きな値）
    fitdist.lower <- c(0.01, 0.01)
  }
  
  #Extended Burr type XII
  if(distr == "BurrXII"){
    
    #最小値がゼロ未満だとエラー
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }
    
    
    #L-momentsから初期値を計算
    start.burr <- function(data){
      
      #L-momentsを計算
      data.lmom <- lmom.ub(data)
      
      #初期値を計算
      start.burr <- parBurrXII.approx(
        L1 = data.lmom$L1, tau = data.lmom$LCV, tau3 = data.lmom$TAU3)
      
      #返り値
      return(start.burr)
    }
    
    #初期値計算
    start.b <- try(start.burr(data), silent = FALSE)
    
    #エラー処理
    if(class(start.b)[1] == "try-error" || start.b %>% na.omit() %>% length() < 3){
      fitdist.start <- list(lambda = 1, k = -1, c = 1)
    }else{
      #初期値ベクトルについている名前を消す
      names(start.b) <- NULL
      
      #初期値に代入
      fitdist.start <- list(
        lambda = start.b[1], 
        k = start.b[2], 
        c = start.b[3]
      )
    }
    
    #最小値と最大値
    fitdist.lower <- c(0, -Inf, 0)
    fitdist.upper <- c(Inf, 0, Inf)
  }
  
  
  #Johnson SU分布
  if(distr == "johnsonSU"){
   
    #パラメータ推定
    su.param <- try.null(eJohnsonSU(sample.data))
    
    if(is.null(su.param)){
      fitdist.start <- list(gamma = 1, 
                            delta = 1, 
                            xi = 1, 
                            lambda = 1)
      
    }else{
      fitdist.start <- list(gamma = su.param$gamma, 
                            delta = su.param$delta, 
                            xi = su.param$xi, 
                            lambda = su.param$lambda
      )
      
    }

    
    fitdist.lower <- c(-Inf, 1e-10, -Inf, 1e-10)
  }
  
  #Johnson SB分布
  if(distr == "johnsonSB"){
    
    #http://www.ntrand.com/jp/johnson-sb-distribution/
    #xiとlambdaの値を設定
    
    #パラメータ推定
    sb.param <- try.null(eJohnsonSB(data))
    

    if(is.null(sb.param)){
      #パラメータ推定に失敗した場合は初期値に全部1をとりあえず入れる
      fitdist.start <- list(gamma = 1, 
                            delta = 1, 
                            xi = 1, 
                            lambda = 1)
      
      fitdist.lower <- c(-Inf, 1e-10, -Inf, 1e-10)
      
    }else{
      
      #初期値にそのまま推定値を入れる
      fitdist.start <- list(gamma = sb.param$gamma, 
                            delta = sb.param$delta, 
                            xi = sb.param$xi, 
                            lambda = sb.param$lambda
                            )
      
      #下限と上限を推定値の近傍で設定する。（±0.01)
      start.vec <- c(sb.param$gamma, sb.param$delta, sb.param$xi, sb.param$lambda)
      fitdist.lower <- start.vec - 0.01
      fitdist.upper <- start.vec + 0.01
      
    }
    
  }
  
  
  #バーンバウム　サンダース分布の初期値
  if(distr == "fatigue"){
    fitdist.start <- list(alpha = 0.5, beta = 1, mu = 0)
    fitdist.lower <- c(0, 0, -Inf)
  }
  
  #ラプラス分布の初期値
  if(distr == "Laplace"){
    fitdist.start <- list(mu = 0, sigma = 1)
    fitdist.lower <- c(-Inf, 0)
  }
  
  #gompertz分布の初期値
  if(distr == "gompertz"){
    fitdist.start <- list(a = 1, b = 1)
    fitdist.lower <- c(1e-10, 1e-10)
  }
  
  #Muth分布の初期値
  if(distr == "muth"){
    
    #最小値がゼロより小さい場合は空の結果を返す
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(alpha = 0.5)
    fitdist.lower <- c(1e-10)
    fitdist.upper <- c(1)
  }
  
  #ロジスティック分布の初期値
  if(distr == "llogis"){
    fitdist.lower <- c(0, 0)
  }

  #双曲線正割分布の初期値
  if(distr == "hs"){
    fitdist.start <- list(mu = 0, sigma = 1)
    fitdist.lower <- c(-Inf, 1e-10)
  }
  
  #アーラン分布の初期値
  if(distr == "erlang"){
    
    #最小値がゼロより小さい場合は空の結果を返す
    if(min(data) < 0){
      return(error.ret(Sys.time()))
    }

    #対数尤度の計算
    derlang.ll <- function(x, k, mu){
      if(min(k, mu) <= 0){return(-Inf)}
      ret <- sum(log(derlang(x = x, k = k, mu = mu)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    derlang.opt <- function(x){
      ret <- derlang.ll(x = data, k = x[1], mu = x[2])
      return(ret)
    }
    
    #対数尤度を最大化
    erlang.opt <- optim(par = c(1, 1), 
                        fn = derlang.opt, control = list(fnscale  = -1))
    
    #kの最適値と推測した値
    k.int <- round(erlang.opt$par[1])
    

    fitdist.start <- list(k = k.int, mu = erlang.opt$par[2])
    
    fitdist.lower <- c(k.int-1e-10, 1e-10)
    fitdist.upper <- c(k.int+1e-10, Inf)

  }
  
  #Voigt分布の初期値
  if(distr == "voigt"){
    
    fitdist.start <- list(
      x0 = mean(data),
      sigma = sd(data)/2,
      gamma = sd(data)/2
    )
    
    fitdist.lower <- c(-Inf, 0, 0)
    
  }
  
  #レヴィ分布
  if(distr == "Levy"){
    
    fitdist.start <- list(
      m = min(data) - 1,
      s = 1
    )
    
    fitdist.lower <- c(-Inf, 0)
    fitdist.upper <- c(min(data), Inf)
    
  }
  
  #ベータ分布の初期値
  if(distr == "beta"){
    
    
    #[0, 1]に入らない場合は空の結果を返す
    if(min(data) < 0 || max(data) > 1){
      return(error.ret(Sys.time()))
    }
    
    
    
    fitdist.start <- list(shape1 = 1, shape2 = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #カイ二乗分布
  if(distr == "chi2"){
    fitdist.start <- list(df = 1)
    fitdist.lower <- c(1e-10)
  }
  
  #非心カイ二乗分布
  if(distr == "ncchi2"){
    fitdist.start <- list(df = 1, ncp = 1)
    fitdist.lower <- c(1e-10, 0)
  }
  
  #t分布の場合の初期値
  if(distr == "t2"){
    fitdist.lower <- c(1e-10)
    fitdist.start <- list(df = 1)
  }
  
  #非心t分布の場合の初期値
  if(distr == "nct"){
    
    #対数尤度の計算
    dnct.ll <- function(x, df, ncp){
      if(df <= 0){return(-Inf)}
      ret <- sum(log(dnct(x = x, df = df, ncp = ncp)))
      if(is.nan(ret)){ret <- -Inf}
      return(ret)
    }
    
    #パラメータから対数尤度を求める関数
    dnct.opt <- function(x){
      ret <- dnct.ll(x = data, df = x[1], ncp = x[2])
      return(ret)
    }
    
    #対数尤度を最大化
    nct.opt <- optim(par = c(2, mean(data)), 
                        fn = dnct.opt, control = list(fnscale  = -1))
    

    fitdist.start <- list(df = nct.opt$par[1], ncp = nct.opt$par[2])
    
    
    fitdist.lower <- c(1e-10, -Inf)
  }
  
  #F分布の場合の初期値
  if(distr == "F"){
    fitdist.start <- list(df1 = 1, df2 = 1)
    fitdist.lower <- c(1e-10, 1e-10)
  }
  
  #非心F分布の場合の初期値
  if(distr == "ncF"){
    fitdist.start <- list(df1 = 1, df2 = 1, ncp = 0)
    fitdist.lower <- c(1e-10, 1e-10, -Inf)
  }
  

  #2変量混合正規分布の場合の初期値
  if(distr == "normmix2"){
    
    #EMアルゴリズム
    normalmixEM.res <- try(normalmixEM(data, k = 2), silent = TRUE)
    
    #エラーの場合は止める
    if(class(normalmixEM.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    
    fitdist.start <- list(
      mean1 = normalmixEM.res$mu[1], 
      sd1 = normalmixEM.res$sigma[1],
      rate1 = normalmixEM.res$lambda[1],
      
      mean2 = normalmixEM.res$mu[2],
      sd2 = normalmixEM.res$sigma[2],
      rate2 = normalmixEM.res$lambda[2])
    
    
    fitdist.lower <- c(-Inf, 0, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, 1, Inf, Inf, 1)
    
    
  }
  
  #3変量混合正規分布の場合の初期値
  if(distr == "normmix3"){
    
    #EMアルゴリズム
    normalmixEM.res <- try(normalmixEM(data, k = 3), silent = TRUE)
    
    #エラーの場合は止める
    if(class(normalmixEM.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    
    fitdist.start <- list(
      mean1 = normalmixEM.res$mu[1], 
      sd1 = normalmixEM.res$sigma[1],
      rate1 = normalmixEM.res$lambda[1],
      
      mean2 = normalmixEM.res$mu[2],
      sd2 = normalmixEM.res$sigma[2],
      rate2 = normalmixEM.res$lambda[2],
      
      mean3 = normalmixEM.res$mu[3],
      sd3 = normalmixEM.res$sigma[3],
      rate3 = normalmixEM.res$lambda[3])
    
    
    fitdist.lower <- c(-Inf, 0, 0, -Inf, 0, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, 1, Inf, Inf, 1, Inf, Inf, 1)
    
    
  }
  
  #4変量混合正規分布の場合の初期値
  if(distr == "normmix4"){
    
    #EMアルゴリズム
    normalmixEM.res <- try(normalmixEM(data, k = 4), silent = TRUE)
    
    #エラーの場合は止める
    if(class(normalmixEM.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(
      mean1 = normalmixEM.res$mu[1], 
      sd1 = normalmixEM.res$sigma[1],
      rate1 = normalmixEM.res$lambda[1],
      
      mean2 = normalmixEM.res$mu[2],
      sd2 = normalmixEM.res$sigma[2],
      rate2 = normalmixEM.res$lambda[2],
      
      mean3 = normalmixEM.res$mu[3],
      sd3 = normalmixEM.res$sigma[3],
      rate3 = normalmixEM.res$lambda[3],
      
      mean4 = normalmixEM.res$mu[4],
      sd4 = normalmixEM.res$sigma[4],
      rate4 = normalmixEM.res$lambda[4]
      
      )
    
    
    fitdist.lower <- c(-Inf, 0, 0, -Inf, 0, 0, -Inf, 0, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, 1, Inf, Inf, 1, Inf, Inf, 1, Inf, Inf, 1)
    
    
  } 
  
  
  #2変量混合対数正規分布の場合の初期値
  if(distr == "lnormmix2"){
    
    #EMアルゴリズム
    normalmixEM.res <- try(normalmixEM(log(data), k = 2), silent = TRUE)
    
    #エラーの場合は止める
    if(class(normalmixEM.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    
    fitdist.start <- list(
      meanlog1 = normalmixEM.res$mu[1], 
      sdlog1 = normalmixEM.res$sigma[1],
      rate1 = normalmixEM.res$lambda[1],
      
      meanlog2 = normalmixEM.res$mu[2],
      sdlog2 = normalmixEM.res$sigma[2],
      rate2 = normalmixEM.res$lambda[2])
    
    
    fitdist.lower <- c(-Inf, 0, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, 1, Inf, Inf, 1)
    
    
  }
  
  #切断正規分布の場合の初期値
  if(distr == "tnorm"){
    
    fitdist.start <- list(mean = mean(data), sd = sd(data), a= min(data) - sd(data), b = max(data) + sd(data))
    fitdist.lower <- c(-Inf, 0, -Inf, -Inf)

    
    
  } 
  
  #切断正規分布の場合の初期値
  if(distr == "zmnorm"){
    
    fitdist.start <- list(mean = mean(data), sd = sd(data), p.zero = 0.5)
    fitdist.lower <- c(-Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf,1)
 
  } 
  
  
  #ポアソン分布の場合
  if(distr == "pois"){
    #整数ではないとだとエラー
    if(!is_integer(data) || min(data) < 0){
      return(error.ret(Sys.time()))
    }
  }

  #幾何分布の場合
  if(distr == "geom"){
    #整数ではないとだとエラー
    if(!is_integer(data) || min(data) < 0){
      return(error.ret(Sys.time()))
    }
  }
  
  #負の二項分布の場合
  if(distr == "nbinom"){
    #整数ではないとだとエラー
    if(!is_integer(data) || min(data) < 0){
      return(error.ret(Sys.time()))
    }
  }
  
  #三角分布の場合
  if(distr == "triangle"){
    fitdist.start <- list(a = min(data)-sd(data), b = max(data) + sd(data), c = (min(data) + max(data))/2)
    
  }
  
  
  
  #計算中の分布関数を表示
  print(distr)
  
  #フィッティングを関数に
  fitdist.fn <- function(method = method){
    
    #フィッティング
    ret <- suppressWarnings(
      try(
        withTimeout({
          fitdist(data = data, distr = distr, method = method, 
                  start = fitdist.start, fix.arg = fix.arg,
                  lower = fitdist.lower, upper = fitdist.upper,
                  optim.method = optim.method)
          
        }, timeout = timeout, onTimeout = "error")
        
        
        , silent = FALSE)
    )
    
    
    #もし、対数尤度が閾値を超えたらエラーとする
    if(class(ret)[1] == "fitdist"){
      
      if(is.numeric(ret$loglik)){
        
        if(ret$loglik >= loglik.th){
          
          
          ret <- list()
          class(ret) <- "try-error"
        }
        
      }
      
      
    }
    

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
  
  #閾値の以下の値をNAに変換
  th_na <- function(vec, th = -Inf){
    
    #戻り値にまず初期値を入れる
    ret <- vec
    
    #等号を入れたのはth = -Infのとき、-Infのみを検出するため。
    ret[vec <= th] <- NA
    
    #戻り値
    return(ret)
    
  }
  
  
  #AICでおかしいと判断するAICの値
  aic.th <- -Inf
  
  #空のtibble
  df0 <- tibble()

  #各結果をデータフレームに変換
  for(i in 1:length(data)){

    print(names(data)[i][1])
    
    #適合度の統計量計算
    gofstat.res <- try(gofstat(data[[i]]), silent = FALSE)
    
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
        "Chi-squared" = null.na(gofstat.res$chisq)[1],
        "Chi-squared p-value" = null.na(gofstat.res$chisqpvalue)[1],
        
        #計算時間
        "Calculation time" = data[[i]]$CalculationTime[1],
        
        #解の存在
        "Solution" = TRUE
      )

    }else{
      #エラーの場合
      
      #データフレーム作成
      df0 <- tibble(
        distr = names(data)[i][1],
        name = data[[i]]$name[1],
        Solution = FALSE
        )

    }

    #結合
    if(!is.null(ret)){
      ret <- dplyr::bind_rows(ret, df0)
    }else{
      ret <- df0
    }

  }
  
  #並び替え用のAIC
  ret <- ret %>% dplyr::mutate(AIC2 = th_na(AIC, th = aic.th))
  
  #並び替え用AICの小さな順に並び替え
  ret <- ret %>% dplyr::arrange(AIC2)
  
  #並び替え用AICを削除
  ret[, c("AIC2")] <- NULL
  
  
  #重複を削除
  ret <- dplyr::distinct(ret)

  #戻り値
  return(ret)
}

#文字列末尾確認
chk.end <- function(chk, chr){
  
  chr.pos <- grep(paste0("\\", chk, "$"), chr)
  if(length(chr.pos) == 0){return(FALSE)}else{return(TRUE)}
  
}

#data.frameの行数を拡張
df.long <- function(df, length){
  
  #元のデータフレームの長さ
  df.nrow <- nrow(df)
  df.ncol <- ncol(df)
  
  #元の長さがlengthより長ければエラー
  if(df.ncol > length){stop("length is too small.")}
  
  #長さが同じならそのまま返す
  if(df.ncol == length){return(df)}
  
  
  #不足しているデータフレームを作成
  add.df.nrow <- length - df.nrow
  
  add.df <- matrix(rep(NA, add.df.nrow * df.ncol), ncol = df.ncol) %>% 
    as.data.frame() %>% as_tibble()
  colnames(add.df) <- colnames(df)
  
  #結合
  ret <- df %>% rbind(add.df)
  
  #戻り値
  return(ret)
  
}

#データ読み込み
read.data <- function(file){

  #csvの場合
  if(chk.end(".csv", file)){
    
    #csvが数値だけのデータの場合はシート全体をベクトル化
    #1行目からデータを読み込む
    ret.chk.raw <- suppressMessages(read_csv(file, col_names = FALSE))  
    
    #数値のみ選択
    ret.chk <- ret.chk.raw %>% select_if(is.numeric) 
    
    #文字の列が全くない場合
    if(ncol(ret.chk.raw) == ncol(ret.chk)){
      
      #ベクトル化
      ret <- tibble(Data = ret.chk %>% as.vec())
      
      #戻り値
      return(ret)
      
    }
    
    #文字の列が存在する阿合は、項目名から再度読み込んで各列を数値化
    ret <- suppressMessages(read_csv(file))   %>%
      select_if(is.numeric) #数値のみ選択
    
    #戻り値
    return(ret)
  }

  #エクセル形式の場合
  if(chk.end(".xlsx", file) || chk.end(".xls", file) || chk.end(".xlsm", file)){

    #空のデータフレームを準備
    ret <- tibble()
    
    #シート名を抜き出す
    excel.sheets <- excel_sheets(file)

    #各シートから抜き出す
    for(i in 1:length(excel.sheets)){
      
      #シートがすべて数値になっているかチェック
      df.i.chk.raw <- try(read_excel(file, sheet = i, col_names = FALSE), 
                          silent = FALSE) #項目名はなし
      df.i.chk <- try(df.i.chk.raw  %>% select_if(is.numeric), silent = FALSE)
      df.i.chk.error <- try(ncol(df.i.chk.raw) == ncol(df.i.chk), silent =FALSE)
      
      #数値のみ抜き出したデータが、元のデータ数と一致するかチェック
      if(df.i.chk.error == TRUE){
        
        #一致する場合
        df.i <- data.frame(Data = df.i.chk %>% as.vec())
        colnames(df.i) <- excel.sheets[i]
        
      }else{
        
        #シートから読み込みし数値列のみ抜き出す
        df.i <- try(read_excel(file, sheet = i) %>% select_if(is.numeric), silent = FALSE)
        
        
      }
      


      #エラーが起きない場合の処理
      if(class(df.i)[1] != "try-error"){
        
        #データが存在する場合
        if(nrow(df.i) > 0 && ncol(df.i) > 0){
          
          #戻り値に何かデータが入っているかチェック
          if(nrow(ret) > 0){
            
            #入っている場合
            
            #長い方の行数
            max.nrow <- max(nrow(df.i), nrow(ret))
            
            #行数を拡張したうえで結合
            ret <- tibble(df.long(ret, max.nrow), df.long(df.i, max.nrow))
            
          }else{
            
            #入ってない場合
            ret <- df.i
          }
          
        }
      }
    }
    
    #結合したデータを返す
    return(ret)

  }
  
  #該当がなかった場合、NAを返す
  return(NULL)
}

#確率紙にプロット
plot_paper <- function(result, rank = "median", method = "norm"){
  
  #参考
  #http://aoki2.si.gunma-u.ac.jp/R/npp2.html
  #http://aoki2.si.gunma-u.ac.jp/R/weibull.html
  
  #エラーチェック
  if(is.null(result)){return(result)}
  if(class(result)[1] != "fitdist"){return(NULL)}
  if(is.null(rank) || is.null(method)){return(NULL)}
  
  #パラメータ推定に失敗した場合
  if(result$estimate %>% na.omit() %>% length() == 0){
    return(NULL)
  } 
  
  #小さい順に並んだデータを抜き出す
  data <- result$data %>% sort
  
  #メジアンランクの場合
  if(rank == "median"){
    mol_plus <- -0.3  #分子
    den_plus <- 0.4   #分母
  }
  
  #平均ランクの場合
  if(rank == "mean"){
    mol_plus <- 0  #分子
    den_plus <- 1   #分母
  }
  
  #不信頼度
  rank_i <- c(1:length(data))
  rank_n <- length(data)
  fi <- (rank_i + mol_plus)/(rank_n + den_plus)
  
  #縦軸のFのベクトルを作成する関数
  make.f.vec <- function(f.min, length = NULL){
    
    #スケールケーズの最小値の自然対数
    f.min.log <- floor(log10(f.min))
    
    #スケールに使用する数値作成
    if(is.null(length)){
      f.vec.small <- 10^c(f.min.log : -1)
    }else{
      f.vec.small <- 10^seq(f.min.log, -1, length = length)
    }
    
    
    f.vec.large <- 1 - f.vec.small
    f.vec.std <- c(0.2, 0.5, 0.632)
    f.vec <- sort(c(f.vec.small, f.vec.std, f.vec.large))
    
    #戻り値
    return(f.vec)
  }
  
  #縦軸のベクトルを作成
  probs <- make.f.vec(min(fi))
  
  #ワイブル分布の尺度に変更
  weib <- function(p) log10(log10(1/(1-p)))
  
  #指数確率紙の尺度に変更
  #http://www.fml.t.u-tokyo.ac.jp/~sakai/kougi/ProbSystem/ProbPaper/probpaper.htm
  expp <- function(p) log(1/(1-p))
  
  
  #対数正規確率かワイブルの場合(x軸が対数)
  if(method == "lnorm" || method == "weibull"){
    
    #対数スケールだと0以下はプロットできないので、NULLを返す
    if(min(data) <= 0){
      return(NULL)
    }
    
    #x軸を対数スケール
    plot.log <- "x"
    
    #q.vecを作る
    q.vec <- exp(seq(log(min(data)) - 1, log(max(data)) + 1, length = 1000))
    
  }else{
    #x軸をリニア
    plot.log <- ""
    
    #q.vecを作る
    q.vec <- seq(min(data) - 3*sd(data), max(data) + 3*sd(data), length = 1000)
  }
  
  #fitdistの結果とqベクトルから、pベクトルを作る関数定義
  pfit <- function(result, q){
    
    #エラーチェック
    if(class(result)[1] != "fitdist"){return(NULL)}
    
    
    #関数を作成
    p.func <- paste0("p", result$distname)
    
    #最適値を表示
    est.vec <- result$estimate
    
    #パラメータを作る
    est.param <- NULL
    for(i in 1:length(est.vec)){
      #iのパラメータを追加
      est.param <- paste0(est.param, names(est.vec)[i], "=", est.vec[i])
      
      #最後でなければコンマを加える
      if(i < length(est.vec)){est.param <- paste0(est.param, ", ")}
    }
    
    
    eval(parse(text = paste0(
      "ret <- ", p.func, "(q, ", est.param, ")"
    )))
    
    return(ret)
    
  }
  
  #pベクトルを作る
  p.vec <- pfit(result, q.vec)
  
  #正規確率か対数正規確率の場合
  if(method == "norm" || method == "lnorm"){
    plot.y <- qnorm(c(min(probs), max(probs)))
    point.y <- qnorm(fi)
    axis.y <- qnorm(probs)
    p.vec.y <- qnorm(p.vec)
  }
  
  #ワイブルかガンベルの場合
  if(method == "weibull" || method == "gumbel"){
    plot.y <- weib(c(min(probs), max(probs)))
    point.y <- weib(fi)
    axis.y <- weib(probs)
    p.vec.y <- weib(p.vec)
  }
  
  #指数の場合
  if(method == "exp"){
    plot.y <- expp(c(min(probs), max(probs)))
    point.y <- expp(fi)
    axis.y <- expp(probs)
    p.vec.y <- expp(p.vec)
  }
  
  #累積確率分布の場合
  if(method == "cdf"){

    plot.y <- (c(0, 1))
    point.y <- (fi)
    
    axis.y <- seq(0, 1, by = 0.1)
    p.vec.y <- (p.vec)
    
  }
  


  #図を作成
  plot(c(data[1], data[rank_n]), plot.y, 
       log = plot.log,
       type = "n", yaxt = "n",
       xlab = "Data", ylab = "Probability")
  
  #データ点を打つ
  points(data, point.y)
  
  #縦軸の確率を表示
  axis(2, axis.y, probs*100)
  

  #フィッテングした関数の線を引く
  lines(q.vec, p.vec.y, col = "Red")
  
}

#ロジット
logit <- function(x) log(x / (1 - x))

#数値のまとめ表示
vec.summary <- function(vec){
  
  #エラーチェック
  if(is.null(vec)){return(NULL)}
  if(!is.vector(vec)){return(NULL)}
  if(!is.numeric(vec)){return(NULL)}
  
  #NA除去
  vec <- na.omit(vec)
  
  #結果を入れる入れ物
  res <- list()
  
  #結果を格納
  res$n <- length(vec)
  res$Mean <- mean(vec)
  res$SD <- sd(vec)
  res$VAR <- var(vec)
  res$Skewness <- EnvStats::skewness(vec)
  res$Kurtosis <- EnvStats::kurtosis(vec, excess = FALSE)
  res$ExcessKurtosis <- EnvStats::kurtosis(vec, excess = TRUE) #過剰尖度
  res$Median <- median(vec)
  res$Max <- max(vec)
  res$Min <- min(vec)
  
  #結果を結合
  ret <- ""
  
  for(i in 1:length(res)){
    ret0 <- paste(names(res)[i], "=", 
                  signif(res[[i]], digits = 4),
                  ", ")
    ret <- paste0(ret, ret0)
  }
  
  #戻り値
  return(ret)
}

#統計量とパラメータ表示
fitdist_summary <- function(result){
  
  #エラーチェック
  if(is.null(result)){return(NULL)}
  if(class(result)[1] != "fitdist"){return(NULL)}
  
  #パラメータ推定に失敗した場合
  if(result$estimate %>% na.omit() %>% length() == 0){
    return("Parameter estimation failed.")
  } 
  
  #適合度の統計量計算
  gofstat.res <- try(gofstat(result), silent = TRUE)
  
  #パラメータを整理
  ret1 <- "Parameters :"
  for(i in 1:length(result$estimate)){
    ret1 <- paste0(ret1, "\n", names(result$estimate)[i], " = ", result$estimate[i] %>% signif(4))
  }
  ret1 <- paste0(ret1, "\n")
  
  #retベクトルにパラメータと統計量を格納
  ret <- NULL
  ret[1] <- ret1[1]
  
  ret[2] <- paste0("n = ", result$n, ", ",
                   "AIC = ", null.na(result$aic)[1] %>% round(2),  ", ",
                   "BIC = ", null.na(result$bic)[1] %>% round(2), ", ",
                   "log likelihood = ", null.na(result$loglik)[1] %>% round(2))
  
  ret[3] <- paste0("Kolmogorov-Smirnov, ",
                   "D = ", null.na(gofstat.res$ks)[1] %>% signif(4), ", ",
                   "p-value = ", kspval(result$n, gofstat.res$ks)[1] %>% signif(4))
  
  ret[4] <- paste0("Cramer-von Mises, ",
                   "omega2 = ", null.na(gofstat.res$cvm)[1] %>% signif(4), ", ",
                   "p-value = ", cvmpval(result$n, gofstat.res$cvm)[1] %>% signif(4))
  
  ret[5] <- paste0("Anderson-Darling, ",
                   "An = ", null.na(gofstat.res$ad)[1] %>% signif(4), ", ",
                   "p-value = ", adpval(result$n, gofstat.res$ad)[1] %>% signif(4))
  
  ret[6] <- paste0("Chi-squared = ",null.na(gofstat.res$chisq)[1] %>% signif(4), ", ",
                   "p-value = ", null.na(gofstat.res$chisqpvalue)[1] %>% signif(4))
  
  
  #改行で結合
  ret <- paste0(ret, sep = "\n")
  
  #cat(ret)
  return(ret)
}

#qから累積確率を計算
pdist <- function(result, q){
  
  #エラーチェック
  if(is.null(result)){return(result)}
  if(class(result)[1] != "fitdist"){return(NULL)}
  if(is.null(q)){return(NULL)}
  
  #パラメータ推定に失敗した場合
  if(result$estimate %>% na.omit() %>% length() == 0){
    return(NULL)
  } 
  
  #fitdistの結果とqベクトルから、pベクトルを作る関数定義
  pfit <- function(result, q){
    
    #エラーチェック
    if(class(result)[1] != "fitdist"){return(NULL)}
    
    
    #関数を作成
    p.func <- paste0("p", result$distname)
    
    #最適値を表示
    est.vec <- result$estimate
    
    #パラメータを作る
    est.param <- NULL
    for(i in 1:length(est.vec)){
      #iのパラメータを追加
      est.param <- paste0(est.param, names(est.vec)[i], "=", est.vec[i])
      
      #最後でなければコンマを加える
      if(i < length(est.vec)){est.param <- paste0(est.param, ", ")}
    }
    
    
    eval(parse(text = paste0(
      "ret <- ", p.func, "(q, ", est.param, ")"
    )))
    
    return(ret)
    
  }
  
  #pベクトルを作る
  p.vec <- pfit(result, q)
  
}

#qから累積確率を計算
qdist <- function(result, p){
  
  #エラーチェック
  if(is.null(result)){return(result)}
  if(class(result)[1] != "fitdist"){return(NULL)}
  if(is.null(p)){return(NULL)}
  
  #パラメータ推定に失敗した場合
  if(result$estimate %>% na.omit() %>% length() == 0){
    return(NULL)
  } 
  
  #fitdistの結果とpベクトルから、qベクトルを作る関数定義
  qfit <- function(result, p){
    
    #エラーチェック
    if(class(result)[1] != "fitdist"){return(NULL)}
    
    
    #関数を作成
    q.func <- paste0("q", result$distname)
    
    #最適値を表示
    est.vec <- result$estimate
    
    #パラメータを作る
    est.param <- NULL
    for(i in 1:length(est.vec)){
      #iのパラメータを追加
      est.param <- paste0(est.param, names(est.vec)[i], "=", est.vec[i])
      
      #最後でなければコンマを加える
      if(i < length(est.vec)){est.param <- paste0(est.param, ", ")}
    }
    
    
    eval(parse(text = paste0(
      "ret <- ", q.func, "(p, ", est.param, ")"
    )))
    
    return(ret)
    
  }
  
  #qベクトルを作る
  q.vec <- qfit(result, p)
  
}

#chrを先頭につけてnumの数値を表示
chr.num <- function(num, chr){
  
  if(is.null(chr) || is.null(num)){return(NULL)}
  
  ret <- paste0(chr, num)
  
  return(ret)
  
}

#エラーの場合NULLを返す
try.null <- function(res){
  
  ret <- try(res, silent = TRUE)
  if(class(ret)[1] == "try-error"){return(NULL)}
  return(ret)
  
}

#データフレームでなかったらNULLを返す
is.data.frame.null <- function(obj){
  
  if(is.data.frame(obj)){
    return(obj)
  }else{
    return(NULL)
  }
  
  
  
  
}


