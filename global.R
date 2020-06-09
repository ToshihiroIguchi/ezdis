#ライブラリ読み込み
library(ggplot2)
library(fitdistrplus)
library(ismev)

library(evd) #evdはactuarより先に読み込ませて、evd::dgumbelなどをマスクする

library(actuar)


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


#一般化パレート分布の最尤推定を Sifferらの方法で行う関数
gp.fit_siffer <- function(x, K = 10, epsilon = 10^-8) {
  #https://hoxo-m.hatenablog.com/entry/2018/12/25/220217
  
  #' @param x データを数値ベクトルとして入力
  #' @param K デフォルトは 10
  #' @param epsilon 境界を探索範囲から外すための微小値。デフォルトは 10^-8
  #' @return パラメータの最尤推定値 c(sigma, gamma) 
  
    # 式 (4) の方程式
  u <- function(theta) mean(1 / (1 + theta * x))
  v <- function(theta) 1 + mean(log(1 + theta * x))
  h <- function(theta) u(theta) * v(theta) - 1
  
  # 最適化目的関数
  fn <- function(theta_K, logged=FALSE) {
    if(logged) theta_K <- exp(theta_K)
    sum(vapply(theta_K, h, double(1))^2)
  }
  
  # 解の範囲
  lower_bound <- -1 / max(x)
  upper_bound <- 2 * (mean(x) - min(x)) / min(x)^2
  
  theta_range1 <- c(lower_bound + epsilon, -epsilon)
  theta_range2 <- c(epsilon, upper_bound)
  
  # 初期値の設定
  theta_init1 <- seq(lower_bound, 0, length.out = K+2)[-c(1, K+2)]
  theta_init2 <- seq(log(epsilon), log(upper_bound), length.out = K+2)[-c(1, K+2)]
  
  # 最適化の実行
  opt1 <- optim(theta_init1, fn = fn, method = "L-BFGS-B",
                lower = theta_range1[1], upper = theta_range1[2])
  opt2 <- optim(theta_init2, fn = fn, logged = TRUE, method = "L-BFGS-B",
                lower = log(theta_range2[1]), upper = log(theta_range2[2]))
  
  # 最適解 (推定値の候補)
  opt_theta_K <- unique(c(opt1$par, 0, exp(opt2$par)))
  
  # プロファイル対数尤度
  n <- length(x)
  profileLogLik <- function(theta) {
    -n * log(mean(log(1 + theta * x))/theta) - n * mean(log(1 + theta * x)) - n
  }
  
  # 最適解のうちプロファイル対数尤度が最大のものを推定値とする
  pLLs <- vapply(opt_theta_K, profileLogLik, double(1))
  est_theta <- opt_theta_K[which.max(pLLs)]
  
  # パラメータの引き戻し
  est_gamma <- mean(log(1 + est_theta * x))
  est_sigma <- est_gamma / est_theta
  
  # 結果の返却
  c(est_sigma, est_gamma)
}


#分布関数にフィッティング
fit.dist <- function(data, distr = "norm", method = "mle"){
  
  #エラーチェック
  if(is.null(data)){return(NULL)}
  
  #初期時間
  t0 <- Sys.time()
  
  #初期値
  fitdist.start <- NULL
  
  #各分布関数の初期値を指定
  
  #ガンベル分布の場合の初期値
  if(distr == "gumbel"){
    gum.res <- gum.fit(data)
    fitdist.start <- list(alpha = gum.res$mle[1], scale = gum.res$mle[2])
  }
  
  #一般化極値分布の場合の初期値
  if(distr == "gev"){
    gev.res <- gev.fit(data)
    fitdist.start <- list(loc = gev.res$mle[1], scale = gev.res$mle[2], shape = gev.res$mle[3])
  }

  #一般化極値分布の場合の初期値
  if(distr == "genpareto"){
    gen.pareto.res <- gp.fit_siffer(data)
    fitdist.start <- list(shape1 = gen.pareto.res[1], shape2 = gen.pareto.res[2])
    print(fitdist.start)
  }
  
  
  #計算中の分布関数を表示
  print(distr)
  
  #フィッティング
  ret <- try(fitdist(data = data, distr = distr, method = method, start = fitdist.start), silent = TRUE)
  
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



