#ライブラリ読み込み
library(ggplot2)
library(tibble)

library(readr)
library(readxl)

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

#多変量混合正規分布
source("normmixn.R")

#パレート分布
source("Pareto.R")

#多重モードワイブル分布
source("multiweibull.R")


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
  
  #ワイブル分布の場合の初期値
  if(distr == "weibull"){
    fitdist.lower <- c(0, 0)
  }
  
  #3母数ワイブル分布の初期値
  if(distr == "weibull3"){
    
    fitdist.start <- list(shape = 1, scale = 1, thres = min(data) - 1)
    fitdist.lower <- c(0, 0, -Inf)
    
  }
  
  #逆ワイブル分布の場合の初期値
  if(distr == "rweibull"){
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
  if(distr == "Pareto"){
    fitdist.start <- list(shape = 1, location = 1)
    fitdist.lower <- c(0, 0)
  }
  
  #タイプ2パレート分布の場合の初期値
  if(distr == "pareto2"){
    fitdist.start <- list(min = 1, shape = 1, scale = 1)
    fitdist.lower <- c(-Inf, 0, 0)
  }
  
  #タイプ3パレート分布の場合の初期値
  if(distr == "pareto3"){
    fitdist.start <- list(min = 1, shape = 1, scale = 1)
    fitdist.lower <- c(-Inf, 0, 0)
  }
  
  #タイプ4パレート分布の場合の初期値
  if(distr == "pareto4"){
    fitdist.start <- list(min = 1, shape1 = 1, shape2 = 1, scale = 1)
    fitdist.lower <- c(-Inf, 0, 0, 0)
  }
  
  #バーンバウム　サンダース分布の初期値
  if(distr == "fatigue"){
    fitdist.start <- list(alpha = 0.5, beta = 1, mu = 0)
    fitdist.lower <- c(0, 0, -Inf)
  }
  
  #ラプラス分布の初期値
  if(distr == "laplace"){
    fitdist.start <- list(mu = 0, sigma = 1)
    fitdist.lower <- c(-Inf, 0)
  }
  
  #ロジスティック分布の初期値
  if(distr == "llogis"){
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
    
    #EMアルゴリズム
    normalmixEM.res <- try(normalmixEM(data), silent = TRUE)
    
    #エラーの場合は止める
    if(class(normalmixEM.res)[1] == "try-error"){
      return(error.ret(Sys.time()))
    }
    
    fitdist.start <- list(
      mean1 = normalmixEM.res$mu[1], 
      sd1 = normalmixEM.res$sigma[1],
      mean2 = normalmixEM.res$mu[2],
      sd2 = normalmixEM.res$sigma[2],
      p.mix = 1-normalmixEM.res$lambda[1])
    
    
    fitdist.lower <- c(-Inf, 0, -Inf, 0, 0)
    fitdist.upper <- c(Inf, Inf, Inf, Inf, 1)
    

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
  
  #4変量混合正規分布の場合の初期値
  if(distr == "tnorm"){
    
    fitdist.start <- list(mean = mean(data), sd = sd(data), a= min(data), b = max(data))
    fitdist.lower <- c(-Inf, 0, -Inf, -Inf)

    
    
  } 
  
  
  #一様分布の場合の初期値
  if(distr == "unif"){
    fitdist.start <- list(min = min(data), max = max(data))
    fitdist.lower <- c(-Inf, mean(data))
    fitdist.upper <- c(mean(data),Inf)
    
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
    ret.chk.raw <- read_csv(file, col_names = FALSE)
    
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
    ret <- read_csv(file) %>%
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
  
  #該当がなかった場合、空のデータを返す
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
  
  #fitdistの結果とpベクトルから、qベクトルを作る関数定義
  qfit <- function(result, p){
    
    #エラーチェック
    if(class(result)[1] != "fitdist"){return(NULL)}
    
    #pの値をチェック
    if(min(p) <= 0 || max(p) >=1){
      stop("p must be greater than 0 and less than 1.")
    }
    
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
  
  #pベクトルを作る
  p.vec <- make.f.vec(min(fi), length = 100)
  
  #qベクトルを作る
  q.vec <- qfit(result, p.vec)
  
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
  
  
  #対数正規確率かワイブルの場合(x軸が対数)
  if(method == "lnorm" || method == "weibull"){
    
    #対数スケールだと0以下はプロットできないので、NULLを返す
    if(min(data) <= 0){
      return(NULL)
    }
    
    #x軸を対数スケール
    plot.log <- "x"
  }else{
    #x軸をリニア
    plot.log <- ""
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




