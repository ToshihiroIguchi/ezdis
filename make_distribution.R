
#乱数種の設定
set.seed(108)

#データ作成
dis.data <- data.frame(
  Normal = rnorm(1000),
  Weibull = rweibull(1000, 2, 1),
  Poisson = rpois(1000, lambda = 1),
  MixNormal = sample(c(rnorm(500, 10,2), rnorm(500, 1, 1)))
)

#データ書き込み
write.csv(dis.data, "sample.csv")




