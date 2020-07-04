#ライブラリ読み込み
library(actuar)

#actuarパッケージのgumbelを別名にコピー
dGumbel <- function(x, alpha, scale) {actuar::dgumbel(x, alpha, scale)}
pGumbel <- function(q, alpha, scale) {actuar::pgumbel(q, alpha, scale)}
qGumbel <- function(p, alpha, scale) {actuar::qgumbel(p, alpha, scale)}

