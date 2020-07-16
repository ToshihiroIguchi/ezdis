library(EnvStats)

#EnvStatsパッケージのparetoを別名にコピー
dPareto <- function(x, shape, location){EnvStats::dpareto(x, shape, location)}
pPareto <- function(q, shape, location){EnvStats::ppareto(q, shape, location)}
qPareto <- function(p, shape, location){EnvStats::qpareto(p, shape, location)}


