
#カイ二乗分布
dchi2 <- function(x, df) dchisq(x, df)
pchi2 <- function(q, df) pchisq(q, df)
qchi2 <- function(p, df) qchisq(p, df)


#非心カイ二乗分布
dncchi2 <- function(x, df, ncp) dchisq(x, df, ncp)
pncchi2 <- function(q, df, ncp) pchisq(q, df, ncp)
qncchi2 <- function(p, df, ncp) qchisq(p, df, ncp)

