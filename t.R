
#t分布
dt2 <- function(x, df) dt(x, df, ncp = 0)
pt2 <- function(q, df) pt(q, df, ncp = 0)
qt2 <- function(p, df) qt(p, df, ncp = 0)


#非心t分布
dnct <- function(x, df, ncp) dt(x, df, ncp)
pnct <- function(q, df, ncp) pt(q, df, ncp)
qnct <- function(p, df, ncp) qt(p, df, ncp)

