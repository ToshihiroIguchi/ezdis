
#Pearson type IVの確率密度関数
dpearson4 <- function(x, m, nu, location, scale){
  ret <- dpearsonIV(x = x, m = m, nu = nu, location = location, scale = scale)
  return(ret)
}


#Pearson type IVの累積分布関数
ppearson4 <- function(q, m, nu, location, scale){
  ret <- ppearsonIV(q = q, m = m, nu = nu, location = location, scale = scale)
  return(ret)
}



#Pearson type IVの累積分布関数の逆関数
qpearson4 <- function(p, m, nu, location, scale){
  ret <- qpearsonIV(p = p, m = m, nu = nu, location = location, scale = scale)
  return(ret)
}