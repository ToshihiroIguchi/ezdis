
#aka Beta
#Pearson type Iの確率密度関数
dpearson1 <- function(x, a, b, location, scale){
  ret <- dpearsonI(x = x, a = a, b = b, location = location, scale = scale)
  return(ret)
}

#Pearson type Iの累積分布関数
ppearson1 <- function(q, a, b, location, scale){
  ret <- ppearsonI(q = q, a = a, b = b, location = location, scale = scale)
  return(ret)
}

#Pearson type Iの累積分布関数の逆関数
qpearson1 <- function(p, a, b, location, scale){
  ret <- qpearsonI(p = p, a = a, b = b, location = location, scale = scale)
  return(ret)
}


#aka Symmetric Beta
#Pearson type IIの確率密度関数
dpearson2 <- function(x, a, location, scale){
  ret <- dpearsonII(x = x, a = a, location = location, scale = scale)
  return(ret)
}

#Pearson type IIの累積分布関数
ppearson2 <- function(q, a, location, scale){
  ret <- ppearsonII(q = q, a = a, location = location, scale = scale)
  return(ret)
}

#Pearson type IIの累積分布関数の逆関数
qpearson2 <- function(p, a, location, scale){
  ret <- qpearsonII(p = p, a = a, location = location, scale = scale)
  return(ret)
}


#aka Gamma
#Pearson type IIIの確率密度関数
dpearson3 <- function(x, shape, location, scale){
  ret <- dpearsonIII(x = x, shape = shape, location = location, scale = scale)
  return(ret)
}

#Pearson type IIIの累積分布関数
ppearson3 <- function(q, shape, location, scale){
  ret <- ppearsonIII(q = q, shape = shape, location = location, scale = scale)
  return(ret)
}

#Pearson type IIIの累積分布関数の逆関数
qpearson3 <- function(p, shape, location, scale){
  ret <- qpearsonIII(p = p, shape = shape, location = location, scale = scale)
  return(ret)
}


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


#aka Inverse Gamma
#Pearson type Vの確率密度関数
dpearson5 <- function(x, shape, location, scale){
  ret <- dpearsonV(x = x, shape = shape, location = location, scale = scale)
  return(ret)
}

#Pearson type Vの累積分布関数
ppearson5 <- function(q, shape, location, scale){
  ret <- ppearsonV(q = q, shape = shape, location = location, scale = scale)
  return(ret)
}

#Pearson type Vの累積分布関数の逆関数
qpearson5 <- function(p, shape, location, scale){
  ret <- qpearsonV(p = p, shape = shape, location = location, scale = scale)
  return(ret)
}


#aka Beta Prime
#Pearson type VIの確率密度関数
dPearson6 <- function(x, a, b, location, scale){
  ret <- dpearsonVI(x = x, a = a, b = b, location = location, scale = scale)
  return(ret)
}

#Pearson type VIの累積分布関数
pPearson6 <- function(q, a, b, location, scale){
  ret <- ppearsonVI(q = q, a = a, b = b, location = location, scale = scale)
  return(ret)
}

#Pearson type VIの累積分布関数の逆関数
qPearson6 <- function(p, a, b, location, scale){
  ret <- qpearsonVI(p = p, a = a, b = b, location = location, scale = scale)
  return(ret)
}


#aka Student's t
#Pearson type VIIの確率密度関数
dpearson7 <- function(x, df, location, scale){
  ret <- dpearsonVII(x = x, df = df, location = location, scale = scale)
  return(ret)
}

#Pearson type VIIの累積分布関数
ppearson7 <- function(q, df, location, scale){
  ret <- ppearsonVII(q = q, df = df, location = location, scale = scale)
  return(ret)
}

#Pearson type VIIの累積分布関数の逆関数
qpearson7 <- function(p, df, location, scale){
  ret <- qpearsonVII(p = p, df = df, location = location, scale = scale)
  return(ret)
}












