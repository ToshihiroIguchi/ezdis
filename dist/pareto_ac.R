#pareto_ac
#パレート分布
library(actuar)

dpareto_ac <- function(x, shape, scale) actuar::dpareto(x = x, shape = shape, scale = scale) 
ppareto_ac <- function(q, shape, scale) actuar::ppareto(q = q, shape = shape, scale = scale) 
qpareto_ac <- function(p, shape, scale) actuar::qpareto(p = p, shape = shape, scale = scale) 

