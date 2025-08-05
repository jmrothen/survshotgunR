# Trimmed out portions of code
Rcpp::sourceCpp('dist.cpp')




### BURRT12 / Singh Maddala
# actuar has a burr distrib, but slightly different from burr12?
# VGAM has dsinmad, possible alternative?
# we'll define the p and d ourself here

dBurrT12 <- function(x,c,k,l=1, log=F){
  # c,k = shape, scale=l
  p <- c*k*(1/l)*((x/l)^(c-1)) * ((1 + (x/l)^c)^(-1*k - 1))
  return(
    ifelse(log, log(p), p)
  )
}
dBurrT12 <- Vectorize(dBurrT12)
pBurrT12 <- function(q,c,k,l=1, lower.tail=T,  log.p=F){
  p <- 1-(1+(q/l)^c)^(-k)
  if(!lower.tail){p <- 1-p}
  return(ifelse(log.p,log(p),p))
}
pBurrT12 <- Vectorize(pBurrT12)
flexsurv_burr <- list(
  name='BurrT12',
  pars = c('c', 'k', 'l'),
  location='l',
  transforms = c(log,log,log),
  inv.transforms = c(exp,exp,exp),
  inits= function(t){c(1,2,1)},
  p=pBurrT12,
  d=dBurrT12,
  fullname='burr_type_12'
)

### Erlang 
# custom
derlang_old <- function(x, k, l, log=F){
  f <- ((l^k)*(x^(k-1))*exp(l*-x))/factorial(k-1)
  return(ifelse(log, log(f),f))
}
derlang_old <- Vectorize(derlang_old)
perlang_old <- function(q, k, l, lower.tail=T, log.p=F){
  s <- 0
  for(n in 0:(k-1)){
    s2 <- (1/(factorial(n)))*exp(-l*q)*((l*q)^n)
    s <- s+s2
  }
  p <- 1 - s
  if(!lower.tail){p <- 1-p}
  return(ifelse(log.p, log(p),p))
}
perlang_old <- Vectorize(perlang_old)

flexsurv_erlang <- list(
  name='erlang',
  pars=c('k','l'),
  location='l',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(median(t), 1/median(t))},
  d=derlang,
  p=perlang,
  fullname='erlang'
)



# Gamma Gompertz

dgamgomp_old<- function(x, b, s, beta, log=F){
  f <- (b*s*exp(b*x)*(beta^s))/ ((beta-1+exp(b*x))^(s+1))
  return(ifelse(log, log(f), f))
}
dgamgomp_old <- Vectorize(dgamgomp_old)
pgamgomp_old<- function(q, b, s, beta, lower.tail=T, log.p=F){
  p <- 1 - ((beta^s)/((beta-1+exp(b*q))^s))
  if(!lower.tail){p <- 1-p}
  return(ifelse(log.p,log(p),p))
}
pgamgomp_old <- Vectorize(pgamgomp_old)

# log cauchy
dlogcauchy_old <- function(x, mu, sigma, log=F){
  p <- ((log(x) - mu)^2) + (sigma^2)
  p2 <- (1/(x*pi))*(sigma/p)
  return(ifelse(log, log(p2),p2))
}
dlogcauchy_old <- Vectorize(dlogcauchy_old)
plogcauchy_old <- function(q, mu, sigma, lower.tail=T,log.p=F){
  p <- (1/pi)*atan((log(q)-mu)/sigma) +.5
  if(!lower.tail){p <- 1-p}
  return(ifelse(log.p,log(p),p))
}
plogcauchy_old <- Vectorize(plogcauchy_old)



# repeitive pareto option from extraDistr
 
# flexsurv_pareto <- list(
#   name='pareto',
#   pars=c('a','b'),
#   location='b',
#   transforms=c(log,log),
#   inv.transforms=c(exp,exp),
#   inits=function(t){c(10,10)},
#   d=extraDistr::dpareto,
#   p=extraDistr::ppareto,
#   fullname='pareto'
# )




### Generalized beta (prob too crazy)
# vgam dgenbeta (GB2 specifically aka gen beta prime aka gen F)
# flexsurv_genbeta<- list(
#   name='genbetaII',
#   pars=c(''),
#   location='',
#   transforms=c(),
#   inv.transforms=c(),
#   inits=function(t){c()},
#   d=,
#   p=
# )



library(microbenchmark)

microbenchmark(derlang(1,1,1), derlang_raw(1,1,1), derlang_old(1,1,1))
microbenchmark(perlang(1,1,1), perlang_raw(1,1,1), perlang_old(1,1,1))
microbenchmark(
  derlang(c(1,2,3,4), 5, c(6,7)), 
  # derlang_raw(c(1,2,3,4), 5, c(6,7)), # doesn't work, as not vectorized!
  derlang_old(c(1,2,3,4), 5, c(6,7))
)



microbenchmark(derlang(1,1,1), derlang_old(1,1,1)) 

microbenchmark(perlang_raw(1,1,1),perlang(1,1,1), perlang_old(1,1,1), perlang_v(1,1,1)) 

microbenchmark(dgamgomp(1,.4,1,3), dgamgomp_old(1,.4,1,3)) 

microbenchmark(pgamgomp_raw(1,.4,1,3),pgamgomp(1,.4,1,3), pgamgomp_old(1,.4,1,3), pgamgomp_v(1,.4,1,3))

microbenchmark(dlogcauchy(1,1,1), dlogcauchy_old(1,1,1))

microbenchmark(plogcauchy_raw(1,1,1),plogcauchy(1,1,1), plogcauchy_old(1,1,1), plogcauchy_v(1,1,1))

microbenchmark(pinvlind_v(1,1), pinvlind_raw(1,1), pinvlind(1,1))  
