# setup for stuff to put into basic survival package
mycauchy <- list(name='Cauchy',
                 init= function(x, weights, ...) 
                   c(median(x), mad(x)),
                 density= function(x, parms) {
                   temp <- 1/(1 + x^2)
                   cbind(.5 + atan(x)/pi, .5+ atan(-x)/pi,
                         temp/pi, -2 *x*temp, 2*temp*(4*x^2*temp -1))
                 },
                 quantile= function(p, parms) tan((p-.5)*pi),
                 deviance= function(...) stop('deviance residuals not defined')
)


library(checkmate)

dbetap <- function(x, a,b){
  # a, b > 0 
  return((x^(a-1) * (1+x)^(-1*(a-b)))/beta(a,b))
}



dBirnSaun <- function(x, u, y, b){
  # u < x
  # g, b > 0
  chunk1 <- sqrt((x-u)/b)
  chunk2 <- sqrt(b/(x-u))
  return(
    ((chunk1 + chunk2)/ (2*y*(x-u))) * dnorm((chunk1-chunk2)/y)
  )
}








# sample code

s <- Surv(aml$time, aml$status)

m1_x <- survreg(s ~ 1, dist='extreme')
m1_l <- survreg(s ~ 1, dist='logistic')
m1_w <- survreg(s ~ 1, dist='weibull')
m1_e <- survreg(s ~ 1, dist='exponential')
m1_r <- survreg(s ~ 1, dist='rayleigh')
m1_lg <- survreg(s ~ 1, dist='loggaussian')
m1_ln <- survreg(s ~ 1, dist='lognormal')
m1_ll <- survreg(s ~ 1, dist='loglogistic')
m1_t <- survreg(s ~ 1, dist='t')
m1_c <- survreg(s ~ 1, dist=mycauchy)


mods <- list(
  extr= m1_x,
  logi = m1_l,
  weib = m1_w,
  expo = m1_e,
  rayl = m1_r,
  logg = m1_lg,
  logn = m1_ln,
  logl = m1_ll,
  tdis = m1_t,
  cauc = m1_c
)

library(modelsummary)

modelsummary::modelsummary(mods)



library(flexsurv)

m2_gg  <- flexsurvreg(s ~ 1, dist ='gengamma')
m2_ggo <- flexsurvreg(s ~ 1, dist ='gengamma.orig')
m2_gf  <- flexsurvreg(s ~ 1, dist ='genf')
m2_gfo <- flexsurvreg(s ~ 1, dist ='genf.orig')
m2_w   <- flexsurvreg(s ~ 1, dist ='weibull')
m2_g   <- flexsurvreg(s ~ 1, dist ='gamma')
m2_e   <- flexsurvreg(s ~ 1, dist='exp')
m2_ll  <- flexsurvreg(s ~ 1, dits='llogis')
m2_ln  <- flexsurvreg(s ~ 1, dist='lnorm')
m2_go  <- flexsurvreg(s ~ 1, dist='gompertz')


m2_sp1 <- flexsurvspline(s ~ 1, k=1, scale='hazard')
m2_sp2 <- flexsurvspline(s ~ 1, k=1, scale='odds')
m2_sp3 <- flexsurvspline(s ~ 1, k=1, scale='normal')
m2_sp4 <- flexsurvspline(s ~ 1, k=3, scale='odds')
m2_sp5 <- flexsurvspline(s ~ 1, k=10, scale='odds')
## https://cran.r-project.org/web/packages/flexsurv/vignettes/flexsurv.pdf



fs_betapr <- list(
  name = 'betapr',
  pars= c('shape1','shape2','scale'),
  location='scale', # scale or location param
  transforms= c(log,log,log),
  inv.transforms= c(exp,exp,exp),
  inits= function(t){c(mean(t),1/mean(t) + 1,1)}
)



flexsurvreg(s~1, dist=fs_betapr)

# flexsurvmix for competing risk

# flexsurvtrunc for double events (event + death, etc)

mlist2 <- list(
  gengam = m2_gg, 
  gengam.o = m2_ggo,
  genf = m2_gf,
  gen.o = m2_gfo,
  weib = m2_w,
  gamm = m2_g,
  expo = m2_e,
  logl = m2_ll,
  logn = m2_ln,
  gomp = m2_go,
  spl_h = m2_sp1,
  spl_o = m2_sp2,
  spl_n = m2_sp3,
  spl_o3 = m2_sp4,
  spl_o10 = m2_sp5
)