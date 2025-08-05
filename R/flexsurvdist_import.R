library(dplyr)
library(survival)
library(flexsurv)

# Rcpp::sourceCpp('dist.cpp')


#### Custom distribution shell
if(F){
  custom_model <- list(
    name = '',               # this should be the shortened name that appears after d/p in the distribution functions i.e. "norm" in dnorm
    pars = c('','',''),      # a vector of each parameter name, exactly as they are stated in the d/p functions i.e c('mean','sd') for norm
    location= c(),           # if there is a location or scale parameter, you can place it here. It's used in the scaling on variable in more complicated models
    transforms = c(),        # a vector of functions, one for each parameter, that transform the respective parameter into the real line. If param1 must be >0, log(param1) would be real
    inv.transforms=c(),      # a vector of functions which transform the real values back to the original scale. exp(log(param1) returns to the original range of param1
    inits = function(t){},   # a declared function which takes the time variable, and outputs initial estimates for the parameters. Ideally, estimates are based off of t somehow, but can be simple ints
    p = pdist,               # by default, flexsurv looks for the d/p functions using the <name> object above, but the shotgun function will use these declared p/d's instead
    d = ddist,               #
    fullname = ''            # this longer name is used in shotgun to better label the distribution in outputs and console tracking
  )
}











### Beta Prime
# extraDistr has as betapr
flexsurv_betapr <- list(
  name = 'betapr',
  pars= c('shape1','shape2','scale'),
  location='scale',
  transforms= c(log,log,log),
  inv.transforms= c(exp,exp,exp),
  inits= function(t){c(3, 2, 1)},
  p = extraDistr::pbetapr,
  d = extraDistr::dbetapr,
  fullname='beta_prime'
)

### Birnbaum - Saunders
# extraDistr has this as fatigue
flexsurv_fatigue <- list(
  name='fatigue',
  pars = c('alpha','beta'), # 'mu' dropped
  location= 'alpha',  # previously mu
  transforms = c(log,log),
  inv.transforms = c(exp,exp),
  inits = function(t){c(mean(t)/2,1)},
  p = extraDistr::pfatigue,
  d = extraDistr::dfatigue,
  fullname='birnbaum_saunders'
)



# Singh Madalla / burr t12
flexsurv_sinmad  <- list(
  name='sinmad',
  pars = c('scale','shape1.a','shape3.q'),
  location = ('scale'),
  transforms = c(log,log,log),
  inv.transforms = c(exp,exp,exp),
  inits=function(t){c(1,1,1)},
  d=VGAM::dsinmad,
  p=VGAM::psinmad,
  fullname='singh_maddala'
)




### CAUCHY (will suck)
# basic stats package has,
flexsurv_cauchy <- list(
  name='cauchy',
  pars= c('location','scale'),
  location= c('location'),
  transforms = c(identity, log),
  inv.transforms = c(identity, exp),
  inits = function(t){c(median(t), mad(t)+1e-8)},
  d= dcauchy,
  p= pcauchy,
  fullname='cauchy'
)

### CHI SQUARED
# basic inclusion
flexsurv_chisq <- list(
  name='chisq',
  pars=c('df'),
  location=c('df'),
  transforms=c(log),
  inv.transforms=c(exp),
  inits = function(t){c(mean(t))},
  d = dchisq,
  p = pchisq,
  fullname='chi_squared'
)

# non central
flexsurv_ncchisq <- list(
  name='chisq',
  pars=c('df','ncp'),
  location=c('ncp'),
  transforms=c(log, log),
  inv.transforms=c(exp,exp),
  inits = function(t){c(mean(t), 1)},
  d = dchisq,
  p = pchisq,
  fullname='non_central_chi_squared'
)



### Dagum
# from vgam
flexsurv_dagum <- list(
  name='dagum',
  pars= c('scale','shape1.a','shape2.p'),
  location=c('scale'),
  transforms = c(log,log,log),
  inv.transforms = c(exp,exp,exp),
  inits = function(t){c(1,1,1)},
  p = VGAM::pdagum,
  d = VGAM::ddagum,
  fullname='dagum'
)





### Exponential logarithmic
# VGAM
logit <- function(p){log(p/(1-p))}
inv_logit <- function(x){return(1/(1+exp(-x)))}
flexsurv_explog <- list(
  name = 'explog',
  pars = c('scale','shape'),
  location = 'scale',
  transforms = c(log, logit),
  inv.transforms = c(exp, inv_logit),
  inits=function(t){c(median(t),.5)},
  p=VGAM::pexplog,
  d=VGAM::dexplog,
  fullname='exponential_logarithmic'
)


### Extreme value (generalized)
flexsurv_extval <- list(
  name='gev',
  pars = c('mu','sigma','xi'),
  location='mu',
  transforms= c(identity, log, log),  # trying to make xi log now, as the bounding is bad if <=0
  inv.transforms = c(identity, exp, exp),
  inits=function(t){c(1,1,1)},
  d=extraDistr::dgev,
  p=extraDistr::pgev,
  fullname='generalized_extreme_value'
)


### F distrib (idk)
flexsurv_f <- list(
  name='f',
  pars = c('df1','df2'),
  location = ('df1'),
  transforms = c(log,log),
  inv.transforms = c(exp,exp),
  inits=function(t){c(median(t),2)},
  d=df,
  p=pf,
  fullname='f'
)

#noncentral
flexsurv_ncf <- list(
  name='f',
  pars = c('df1','df2', 'ncp'),
  location = ('ncp'),
  transforms = c(log,log,log),
  inv.transforms = c(exp,exp,exp),
  inits=function(t){c(median(t),2,0.1)},
  d=df,
  p=pf,
  fullname='noncentral_f'
)



### folded normal?
flexsurv_fnorm <- list(
  name= 'foldnorm',
  pars = c('mean','sd','a1','a2'),
  location = c('mean'),
  transforms = c(identity, log, identity, identity),
  inv.transforms = c(identity, exp, identity, identity),
  inits = function(t){c(mean(t),sd(t),1,1)},
  d= VGAM::dfoldnorm,
  p= VGAM::pfoldnorm,
  fullname='folded_normal'
)


### Frechet
flexsurv_frech <- list(
  name='frechet',
  pars = c('lambda','mu','sigma'), # shape, location, scale
  location='mu',
  transforms= c(log,identity,log),
  inv.transforms = c(exp,identity,exp),
  inits = function(t){c(1.5,min(t)-1,1)},
  d= extraDistr::dfrechet,
  p= extraDistr::pfrechet,
  fullname='frechet'
)



### Gamma-Gompertz
flexsurv_gamgomp <- list(
  name='gamgomp',
  pars=c('b','s','beta'),
  location='b',
  transforms=c(log,log,log),
  inv.transforms=c(exp,exp,exp),
  inits=function(t){c(log(2)/median(t),1,IQR(t)/median(t))},
  d=dgamgomp,
  p=pgamgomp,
  fullname='gamma_gompertz'
)





### Gumbel (basically log weibull)
# vgam extraDistr
flexsurv_gumbel <- list(
  name='gumbel',
  pars=c('mu','sigma'),
  location='mu',
  transforms=c(identity, log),
  inv.transforms=c(identity, exp),
  inits=function(t){c(mean(t),(mean(t)/.577))},
  d=extraDistr::dgumbel,
  p=extraDistr::pgumbel,
  fullname='gumbel'
)



### Inv Chi-sq
# extradistr
flexsurv_invchisq <- list(
  name='invchisq',
  pars=c('nu'),
  location='nu',
  transforms=c(log),
  inv.transforms=c(exp),
  inits=function(t){c(3)},
  d=extraDistr::dinvchisq,
  p=extraDistr::pinvchisq,
  fullname='inverse_chi_squared'
)


# scaled version too
flexsurv_sinvchisq <- list(
  name='invchisq',
  pars=c('nu','tau'),
  location='tau',
  transforms=c(log, log),
  inv.transforms=c(exp, exp),
  inits=function(t){c(3, 1)},
  d=extraDistr::dinvchisq,
  p=extraDistr::pinvchisq,
  fullname='scaled_inverse_chi_squared'
)


### Inv Gamma
# extradistr
flexsurv_invgam <- list(
  name='invgamma',
  pars=c('alpha','beta'),
  location='beta',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(1,1)},
  d=extraDistr::dinvgamma,
  p=extraDistr::pinvgamma,
  fullname='inverse_gamma'
)


### Inv Gaussian
# VGAM
flexsurv_invgaus <- list(
  name='inv.gaussian',
  pars=c('mu','lambda'),
  location='mu',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(median(t)/10, 2)},
  d=VGAM::dinv.gaussian,
  p=VGAM::pinv.gaussian,
  fullname='inverse_gaussian'
)


### Laplace
# extraDistr
flexsurv_laplace <- list(
  name='laplace',
  pars=c('mu','sigma'),
  location='mu',
  transforms=c(identity, log),
  inv.transforms=c(identity, exp),
  inits=function(t){c(mean(t), sqrt(var(t)/2))},
  d=extraDistr::dlaplace,
  p=extraDistr::plaplace,
  fullname='laplace'
)


### Levy?
# VGAM
flexsurv_levy <- list(
  name='levy',
  pars=c('location','scale'),
  location='location',
  transforms=c(identity, log),
  inv.transforms=c(identity, exp),
  inits=function(t){c(0 - 2e-4, median(t))},
  d=VGAM::dlevy,
  p=VGAM::plevy,
  fullname='levy'
)


### Log cauchy
# could be interesting

flexsurv_logcauchy <- list(
  name='logcauchy',
  pars=c('mu','sigma'),
  location='mu',
  transforms=c(identity, log),
  inv.transforms=c(identity, exp),
  inits=function(t){c(1,1)},
  d=dlogcauchy,
  p=plogcauchy,
  fullname='log_cauchy'
)


### Lomax
flexsurv_lomax <- list(
  name='lomax',
  pars=c('lambda','kappa'),
  location='lambda',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(median(t),1)},
  d=extraDistr::dlomax,
  p=extraDistr::plomax,
  fullname='lomax'
)


### nakagami
flexsurv_nakagami <- list(
  name='naka',
  pars=c('scale','shape'),
  location='scale',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(mean(t^2),(mean(t^2)^2)/var(t^2))},
  d=VGAM::dnaka,
  p=VGAM::pnaka,
  fullname='nakagami'
)


### PARETO, THERES A FEW HERE
flexsurv_pareto1 <- list(
  name='paretoI',
  pars=c('scale','shape'),
  location='scale',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(min(t)*.95,1)},
  d=VGAM::dparetoI,
  p=VGAM::pparetoI,
  fullname='pareto_type_1'
)

flexsurv_pareto2 <- list(
  name='paretoII',
  pars=c('location','scale','shape'),
  location='location',
  transforms=c(identity,log,log),
  inv.transforms=c(identity,exp,exp),
  inits=function(t){c(0,min(t)*.95,1)},
  d=VGAM::dparetoII,
  p=VGAM::pparetoII,
  fullname='pareto_type_2'
)

flexsurv_pareto3 <- list(
  name='paretoIII',
  pars=c('location','scale','inequality'),
  location='location',
  transforms=c(identity,log,log),
  inv.transforms=c(identity,exp,exp),
  inits=function(t){c(0,min(t)*.95,1)},
  d=VGAM::dparetoIII,
  p=VGAM::pparetoIII,
  fullname='pareto_type_3'
)

flexsurv_pareto4 <- list(
  name='paretoIV',
  pars=c('location','scale','inequality','shape'),
  location='location',
  transforms=c(identity,log,log,log),
  inv.transforms=c(identity,exp,exp,exp),
  inits=function(t){c(0,min(t)*95,1,2)},
  d=VGAM::dparetoIV,
  p=VGAM::pparetoIV,
  fullname='pareto_type_4'
)

flexsurv_fpareto <- list(
  name='fpareto',
  pars=c('min','shape1','shape2','shape3','rate'),
  location='min',
  transforms=c(identity,log,log,log,log),
  inv.transforms=c(identity,exp,exp,exp,exp),
  inits=function(t){c(0,1,1,1,1)},
  d=actuar::dfpareto,
  p=actuar::pfpareto,
  fullname='feller_pareto'
)





flexsurv_tpareto <- list(
  name='truncpareto',
  pars=c('lower','upper','shape'),
  location='lower',
  transforms=c(identity,identity, log),
  inv.transforms=c(identity, identity, exp),
  inits=function(t){c(min(t)-.0001,max(t)+1, 1)},
  d=VGAM::dtruncpareto,
  p=VGAM::ptruncpareto,
  fullname='truncated_pareto'
)


### Rayleigh
flexsurv_rayleigh <- list(
  name='rayleigh',
  pars=c('sigma'),
  location='sigma',
  transforms=c(log),
  inv.transforms=c(exp),
  inits=function(t){c(median(t)/sqrt(2*log(2)))},
  d=extraDistr::drayleigh,
  p=extraDistr::prayleigh,
  fullname='rayleigh'
)

### Rice
flexsurv_rice <- list(
  name='rice',
  pars=c('sigma','vee'),
  location='vee',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(median(t), sd(t))},
  d=VGAM::drice,
  p=VGAM::price,
  fullname='rice'
)


### Shifted Gompertz
flexsurv_sgomp <- list(
  name='sgomp',
  pars=c('b','eta'),
  location='b',
  transforms=c(log,log),
  inv.transforms=c(exp,exp),
  inits=function(t){c(1/median(t),median(t))},
  d=extraDistr::dsgomp,
  p=extraDistr::psgomp,
  fullname='shifted_gompertz'
)


### Type-2 Gumbel
flexsurv_gumbelII <- list(
  name='gumbelII',
  pars=c('scale','shape'),
  location='scale',
  transforms=c(identity, identity),
  inv.transforms=c(identity, identity),
  inits=function(t){c(1,1)},
  d=VGAM::dgumbelII,
  p=VGAM::pgumbelII,
  fullname='gumbel_type_2'
)


# hyper tabastic
flexsurv_hypertab <- list(
  name='hypertab',
  pars= c('a','b'),
  location= 'a', # could also be B, worth trying
  transforms= c(log,log),
  inv.transforms= c(exp,exp),
  inits= function(t){c(1,1/mean(t))},
  d=dhypertab,
  p=phypertab,
  fullname='hypertabastic'
)

# lindley distribution
flexsurv_lindley <- list(
  name= 'lind',
  pars=c('theta'),
  location='theta',
  transforms=c(log),
  inv.transforms=c(exp),
  inits= function(t){c(1/mean(t))},
  d=VGAM::dlind,
  p=VGAM::plind,
  fullname='lindley'
)



# inverse lindley
flexsurv_invlind <- list(
  name= 'invlind',
  pars=c('theta'),
  location='theta',
  transforms=c(log),
  inv.transforms=c(exp),
  inits= function(t){c(1)},
  d=dinvlind,
  p=pinvlind,
  fullname='inverse_lindley'
)



surv_shotgun <- function(formula, data=NA, skip=c('default'), dump_models =F, progress=T, warn=F){
  require(flexsurv)
  require(dplyr)
  require(actuar)
  require(VGAM)
  require(extraDistr)
  require(tictoc)


  if(length(skip)==1 & skip[1]=='default'){
    skip <- c('genf.orig', 'gengamma.orig', 'erlang', 'truncpareto', 'chisq','non_central_chi_squared', 'exponential','lognormal')
  }


  message("Loading shotgun...")
  tic("Kapow")

  vars <- all.vars(formula)
  in_global <- sapply(vars, exists, envir=.GlobalEnv)
  data_req <- F
  if(!all(in_global)){
    if(is.data.frame(data)){
      if(all(vars %in% names(data))){
        data_req <- T
      }else{
        message("it broke")
        stop("you")
      }
    }else{
      stop("not in the data")
    }
  }

  my_dist_list <- list(
    flexsurv_betapr, flexsurv_fatigue, flexsurv_sinmad, flexsurv_cauchy, flexsurv_chisq,
    flexsurv_ncchisq, flexsurv_dagum,  flexsurv_explog, flexsurv_extval, flexsurv_f,
    flexsurv_ncf, flexsurv_fnorm,flexsurv_frech, flexsurv_gamgomp, flexsurv_gumbel, flexsurv_invchisq,
    flexsurv_sinvchisq, flexsurv_invgam, flexsurv_invgaus, flexsurv_laplace, flexsurv_levy, flexsurv_logcauchy,
    flexsurv_lomax, flexsurv_nakagami, flexsurv_pareto1, flexsurv_pareto2, flexsurv_pareto3, flexsurv_pareto4,
    flexsurv_fpareto, flexsurv_tpareto, flexsurv_rayleigh, flexsurv_rice, flexsurv_sgomp, flexsurv_gumbelII,
    flexsurv_hypertab,flexsurv_lindley, flexsurv_invlind
    # flexsurv_pareto, #flexsurv_burr,#flexsurv_erlang,
  )

  dist_list <- c(flexsurv.dists, my_dist_list)

  dist_summary <- data.frame(
    dist_name = '1',
    dist_ran = T,
    aic = 1,
    bic = 1,
    loglik = 1
  )[-1,]

  iter <- 1
  for(i in dist_list){
    if(i$name %in% skip | coalesce(i$fullname, i$name) %in% skip | names(dist_list)[iter] %in% skip){
      iter <- iter+1
      next
    }

    current_dist <- coalesce(i$fullname, i$name)
    if(progress){
      message('------------------------------------------------------------------------')
      message(current_dist)
      tic(current_dist)
    }

    dist_success<- F
    current_aic <- NA
    current_bic <- NA
    current_ll <- NA

    tryCatch({
        withCallingHandlers({
          if(data_req){
            flexsurvreg(formula, dist=i, data=data, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
          }else{
            flexsurvreg(formula, dist=i, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
          }
          current_aic <- AIC(current_model)
          current_bic <- BIC(current_model)
          current_ll <- current_model$loglik



          dist_success<-T
          },
          warning=function(w){
            if(warn){message(paste('Warning in', current_dist,'model :',w))}
            invokeRestart('muffleWarning')
          })},
        error=function(e){
          message(paste("Error in", current_dist,"model :",e))
        }
    )

    dist_summary %>%
      add_row(
        dist_name = current_dist, dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
      ) -> dist_summary

    if(dump_models & dist_success){
      # flexsurv shotgun
      assign(paste('fssg_',current_dist,sep='',collapse=''), current_model, envir = .GlobalEnv)
    }

    toc(quiet =T)$callback_msg %>% message()
    iter = iter+1
  }


  gc()
  message('------------------------------------------------------------------------')
  message('Final Summary!')


  dist_summary %>% filter(!is.na(aic)) %>% filter(aic== min(aic)) %>% select(dist_name) %>% pull() -> best_aic
  dist_summary %>% filter(!is.na(aic)) %>% filter(bic== min(bic)) %>% select(dist_name) %>% pull() -> best_bic
  dist_summary %>% filter(!is.na(aic)) %>% filter(loglik== max(loglik)) %>% select(dist_name) %>% pull() -> best_ll

  message(paste('Model with best AIC:',best_aic))
  message(paste('Model with best BIC:',best_bic))
  message(paste('Model with best Log Likelihood:',best_ll))
  toc(quiet=T)$callback_msg %>% message()
  return(dist_summary)
}






# raw testing from legacy testing
if(F){


  # quick tests
  formula_test <- Surv(time, status)~1
  surv_shotgun(formula_test, data=aml, dump_models=T, warn = T) -> shotgun_summary
  shotgun_summary %>% arrange(aic, bic, desc(loglik))

  # breaks my custom distribs :c
  formula_test2 <- Surv(time,status)~ x
  surv_shotgun(formula_test2, data=aml, dump_models=T, warn = T) -> shotgun_summary2
  shotgun_summary2 %>% arrange(aic, bic, desc(loglik))
}
