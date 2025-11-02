#### Custom distribution shell
# if(F){
#   custom_model <- list(
#     name = '',               # this should be the shortened name that appears after d/p in the distribution functions i.e. "norm" in dnorm
#     pars = c('','',''),      # a vector of each parameter name, exactly as they are stated in the d/p functions i.e c('mean','sd') for norm
#     location= c(),           # IMPORTANT: specifies which parameters should vary based on the covariates in the model. Can usually omit shape /
#                              #  true-locations (like minimum or boundary parameters). This also is a key part of specifying AFT/PH interpretation of your model
#     transforms = c(),        # a vector of functions, one for each parameter, that transform the respective parameter into the real line. If param1 must be >0, log(param1) would be real
#     inv.transforms=c(),      # a vector of functions which transform the real values back to the original scale. exp(log(param1) returns to the original range of param1
#     inits = function(t){},   # a declared function which takes the time variable, and outputs initial estimates for the parameters. Ideally, estimates are based off of t somehow, but can be simple ints
#     p = pdist,               # by default, flexsurv looks for the d/p functions using the <name> object above, but the fssg function will use these declared p/d's instead
#     d = ddist,               #
#     h = hdist,               # May as well include the hazard functions as well!
#     H = Hdist,               #
#     fullname = ''            # this longer name is used in fssg to better label the distribution in outputs and console tracking
#   )
# }


### Important note: flexsurv by default only varies one model parameter (what is specified in the distributions as location)
# We can make more than one parameter vary though, using the anc parameter in flexsurvreg
# Example: anc = list(shape1 = ~ var1 + var2, shape2 = ~ var3)
# It may be worth adding in some anc options in the dist-list...



### Function-Factory
# Couple of functions that extrapolate the S, h, and H functions based on pdf/cdf
#
# S<dist> = 1-p<dist>
# h<dist> = d<dist> / S<dist>
# H<dist> = -log(S<dist>)
survivify <- function(p_func){
  force(p_func)  # ensure it's evaluated in the parent environment
  function(...) {
    1 - p_func(...)
  }
}

hazardify <- function(d_func, p_func){
  force(d_func)
  force(p_func)
  function(x, ...){
    fx <- d_func(x, ...)
    Fx <- p_func(q=x, ...)
    fx/ (1-Fx)
  }
}

cumhazardify <- function(p_func){
  force(p_func)
  function(x, ...){
    sx <- 1 - p_func(q=x,...)
    -log(sx)
  }
}



#' Function to check if times can be calculated using the distribution with default inits.
#'
#' @param times Surv object or numeric vector.
#' @param distribution A distribution object from fssg_dist_lists.
#' @returns Boolean indicator for success. If true, then all values can be calculated, and life is good.
#' @export
check_inits <- function(times, distribution){
  # accepts a vector or Surv object
  if(methods::is(times,'Surv')){
    times[,1] %>%
      as.numeric() %>%
      unique() %>%
      sort() -> time_vector
  }else{
    times %>%
      as.numeric() %>%
      unique() %>%
      sort() -> time_vector
  }

  if('d' %in% names(distribution)){
    dfunc <- distribution$d
  }else{
    dfunc <- get(paste('d', distribution$name, sep=''))
  }

  inits <- distribution$inits(time_vector)

  dfunc %>%
    formals() %>%
    names() %>%
    setdiff(c('log')) -> arguments


  dataframe <- data.frame(t = time_vector, p = NA, s = F)

  for(i in 1:length(time_vector)){
    as.list(c(time_vector[i], inits)) %>%
      stats::setNames(c(arguments)) %>%
      do.call(what=dfunc, args=.) -> dataframe$p[i]

    dataframe$s[i] <- is.finite(dataframe$p[i])
  }

  if(!any(dataframe$s)){
    message("Found some weird entries")

    dataframe %>%
      dplyr::filter(s=F) %>%
      print()
  }else{
    message("Works at every time point!")
    plot(
      x = dataframe$t[dataframe$s],
      y = dataframe$p[dataframe$s],
      type='o', main = 'Distribution with Default Inits',xlab = 'Times', ylab='Probability Density',sub=substitute(density_function))
  }

  return(any(dataframe$s))
}






#' Compiles list of available distributions
#'
#' @returns a list of all possible distributions
#'
#' @export
fssg_dist_list <- function(){
  # create a list item for each currently available distribution

  ### Beta Prime
  flexsurv_betapr <- list(
    name = 'betapr',
    pars= c('shape1','shape2','scale'),
    location='scale',
    transforms= c(log,log,log),
    inv.transforms= c(exp,exp,exp),
    inits= function(t){c(3, 2, 1)},  # can be improved
    d = extraDistr::dbetapr,
    p = extraDistr::pbetapr,
    h = hazardify(extraDistr::dbetapr,extraDistr::pbetapr),
    H = cumhazardify(extraDistr::pbetapr),
    fullname='beta_prime'
  )

  ### Birnbaum - Saunders, specified as fatigue in extraDistr
  flexsurv_fatigue_shape <- list(
    name='fatigue',
    pars = c('alpha','beta'), # shape, scale,  (mu is omitted but can be location. set to zero by default)
    location= 'alpha',  # shape, PH effect?
    transforms = c(log,log),
    inv.transforms = c(exp,exp),
    inits = function(t){c(mean(t)/2,1)},
    d = extraDistr::dfatigue,
    p = extraDistr::pfatigue,
    h = hazardify(extraDistr::dfatigue,extraDistr::pfatigue),
    H = cumhazardify(extraDistr::pfatigue),
    fullname='birnbaum_saunders_shape'
  )

  flexsurv_fatigue <- list(
    name='fatigue',
    pars = c('alpha','beta'), # shape, scale,  (mu is omitted but can be location. set to zero by default)
    location= 'beta',  # we'll use scale: AFT effect?
    transforms = c(log,log),
    inv.transforms = c(exp,exp),
    inits = function(t){c(mean(t)/2,1)},
    d = extraDistr::dfatigue,
    p = extraDistr::pfatigue,
    h = hazardify(extraDistr::dfatigue,extraDistr::pfatigue),
    H = cumhazardify(extraDistr::pfatigue),
    fullname='birnbaum_saunders'
  )

  flexsurv_fatigue_shape_loc <- list(
    name='fatigue',
    pars = c('alpha','beta', 'mu'), # shape, scale, location
    location= 'alpha',  # we'll use shape: PH Effect?
    transforms = c(log,log, identity),
    inv.transforms = c(exp,exp, identity),
    inits = function(t){c(mean(t)/2,1,0)},
    d = extraDistr::dfatigue,
    p = extraDistr::pfatigue,
    h = hazardify(extraDistr::dfatigue,extraDistr::pfatigue),
    H = cumhazardify(extraDistr::pfatigue),
    fullname='birnbaum_saunders_shape_location'
  )

  flexsurv_fatigue_loc <- list(
    name='fatigue',
    pars = c('alpha','beta', 'mu'), # shape, scale, location
    location= 'beta',  # we'll use scale: AFT effect?
    transforms = c(log,log, identity),
    inv.transforms = c(exp,exp, identity),
    inits = function(t){c(mean(t)/2,1,0)},
    d = extraDistr::dfatigue,
    p = extraDistr::pfatigue,
    h = hazardify(extraDistr::dfatigue,extraDistr::pfatigue),
    H = cumhazardify(extraDistr::pfatigue),
    fullname='birnbaum_saunders_location'
  )

  ### Singh Madalla aka Burr-T12
  flexsurv_sinmad  <- list(
    name='sinmad',
    pars = c('scale','shape1.a','shape3.q'), # scale, shape, shape
    location = ('scale'),
    transforms = c(log,log,log),
    inv.transforms = c(exp,exp,exp),
    inits=function(t){c(1,1,1)},
    d=VGAM::dsinmad,
    p=VGAM::psinmad,
    h = hazardify(VGAM::dsinmad, VGAM::psinmad),
    H = cumhazardify(VGAM::psinmad),
    fullname='singh_maddala'
  )

  ### CAUCHY
  flexsurv_cauchy <- list(
    name='cauchy',
    pars= c('location','scale'), # location, shape
    location= c('scale'),
    transforms = c(identity, log),
    inv.transforms = c(identity, exp),
    inits = function(t){c(stats::median(t), stats::mad(t)+1e-8)},
    d= stats::dcauchy,
    p= stats::pcauchy,
    h = hazardify(stats::dcauchy,  stats::pcauchy),
    H = cumhazardify(stats::pcauchy),
    fullname='cauchy'
  )

  flexsurv_cauchy_loc <- list(  # alternative location based approach?
    name='cauchy',
    pars= c('location','scale'),
    location= c('location'),
    transforms = c(identity, log),
    inv.transforms = c(identity, exp),
    inits = function(t){c(stats::median(t), stats::mad(t)+1e-8)},
    d= stats::dcauchy,
    p= stats::pcauchy,
    h = hazardify(stats::dcauchy,  stats::pcauchy),
    H = cumhazardify(stats::pcauchy),
    fullname='cauchy_location'
  )

  ### CHI SQUARED
  flexsurv_chisq <- list(
    name='chisq',
    pars=c('df'),
    location=c('df'),
    transforms=c(log),
    inv.transforms=c(exp),
    inits = function(t){c(mean(t))},
    d = stats::dchisq,
    p = stats::pchisq,
    h = hazardify(stats::dchisq, stats::pchisq),
    H = cumhazardify(stats::pchisq),
    fullname='chi_squared'
  )

  ### non central chi squared
  flexsurv_ncchisq <- list(
    name='chisq',
    pars=c('df','ncp'),  # shape scale
    location=c('ncp'),
    transforms=c(log, log),
    inv.transforms=c(exp,exp),
    inits = function(t){c(mean(t), 1)},
    d = stats::dchisq,
    p = stats::pchisq,
    h = hazardify(stats::dchisq, stats::pchisq),
    H = cumhazardify(stats::pchisq),
    fullname='non_central_chi_squared'
  )

  ### Dagum
  flexsurv_dagum <- list(
    name='dagum',
    pars= c('scale','shape1.a','shape2.p'), # scale, shape, shape
    location=c('scale'),
    transforms = c(log,log,log),
    inv.transforms = c(exp,exp,exp),
    inits = function(t){c(1,1,1)},
    p = VGAM::pdagum,
    d = VGAM::ddagum,
    h = hazardify(VGAM::ddagum, VGAM::pdagum),
    H = cumhazardify(VGAM::pdagum),
    fullname='dagum'
  )

  ### Exponential logarithmic
  logit <- function(p){log(p/(1-p))}
  inv_logit <- function(x){return(1/(1+exp(-x)))}
  flexsurv_explog <- list(
    name = 'explog',
    pars = c('scale','shape'), # scale, shape
    location = 'scale',
    transforms = c(log, logit),
    inv.transforms = c(exp, inv_logit),
    inits=function(t){c(stats::median(t),.5)},
    p=VGAM::pexplog,
    d=VGAM::dexplog,
    h = hazardify(VGAM::dexplog, VGAM::pexplog),
    H = cumhazardify(VGAM::pexplog),
    fullname='exponential_logarithmic'
  )

  ### Extreme value (generalized)
  flexsurv_extval <- list(
    name='gev',
    pars = c('mu','sigma','xi'), # loc, scale, shape
    location='sigma',
    transforms= c(identity, log, log),  # trying to make xi use log now, as the bounding is bad if <=0
    inv.transforms = c(identity, exp, exp),
    inits=function(t){c(1,1,1)},
    d=extraDistr::dgev,
    p=extraDistr::pgev,
    h = hazardify(extraDistr::dgev, extraDistr::pgev),
    H = cumhazardify(extraDistr::pgev),
    fullname='generalized_extreme_value'
  )

  ### F distribution
  flexsurv_f <- list(
    name='f',
    pars = c('df1','df2'),
    location = ('df1'),
    transforms = c(log,log),
    inv.transforms = c(exp,exp),
    inits=function(t){c(stats::median(t),2)},
    d=stats::df,
    p=stats::pf,
    h = hazardify(stats::df, stats::pf),
    H = cumhazardify(stats::pf),
    fullname='f'
  )

  ### noncentral F
  flexsurv_ncf <- list(
    name='f',
    pars = c('df1','df2', 'ncp'),
    location = ('ncp'),
    transforms = c(log,log,log),
    inv.transforms = c(exp,exp,exp),
    inits=function(t){c(stats::median(t),2,0)},
    d=stats::df,
    p=stats::pf,
    h = hazardify(stats::df, stats::pf),
    H = cumhazardify(stats::pf),
    fullname='noncentral_f'
  )

  ### folded normal parameterized on mean (location)
  flexsurv_fnorm_loc <- list(
    name= 'foldnorm',
    pars = c('mean','sd','a1','a2'),
    location = c('mean'),
    transforms = c(identity, log, identity, identity),
    inv.transforms = c(identity, exp, identity, identity),
    inits = function(t){c(mean(t),stats::sd(t),1,1)},
    d= VGAM::dfoldnorm,
    p= VGAM::pfoldnorm,
    h = hazardify(VGAM::dfoldnorm, VGAM::pfoldnorm),
    H = cumhazardify(VGAM::pfoldnorm),
    fullname='folded_normal_location'
  )

  ### Folded normal, parameterized on scale
  flexsurv_fnorm <- list(
    name= 'foldnorm',
    pars = c('mean','sd','a1','a2'),
    location = c('sd'),
    transforms = c(identity, log, identity, identity),
    inv.transforms = c(identity, exp, identity, identity),
    inits = function(t){c(mean(t),stats::sd(t),1,1)},
    d= VGAM::dfoldnorm,
    p= VGAM::pfoldnorm,
    h = hazardify(VGAM::dfoldnorm, VGAM::pfoldnorm),
    H = cumhazardify(VGAM::pfoldnorm),
    fullname='folded_normal'
  )

  ### Frechet
  flexsurv_frech <- list(
    name='frechet',
    pars = c('lambda','mu','sigma'), # shape, location, scale
    location='sigma',
    transforms= c(log,identity,log),
    inv.transforms = c(exp,identity,exp),
    inits = function(t){c(1.5,min(t)-1,1)},
    d= extraDistr::dfrechet,
    p= extraDistr::pfrechet,
    h = hazardify(extraDistr::dfrechet, extraDistr::pfrechet),
    H = cumhazardify(extraDistr::pfrechet),
    fullname='frechet'
  )

  ### Gamma-Gompertz
  flexsurv_gamgomp <- list(
    name='gamgomp',
    pars=c('b','sigma','beta'), # scale, shape, shape
    location='b',
    transforms=c(log,log,log),
    inv.transforms=c(exp,exp,exp),
    inits=function(t){c(1/max(t), 1, 1)},
    d=dgamgomp,
    p=pgamgomp,
    h = hazardify(dgamgomp, pgamgomp),
    H = cumhazardify(pgamgomp),
    fullname='gamma_gompertz'
  )

  ### Gumbel (basically log weibull)
  flexsurv_gumbel <- list(
    name='gumbel',
    pars=c('mu','sigma'), # location, scale
    location='sigma',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(mean(t),(mean(t)/.577))},
    d=extraDistr::dgumbel,
    p=extraDistr::pgumbel,
    h = hazardify(extraDistr::dgumbel, extraDistr::pgumbel),
    H = cumhazardify(extraDistr::pgumbel),
    fullname='gumbel'
  )

  ### Inverse Chi-squared
  flexsurv_invchisq <- list(
    name='invchisq',
    pars=c('nu'),
    location='nu',
    transforms=c(log),
    inv.transforms=c(exp),
    inits=function(t){c(3)},
    d=extraDistr::dinvchisq,
    p=extraDistr::pinvchisq,
    h = hazardify(extraDistr::dinvchisq, extraDistr::pinvchisq),
    H = cumhazardify(extraDistr::pinvchisq),
    fullname='inverse_chi_squared'
  )

  # scaled version inverse Chi Squared
  flexsurv_sinvchisq <- list(
    name='invchisq',
    pars=c('nu','tau'),
    location='tau',
    transforms=c(log, log),
    inv.transforms=c(exp, exp),
    inits=function(t){c(3, 1)},
    d=extraDistr::dinvchisq,
    p=extraDistr::pinvchisq,
    h = hazardify(extraDistr::dinvchisq, extraDistr::pinvchisq),
    H = cumhazardify(extraDistr::pinvchisq),
    fullname='scaled_inverse_chi_squared'
  )


  ### Inverse Gamma
  flexsurv_invgam <- list(
    name='invgamma',
    pars=c('alpha','beta'), # shape, scale
    location='beta',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(1,1)},
    d=extraDistr::dinvgamma,
    p=extraDistr::pinvgamma,
    h = hazardify(extraDistr::dinvgamma, extraDistr::pinvgamma),
    H = cumhazardify(extraDistr::pinvgamma),
    fullname='inverse_gamma'
  )

  ### Inverse Gaussian
  flexsurv_invgaus_loc <- list(
    name='inv.gaussian',
    pars=c('mu','lambda'), # mean, scale?
    location='mu',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t)/10, 2)},
    d=VGAM::dinv.gaussian,
    p=VGAM::pinv.gaussian,
    h = hazardify(VGAM::dinv.gaussian, VGAM::pinv.gaussian),
    H = cumhazardify(VGAM::pinv.gaussian),
    fullname='inverse_gaussian_location'
  )

  flexsurv_invgaus <- list(
    name='inv.gaussian',
    pars=c('mu','lambda'), # mean, scale?
    location='lambda',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t)/10, 2)},
    d=VGAM::dinv.gaussian,
    p=VGAM::pinv.gaussian,
    h = hazardify(VGAM::dinv.gaussian, VGAM::pinv.gaussian),
    H = cumhazardify(VGAM::pinv.gaussian),
    fullname='inverse_gaussian'
  )

  ### Laplace
  flexsurv_laplace <- list(
    name='laplace',
    pars=c('mu','sigma'), # location, scale
    location='sigma',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(mean(t), sqrt(stats::var(t)/2))},
    d=extraDistr::dlaplace,
    p=extraDistr::plaplace,
    h = hazardify(extraDistr::dlaplace, extraDistr::plaplace),
    H = cumhazardify(extraDistr::plaplace),
    fullname='laplace'
  )

  ### Levy
  flexsurv_levy <- list(
    name='levy',
    pars=c('location','scale'), # location, scale
    location='scale',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(0, stats::median(t))},
    d=VGAM::dlevy,
    p=VGAM::plevy,
    h = hazardify(VGAM::dlevy, VGAM::plevy),
    H = cumhazardify(VGAM::plevy),
    fullname='levy'
  )

  ### Log cauchy
  flexsurv_logcauchy_loc <- list(
    name='logcauchy',
    pars=c('mu','sigma'), # location, scale
    location='mu',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(1,1)},
    d=dlogcauchy,
    p=plogcauchy,
    h = hazardify(dlogcauchy, plogcauchy),
    H = cumhazardify(plogcauchy),
    fullname='log_cauchy_location'
  )

  flexsurv_logcauchy <- list(
    name='logcauchy',
    pars=c('mu','sigma'),
    location='sigma',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(1,1)},
    d=dlogcauchy,
    p=plogcauchy,
    h = hazardify(dlogcauchy, plogcauchy),
    H = cumhazardify(plogcauchy),
    fullname='log_cauchy'
  )

  ### Lomax
  flexsurv_lomax <- list(
    name='lomax',
    pars=c('lambda','kappa'), # scale, shape?
    location='lambda',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t),1)},
    d=extraDistr::dlomax,
    p=extraDistr::plomax,
    h = hazardify(extraDistr::dlomax, extraDistr::plomax),
    H = cumhazardify(extraDistr::plomax),
    fullname='lomax'
  )

  ### Lomax shape parameter
  flexsurv_lomax_shape <- list(
    name='lomax',
    pars=c('lambda','kappa'), # scale, shape?
    location='kappa',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t),1)},
    d=extraDistr::dlomax,
    p=extraDistr::plomax,
    h = hazardify(extraDistr::dlomax, extraDistr::plomax),
    H = cumhazardify(extraDistr::plomax),
    fullname='lomax_shape'
  )

  ### nakagami
  flexsurv_nakagami <- list(
    name='naka',
    pars=c('scale','shape'),
    location='scale',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(mean(t^2),(mean(t^2)^2)/stats::var(t^2))},
    d=VGAM::dnaka,
    p=VGAM::pnaka,
    h=hazardify(VGAM::dnaka, VGAM::pnaka),
    H=cumhazardify(VGAM::pnaka),
    fullname='nakagami'
  )

  ### PARETO I , there's a few different versions of Pareto, all of which have almost never worked in testing
  flexsurv_pareto1 <- list(
    name='paretoI',
    pars=c('scale','shape'),
    location='scale',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(min(t)*.95,1)},
    d=VGAM::dparetoI,
    p=VGAM::pparetoI,
    h = hazardify(VGAM::dparetoI, VGAM::pparetoI),
    H = cumhazardify(VGAM::pparetoI),
    fullname='pareto_type_1'
  )

  ### Pareto II
  flexsurv_pareto2 <- list(
    name='paretoII',
    pars=c('location','scale','shape'),
    location='scale',
    transforms=c(identity,log,log),
    inv.transforms=c(identity,exp,exp),
    inits=function(t){c(0,min(t)*.95,1)},
    d=VGAM::dparetoII,
    p=VGAM::pparetoII,
    h = hazardify(VGAM::dparetoII, VGAM::pparetoII),
    H = cumhazardify(VGAM::pparetoII),
    fullname='pareto_type_2'
  )

  ### Pareto III
  flexsurv_pareto3 <- list(
    name='paretoIII',
    pars=c('location','scale','inequality'),
    location='scale',
    transforms=c(identity,log,log),
    inv.transforms=c(identity,exp,exp),
    inits=function(t){c(0,min(t)*.95,1)},
    d=VGAM::dparetoIII,
    p=VGAM::pparetoIII,
    h = hazardify(VGAM::dparetoIII, VGAM::pparetoIII),
    H = cumhazardify(VGAM::pparetoIII),
    fullname='pareto_type_3'
  )

  ### Pareto IV
  flexsurv_pareto4 <- list(
    name='paretoIV',
    pars=c('location','scale','inequality','shape'),
    location='scale',
    transforms=c(identity,log,log,log),
    inv.transforms=c(identity,exp,exp,exp),
    inits=function(t){c(0,min(t)*95,1,2)},
    d=VGAM::dparetoIV,
    p=VGAM::pparetoIV,
    h = hazardify(VGAM::dparetoIV, VGAM::pparetoIV),
    H = cumhazardify(VGAM::pparetoIV),
    fullname='pareto_type_4'
  )

  ### Feller-Pareto (basically Pareto type five)
  flexsurv_fpareto <- list(
    name='fpareto',
    pars=c('min','shape1','shape2','shape3','rate'),
    location='rate',
    transforms=c(identity,log,log,log,log),
    inv.transforms=c(identity,exp,exp,exp,exp),
    inits=function(t){c(0,1,1,1,1)},
    d=actuar::dfpareto,
    p=actuar::pfpareto,
    h = hazardify(actuar::dfpareto, actuar::pfpareto),
    H = cumhazardify(actuar::pfpareto),
    fullname='feller_pareto'
  )

  ### Truncated Pareto
  flexsurv_tpareto <- list(
    name='truncpareto',
    pars=c('lower','upper','shape'),
    location='shape',
    transforms=c(identity,identity, log),
    inv.transforms=c(identity, identity, exp),
    inits=function(t){c(min(t)-.1,max(t)+1, 1)},
    d=VGAM::dtruncpareto,
    p=VGAM::ptruncpareto,
    h = hazardify(VGAM::dtruncpareto, VGAM::ptruncpareto),
    H = cumhazardify(VGAM::ptruncpareto),
    fullname='truncated_pareto'
  )

  ### Rayleigh
  flexsurv_rayleigh <- list(
    name='rayleigh',
    pars=c('sigma'),
    location='sigma',
    transforms=c(log),
    inv.transforms=c(exp),
    inits=function(t){c(stats::median(t)/sqrt(2*log(2)))},
    d=extraDistr::drayleigh,
    p=extraDistr::prayleigh,
    h = hazardify(extraDistr::drayleigh, extraDistr::prayleigh),
    H = cumhazardify(extraDistr::prayleigh),
    fullname='rayleigh'
  )

  ### Rice
  flexsurv_rice_loc <- list(
    name='rice',
    pars=c('sigma','vee'), # scale location?
    location='vee',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t), stats::sd(t))},
    d=VGAM::drice,
    p=VGAM::price,
    h = hazardify(VGAM::drice, VGAM::price),
    H = cumhazardify(VGAM::price),
    fullname='rice_location'
  )

  flexsurv_rice <- list(
    name='rice',
    pars=c('sigma','vee'), # scale, location?
    location='sigma',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t), stats::sd(t))},
    d=VGAM::drice,
    p=VGAM::price,
    h = hazardify(VGAM::drice, VGAM::price),
    H = cumhazardify(VGAM::price),
    fullname='rice'
  )

  ### Shifted Gompertz
  flexsurv_sgomp <- list(
    name='sgomp',
    pars=c('b','eta'),
    location='b',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(1/stats::median(t),stats::median(t))},
    d=extraDistr::dsgomp,
    p=extraDistr::psgomp,
    h = hazardify(extraDistr::dsgomp, extraDistr::psgomp),
    H = cumhazardify(extraDistr::psgomp),
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
    h = hazardify(VGAM::dgumbelII, VGAM::pgumbelII),
    H = cumhazardify(VGAM::pgumbelII),
    fullname='gumbel_type_2'
  )

  # hyper tabastic
  flexsurv_hypertab_a <- list(
    name='hypertab',
    pars= c('a','b'),
    location= 'a',
    transforms= c(log,log),
    inv.transforms= c(exp,exp),
    inits= function(t){c(1,.1)},
    d=dhypertab,
    p=phypertab,
    h = hazardify(dhypertab, phypertab),
    H = cumhazardify(phypertab),
    fullname='hypertabastic'
  )

  flexsurv_hypertab_b <- list(
    name='hypertab',
    pars= c('a','b'),
    location= 'b',
    transforms= c(log,log),
    inv.transforms= c(exp,exp),
    inits= function(t){c(1,.1)},
    d=dhypertab,
    p=phypertab,
    h = hazardify(dhypertab, phypertab),
    H = cumhazardify(phypertab),
    fullname='hypertabastic'
  )

  ### Lindley distribution
  flexsurv_lindley <- list(
    name= 'lind',
    pars=c('theta'),
    location='theta',
    transforms=c(log),
    inv.transforms=c(exp),
    inits= function(t){c(1/mean(t))},
    d=VGAM::dlind,
    p=VGAM::plind,
    h = hazardify(VGAM::dlind, VGAM::plind),
    H = cumhazardify(VGAM::plind),
    fullname='lindley'
  )

  ### inverse Lindley
  flexsurv_invlind <- list(
    name= 'invlind',
    pars=c('theta'),
    location='theta',
    transforms=c(log),
    inv.transforms=c(exp),
    inits= function(t){c(1)},
    d=dinvlind,
    p=pinvlind,
    h = hazardify(dinvlind, pinvlind),
    H = cumhazardify(pinvlind),
    fullname='inverse_lindley'
  )

  dist_list <-list(
        betaprime               = flexsurv_betapr,
        fatigue                 = flexsurv_fatigue,
        fatigue_location        = flexsurv_fatigue_loc,
        fatigue_shape           = flexsurv_fatigue_shape,
        fatigue_shape_location  = flexsurv_fatigue_shape_loc,
        sinmad                  = flexsurv_sinmad,
        cauchy                  = flexsurv_cauchy,
        cauchy_location         = flexsurv_cauchy_loc,
        chisq                   = flexsurv_chisq,
        noncen_chisq            = flexsurv_ncchisq,
        dagum                   = flexsurv_dagum,
        exp_log                 = flexsurv_explog,
        extr_val                = flexsurv_extval,
        f                       = flexsurv_f,
        noncen_f                = flexsurv_ncf,
        fold_norm               = flexsurv_fnorm,
        fold_norm_location      = flexsurv_fnorm_loc,
        frechet                 = flexsurv_frech,
        gamgomp                 = flexsurv_gamgomp,
        gumbel                  = flexsurv_gumbel,
        inv_chisq               = flexsurv_invchisq,
        scale_inv_chisq         = flexsurv_sinvchisq,
        inv_gamma               = flexsurv_invgam,
        inv_gaussian            = flexsurv_invgaus,
        inv_gaussian_location   = flexsurv_invgaus_loc,
        laplace                 = flexsurv_laplace,
        levy                    = flexsurv_levy,
        log_cauchy              = flexsurv_logcauchy,
        log_cauchy_location     = flexsurv_logcauchy_loc,
        lomax                   = flexsurv_lomax,
        lomax_shape             = flexsurv_lomax_shape,
        nakagami                = flexsurv_nakagami,
        pareto_1                = flexsurv_pareto1,
        pareto_2                = flexsurv_pareto2,
        pareto_3                = flexsurv_pareto3,
        pareto_4                = flexsurv_pareto4,
        feller_pareto           = flexsurv_fpareto,
        trunc_pareto            = flexsurv_tpareto,
        rayleigh                = flexsurv_rayleigh,
        rice                    = flexsurv_rice,
        rice_location           = flexsurv_rice_loc,
        shift_gomp              = flexsurv_sgomp,
        gumbel_T2               = flexsurv_gumbelII,
        hypertab_a              = flexsurv_hypertab_a,
        hypertab_b              = flexsurv_hypertab_b,
        lindley                 = flexsurv_lindley,
        inv_lind                = flexsurv_invlind

        ### Removed for the time being
        # flexsurv_erlang,
      )

  dist_list <- c(flexsurv::flexsurv.dists, dist_list)
  return(dist_list)
}
