#### Custom distribution shell
# if(F){
#   custom_model <- list(
#     name = '',               # this should be the shortened name that appears after d/p in the distribution functions i.e. "norm" in dnorm
#     pars = c('','',''),      # a vector of each parameter name, exactly as they are stated in the d/p functions i.e c('mean','sd') for norm
#     location= c(),           # if there is a location or scale parameter, you can place it here. It's used in the scaling on variable in more complicated models
#     transforms = c(),        # a vector of functions, one for each parameter, that transform the respective parameter into the real line. If param1 must be >0, log(param1) would be real
#     inv.transforms=c(),      # a vector of functions which transform the real values back to the original scale. exp(log(param1) returns to the original range of param1
#     inits = function(t){},   # a declared function which takes the time variable, and outputs initial estimates for the parameters. Ideally, estimates are based off of t somehow, but can be simple ints
#     p = pdist,               # by default, flexsurv looks for the d/p functions using the <name> object above, but the shotgun function will use these declared p/d's instead
#     d = ddist,               #
#     fullname = ''            # this longer name is used in shotgun to better label the distribution in outputs and console tracking
#   )
# }

#' Compiles list of available distributions
#'
#' @returns a list of all possible distributions
#'
#' @export
shotgun_dist_list <- function(){

  # create a list item for each currently available distribution

  ### Beta Prime
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

  ### Birnbaum - Saunders, specified as fatigue in extraDistr
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

  ### Singh Madalla / burr t12
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

  ### CAUCHY
  flexsurv_cauchy <- list(
    name='cauchy',
    pars= c('location','scale'),
    location= c('location'),
    transforms = c(identity, log),
    inv.transforms = c(identity, exp),
    inits = function(t){c(stats::median(t), stats::mad(t)+1e-8)},
    d= stats::dcauchy,
    p= stats::pcauchy,
    fullname='cauchy'
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
    fullname='chi_squared'
  )

  ### non central chi squared
  flexsurv_ncchisq <- list(
    name='chisq',
    pars=c('df','ncp'),
    location=c('ncp'),
    transforms=c(log, log),
    inv.transforms=c(exp,exp),
    inits = function(t){c(mean(t), 1)},
    d = stats::dchisq,
    p = stats::pchisq,
    fullname='non_central_chi_squared'
  )

  ### Dagum
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
  logit <- function(p){log(p/(1-p))}
  inv_logit <- function(x){return(1/(1+exp(-x)))}
  flexsurv_explog <- list(
    name = 'explog',
    pars = c('scale','shape'),
    location = 'scale',
    transforms = c(log, logit),
    inv.transforms = c(exp, inv_logit),
    inits=function(t){c(stats::median(t),.5)},
    p=VGAM::pexplog,
    d=VGAM::dexplog,
    fullname='exponential_logarithmic'
  )

  ### Extreme value (generalized)
  flexsurv_extval <- list(
    name='gev',
    pars = c('mu','sigma','xi'),
    location='mu',
    transforms= c(identity, log, log),  # trying to make xi use log now, as the bounding is bad if <=0
    inv.transforms = c(identity, exp, exp),
    inits=function(t){c(1,1,1)},
    d=extraDistr::dgev,
    p=extraDistr::pgev,
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
    fullname='f'
  )

  ### noncentral F
  flexsurv_ncf <- list(
    name='f',
    pars = c('df1','df2', 'ncp'),
    location = ('ncp'),
    transforms = c(log,log,log),
    inv.transforms = c(exp,exp,exp),
    inits=function(t){c(stats::median(t),2,0.1)},
    d=stats::df,
    p=stats::pf,
    fullname='noncentral_f'
  )

  ### folded normal?
  flexsurv_fnorm <- list(
    name= 'foldnorm',
    pars = c('mean','sd','a1','a2'),
    location = c('mean'),
    transforms = c(identity, log, identity, identity),
    inv.transforms = c(identity, exp, identity, identity),
    inits = function(t){c(mean(t),stats::sd(t),1,1)},
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
    inits=function(t){c(log(2)/stats::median(t),1,stats::IQR(t)/stats::median(t))},
    d=dgamgomp,
    p=pgamgomp,
    fullname='gamma_gompertz'
  )

  ### Gumbel (basically log weibull)
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
    fullname='scaled_inverse_chi_squared'
  )


  ### Inverse Gamma
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

  ### Inverse Gaussian
  flexsurv_invgaus <- list(
    name='inv.gaussian',
    pars=c('mu','lambda'),
    location='mu',
    transforms=c(log,log),
    inv.transforms=c(exp,exp),
    inits=function(t){c(stats::median(t)/10, 2)},
    d=VGAM::dinv.gaussian,
    p=VGAM::pinv.gaussian,
    fullname='inverse_gaussian'
  )

  ### Laplace
  flexsurv_laplace <- list(
    name='laplace',
    pars=c('mu','sigma'),
    location='mu',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(mean(t), sqrt(stats::var(t)/2))},
    d=extraDistr::dlaplace,
    p=extraDistr::plaplace,
    fullname='laplace'
  )

  ### Levy
  flexsurv_levy <- list(
    name='levy',
    pars=c('location','scale'),
    location='location',
    transforms=c(identity, log),
    inv.transforms=c(identity, exp),
    inits=function(t){c(0 - 2e-4, stats::median(t))},
    d=VGAM::dlevy,
    p=VGAM::plevy,
    fullname='levy'
  )

  ### Log cauchy
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
    inits=function(t){c(stats::median(t),1)},
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
    inits=function(t){c(mean(t^2),(mean(t^2)^2)/stats::var(t^2))},
    d=VGAM::dnaka,
    p=VGAM::pnaka,
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
    fullname='pareto_type_1'
  )

  ### Pareto II
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

  ### Pareto III
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

  ### Pareto IV
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

  ### Feller-Pareto (basically Pareto type five)
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

  ### Truncated Pareto
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
    inits=function(t){c(stats::median(t)/sqrt(2*log(2)))},
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
    inits=function(t){c(stats::median(t), stats::sd(t))},
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
    inits=function(t){c(1/stats::median(t),stats::median(t))},
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
    fullname='inverse_lindley'
  )

  dist_list <-list(
        betaprime       = flexsurv_betapr,
        fatigue         = flexsurv_fatigue,
        sinmad          = flexsurv_sinmad,
        cauchy          = flexsurv_cauchy,
        chisq           = flexsurv_chisq,
        noncen_chisq    = flexsurv_ncchisq,
        dagum           = flexsurv_dagum,
        exp_log         = flexsurv_explog,
        extr_val        = flexsurv_extval,
        f               = flexsurv_f,
        noncen_f        = flexsurv_ncf,
        fold_norm       = flexsurv_fnorm,
        frechet         = flexsurv_frech,
        gamgomp         = flexsurv_gamgomp,
        gumbel          = flexsurv_gumbel,
        inv_chisq       = flexsurv_invchisq,
        scale_inv_chisq = flexsurv_sinvchisq,
        inv_gamma       = flexsurv_invgam,
        inv_gaussian    = flexsurv_invgaus,
        inv_lap         = flexsurv_laplace,
        levy            = flexsurv_levy,
        log_cauchy      = flexsurv_logcauchy,
        lomax           = flexsurv_lomax,
        nakagami        = flexsurv_nakagami,
        pareto_1        = flexsurv_pareto1,
        pareto_2        = flexsurv_pareto2,
        pareto_3        = flexsurv_pareto3,
        pareto_4        = flexsurv_pareto4,
        feller_pareto   = flexsurv_fpareto,
        trunc_pareto    = flexsurv_tpareto,
        rayleigh        = flexsurv_rayleigh,
        rice            = flexsurv_rice,
        shift_gomp      = flexsurv_sgomp,
        gumbel_T2       = flexsurv_gumbelII,
        hypertab        = flexsurv_hypertab,
        lindley         = flexsurv_lindley,
        inv_lind        = flexsurv_invlind

        ### Removed for the time being
        # flexsurv_erlang,
      )

  dist_list <- c(flexsurv::flexsurv.dists, dist_list)
  return(dist_list)
}
