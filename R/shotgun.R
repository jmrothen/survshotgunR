## This script contains the primary survival-shotgun function



#' Shotgun
#'
#' @param formula Formula. Should be a survival formula, with a Surv object on the left hand side.
#' @param data If your formula needs a dataset, provide that here.
#' @param skip Vector. If you want to skip any specific models, you can add their names here. By default, some of the repetitive or incredibly niche models are skipped.
#' @param dump_models Logical. If TRUE, each successful model created will be loaded into memory as fssg_<model_name>
#' @param progress Logical. Want progress updates?
#' @param warn Logical. If TRUE, also prints any warnings that appear.
#' @param spline Vector. Should include 'rp' for Royston-Parmar natural cubic spline. Can also include 'wy' for Wang-Yan alternative natural cubic spline. The Wang-Yan version requires the package 'splines2ns'.
#' @param max_knots Integer. Specifies the maximum number of knots considered in spline models.
#' @param coxph Logical. If TRUE, calculates a Cox Proportional Hazard Model as well. Please note that is not recommended to directly compare the AIC/BIC/LogLik of Cox models to Parametric models
#' @returns Data frame summarizing each model, and some general goodness of fit measures.
#'
#' @examples
#' library(survival)
#' surv_shotgun(Surv(time, status)~1, data=aml, dump_models=TRUE, warn = TRUE)
#'
#' @export
surv_shotgun <- function(
    formula,
    data=NA,
    skip=c('default'),
    dump_models=F,
    progress=T,
    warn=F,
    spline=c('rp'),
    max_knots=5,
    coxph=F
){
  # The default shotgun list will exclude the following. Comments describe why
  if(length(skip)==1 & skip[1]=='default'){
    skip <- c(
      'weibullPH',                  # identical to other weibull
      'genf.orig',                  # running genf instead
      'gengamma.orig',              # running gengamma instead
      'erlang',                     # erlang is temporarily removed anyways
      'truncpareto',                # have yet to get this to work
      'chisq',                      # very rarely useful
      'non_central_chi_squared',    # very rarely useful
      'exponential',                # flexsurv's distribution list has 'exp','exponential', which are identical in practice. So we toss this one
      'lognormal'                   # flexsurv's distribution list has 'lognormal','loggaussian', which are identical in practice. So we toss this one
    )
  }

  message("Loading shotgun...")
  tictoc::tic("Kapow") # overall timer

  # parsing / verifying the provided formula
  vars <- all.vars(formula)
  in_global <- sapply(vars, exists, envir=.GlobalEnv)
  data_req <- F

  # if data is provided (as a dataframe)
  if(is.data.frame(data)){

    # if the formula doesnt fit the dataframe
    if(!all(vars %in% names(data))){

      # check if we have the formula variables in global
      if(!all(in_global)){
        stop("Formula does not fit data or global environment.")
      }else{
        message("Formula did not fit with dataset, so global variables are being used instead")
      }
    }

    # if formula DOES fit data frame, mark that we'll use the dataframe
    else{
      data_req <-T
    }
  }

  # if we don't have a data frame
  else{

    # check if we have formula variables in global
    if(!all(in_global)){
      stop("Formula does not fit global environment.")
    }
  }

  # grab list of all available distributions
  dist_list <- shotgun_dist_list()

  # initialize a quick dataframe which we'll update with each iteration
  dist_summary <- data.frame(
    dist_name = '1',
    dist_source= '2',
    dist_ran = T,
    aic = 1,
    bic = 1,
    loglik = 1
  )[-1,]

  # iterate through each distribution, creating the model if possible and storing results
  iter <- 1
  for(i in dist_list){

    current_dist <- names(dist_list)[iter]
    custom_indicator <- !(current_dist %in% names(flexsurv::flexsurv.dists))
    current_source <- ifelse(custom_indicator, 'survshotgun','flexsurv')

    # skip distribution if in our skip list
    if(i$name %in% skip | dplyr::coalesce(i$fullname, i$name) %in% skip | names(dist_list)[iter] %in% skip){
      iter <- iter+1
      next
    }

    # optional progress tracking chunk
    if(progress){
      message('------------------------------------------------------------------------')
      message(current_dist)
      tictoc::tic(current_dist)
    }

    # reset iteration level variables
    dist_success<- F
    current_aic <- NA
    current_bic <- NA
    current_ll <- NA

    # attempt to perform regression
    tryCatch({

      # additional level of obfuscation here to allow for us to continue on warnings
      withCallingHandlers({

        # if working with a native distribution, we'll not specify the <dfns> argument
        if(!custom_indicator)
          if(data_req){
            flexsurv::flexsurvreg(formula, dist=current_dist, data=data) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvreg(formula, dist=current_dist) %>% suppressMessages() -> current_model
          }

        # for custom distributions, we specify the DFNs
        else{
          if(data_req){
            flexsurv::flexsurvreg(formula, dist=i, data=data, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvreg(formula, dist=i, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
          }
        }

        # if model succeeds, collect information
        current_aic <- stats::AIC(current_model)
        current_bic <- stats::BIC(current_model)
        current_ll <- current_model$loglik
        dist_success<-T
      },

      # warnings are very common in flexsurv / optim, so if we get one, we conditionally print and continue
      warning=function(w){
        if(warn){message(paste('Warning in', current_dist,'model :',w))}
        invokeRestart('muffleWarning')
      })
    },

    # If the TryCatch errors, print out the error that occurred and continue
    error=function(e){
      message(paste("Error in", current_dist,"model :",e))
    })

    # add current iteration's results to our progress dataframe
    dist_summary %>%
      dplyr::add_row(
        dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
      ) -> dist_summary

    # if we are model dumping, we assign the model to the global environment with name fssg_<model name>
    if(dump_models & dist_success){
      assign(paste('fssg_',current_dist,sep='',collapse=''), current_model, envir = .GlobalEnv) # fssg = flex surv shot gun
    }

    # close the iteration tracker and print message
    if(progress){
      tictoc::toc(quiet =T)$callback_msg %>% message()
    }

    iter = iter+1
  }

  # Spline section
  if('rp' %in% spline  | 'wy' %in% spline){

    # RP works on it's own, but if we also want to try the 'natural cubic spline', we need splines2
    if('wy' %in% spline){requireNamespace('splines2')}

    # optional progress tracking chunk
    if(progress){
      message('------------------------------------------------------------------------')
      message("Spline")
      tictoc::tic(current_dist)
    }

    n_reps <- sum(('rp' %in% spline), ('wy' %in% spline))

    kvec <- c(rep(rep(1:max_knots), 3*n_reps))
    svec <- c(rep(c(rep('hazard',max_knots), rep('odds',max_knots), rep('normal',max_knots)),n_reps))
    if(n_reps==2){
      mvec <- c(rep("rp",max_knots*3), rep('splines2ns',max_knots*3))
    }else{
      if('rp' %in% spline){
        mvec <- rep("rp",max_knots*3)
      }else{
        mvec <-  rep('splines2ns',max_knots*3)
      }
    }
    current_source <- 'flexsurv'

    for(s in 1:length(kvec)){

      # iteration level name
      current_dist <- paste('spline',ifelse(mvec[s]=='rp','rp','wy'), svec[s], kvec[s], sep='_')

      # reset iteration level variables
      dist_success<- F
      current_aic <- NA
      current_bic <- NA
      current_ll <- NA

      tryCatch({

        # additional level of obfuscation here to allow for us to continue on warnings
        withCallingHandlers({

          # cycle through our spline options
          if(data_req){
            flexsurv::flexsurvspline(formula, data=data, k=kvec[s], scale=svec[s], spline=mvec[s]) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvspline(formula, k=kvec[s], scale=svec[s], spline=mvec[s]) %>% suppressMessages() -> current_model
          }

          # if model succeeds, collect information
          current_aic <- stats::AIC(current_model)
          current_bic <- stats::BIC(current_model)
          current_ll <- current_model$loglik
          dist_success<-T
        },

        # warnings are very common in flexsurv / optim, so if we get one, we conditionally print and continue
        warning=function(w){
          if(warn){message(paste('Warning in', current_dist,'model :',w))}
          invokeRestart('muffleWarning')
        })
      },

      # If the TryCatch errors, print out the error that occurred and continue
      error=function(e){
        message(paste("Error in", current_dist,"model :",e))
      })

      # If spline succeeds, we add, otherwise don't (prevents clutter)
      if(dist_success){
        dist_summary %>%
          dplyr::add_row(
            dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
          ) -> dist_summary

        # if we are model dumping, we assign the model to the global environment with name fssg_<model name>
        if(dump_models & dist_success){
          assign(paste('fssg_',current_dist,sep='',collapse=''), current_model, envir = .GlobalEnv) # fssg = flex surv shot gun
        }
      }
    }

    # close the spline tracker and print message
    if(progress){
      tictoc::toc(quiet =T)$callback_msg %>% message()
    }

  }

  # simple Coxph section
  if(coxph){

    current_dist <- 'coxph'
    current_source <- 'survival'

    # optional progress tracking chunk
    if(progress){
      message('------------------------------------------------------------------------')
      message(current_dist)
      tictoc::tic(current_dist)
    }

    # reset iteration level variables
    dist_success<- F
    current_aic <- NA
    current_bic <- NA
    current_ll <- NA

    tryCatch({

      # additional level of obfuscation here to allow for us to continue on warnings
      withCallingHandlers({

        # cycle through our spline options
        if(data_req){
          survival::coxph(formula, data=data) %>% suppressMessages() -> current_model
        }else{
          survival::coxph(formula) %>% suppressMessages() -> current_model
        }

        # if model succeeds, collect information
        current_aic <- stats::AIC(current_model)
        current_bic <- stats::BIC(current_model)
        current_ll <- current_model$loglik
        dist_success<-T
      },

      # warnings are very common in flexsurv / optim, so if we get one, we conditionally print and continue
      warning=function(w){
        if(warn){message(paste('Warning in', current_dist,'model :',w))}
        invokeRestart('muffleWarning')
      })
    },

    # If the TryCatch errors, print out the error that occurred and continue
    error=function(e){
      message(paste("Error in", current_dist,"model :",e))
    })

    # If spline succeeds, we add, otherwise don't (prevents clutter)
    dist_summary %>%
      dplyr::add_row(
        dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
      ) -> dist_summary

      # if we are model dumping, we assign the model to the global environment with name fssg_<model name>
    if(dump_models & dist_success){
      assign(paste('fssg_',current_dist,sep='',collapse=''), current_model, envir = .GlobalEnv) # fssg = flex surv shot gun
    }

    # close the cox tracker and print message
    if(progress){
      tictoc::toc(quiet =T)$callback_msg %>% message()
    }
  }

  # garbage clean
  gc(verbose = F)

  message('------------------------------------------------------------------------')
  message('Final Summary!')

  # identify models with lowest AIC, BIC, or highest LogLik  <PLANNING TO ADD MORE STATS>
  dist_summary %>% dplyr::filter(!is.na(aic)) %>% dplyr::filter(aic== min(aic))       %>% dplyr::select(dist_name) %>% dplyr::pull() -> best_aic
  dist_summary %>% dplyr::filter(!is.na(aic)) %>% dplyr::filter(bic== min(bic))       %>% dplyr::select(dist_name) %>% dplyr::pull() -> best_bic
  dist_summary %>% dplyr::filter(!is.na(aic)) %>% dplyr::filter(loglik== max(loglik)) %>% dplyr::select(dist_name) %>% dplyr::pull() -> best_ll

  message(paste('Model with best AIC:',best_aic))
  message(paste('Model with best BIC:',best_bic))
  message(paste('Model with best Log Likelihood:',best_ll))

  # end overall timer (aka fire the shotgun)
  tictoc::toc(quiet=T)$callback_msg %>% message()

  # output is a data-frame, one row for each attempted model
  dist_summary %>% dplyr::arrange(aic, bic, desc(loglik)) %>% return()
}


### raw testing chunk from legacy testing
if(F){

  ###
  # coxph(survival::Surv(time, status) ~ 1, data = survival::aml)
  # coxph(survival::Surv(time, status) ~ factor(x), data = survival::aml)

  ### Verify flexsurvreg and survreg work the same for us (same LL, AIC, etc)
  # survreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')
  # flexsurvreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')

  # simple model tests
  surv_shotgun(survival::Surv(time, status) ~ 1, data = survival::aml,    dump_models = T, warn = F, spline=c('rp','wy'))
  surv_shotgun(survival::Surv(time, status) ~ 1, data = survival::cancer, dump_models = T, warn = F)

  # single variable models
  surv_shotgun(survival::Surv(time, status) ~ factor(x),   data = survival::aml,    dump_models = T, warn = F)
  surv_shotgun(survival::Surv(time, status) ~ factor(sex), data = survival::cancer, dump_models = T, warn = F)
}
