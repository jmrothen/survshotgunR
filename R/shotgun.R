## This script contains the primary survival-shotgun function
##



#' Shotgun
#'
#' @param formula Formula. Should be a survival formula, with a Surv object on the left hand side.
#' @param data If your formula needs a dataset, provide that here.
#' @param skip Vector. If you want to skip any specific models, you can add their names here. By default, some of the repetitive or incredibly niche models are skipped.
#' @param dump_models Logical. If TRUE, each successful model created will be loaded into memory as fssg_<model_name>
#' @param progress Logical. Want progress updates?
#' @param warn Logical. If TRUE, also prints any warnings that appear.
#' @returns Data frame summarizing each model, and some general goodness of fit measures.
#'
#' @examples
#' library(survival)
#' surv_shotgun(Surv(time, status)~1, data=aml, dump_models=TRUE, warn = TRUE)
#'
#'
#' @export
surv_shotgun <- function(
    formula,
    data=NA,
    skip=c('default'),
    dump_models =F,
    progress=T,
    warn=F
){
  # The default shotgun list will exclude the following. Comments describe why
  if(length(skip)==1 & skip[1]=='default'){
    skip <- c(
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
    dist_ran = T,
    aic = 1,
    bic = 1,
    loglik = 1
  )[-1,]


  # iterate through each distribution, creating the model if possible and storing results
  iter <- 1
  for(i in dist_list){

    # skip distribution if in our skip list
    if(i$name %in% skip | dplyr::coalesce(i$fullname, i$name) %in% skip | names(dist_list)[iter] %in% skip){
      iter <- iter+1
      next
    }

    current_dist <- dplyr::coalesce(i$fullname, i$name)

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

        # if we were supplied dataset, include that here, otherwise we just use formula
        if(data_req){
          flexsurv::flexsurvreg(formula, dist=i, data=data, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
        }else{
          flexsurv::flexsurvreg(formula, dist=i, dfns=list(d=i$d, p=i$p)) %>% suppressMessages() -> current_model
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
        dist_name = current_dist, dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
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
  # quick tests
  surv_shotgun(survival::Surv(time, status)~1, data=survival::aml, dump_models=T, warn = T)
  surv_shotgun(survival::Surv(time,status) ~1, data = survival::cancer, dump_models = T, warn=T)


  surv_shotgun(survival::Surv(time,status)~ x, data=survival::aml, dump_models=T, warn = T)

  surv_shotgun(survival::Surv(time,status)~ inst, data = survival::cancer, dump_models = T, warn=T)
  surv_shotgun(survival::Surv(time,status)~ sex, data = survival::cancer, dump_models = T, warn=T)
  surv_shotgun(survival::Surv(time,status)~ ph.ecog, data = survival::cancer, dump_models = T, warn=T)

}
