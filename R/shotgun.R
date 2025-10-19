## This script contains the primary survival-shotgun function



#' Shotgun
#'
#' @param formula Formula. Should be a survival formula, with a Surv object on the left hand side.
#' @param data If your formula needs a dataset, provide that here.
#' @param models Vector of strings. If you only want to run specific models, specify them here by their list name in shotgun_dist_list.
#' @param skip Vector. If you want to skip any specific models, you can add their names here. By default, some of the repetitive or incredibly niche models are skipped.
#' @param opt_method String. By default, 'BFGS' is used in flexsurvreg, however some distributions appreciate the more flexible 'Nelder-Mead' method. This is passed to the "optim" function as method = opt_method.
#' @param spline Vector. Should include 'rp' for Royston-Parmar natural cubic spline. Can also include 'wy' for Wang-Yan alternative natural cubic spline. The Wang-Yan version requires the package 'splines2ns'. If set to NA, then the spline step will be skipped.
#' @param max_knots Integer. Specifies the maximum number of knots considered in spline models.
#' @param dump_models Logical. If TRUE, each successful model will be placed into a list and returned by the function invisibly.
#' @param detailed Logical. If True, calculates a number of additional fit statistics for each model.
#' @param ibs Logical. If TRUE, calculate integrated brier score for each model. Please note that this greatly increases run time, and is not recommended for extremely large data.
#' @param progress Logical. Want progress updates?
#' @param warn Logical. If TRUE, also prints any warnings that appear.
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
    models=NA,
    skip=c('default'),
    opt_method = 'BFGS',
    spline=NA,   # c('rp', 'wy')
    max_knots=1,
    dump_models=TRUE,
    detailed=TRUE,
    ibs=FALSE,
    progress=TRUE,
    warn=FALSE
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
  }else{ # if we don't have a data frame

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

  # If we want the additional fit stats, initiate them as well
  if(detailed){
    dist_summary <- dplyr::mutate(
      dist_summary,
      iAUC=1, Cindex=1, Unos.C=1, Brier.Median=1, MAE=1, IAE=1, ISE=1, IBS=1
    )
  }

  # create a simple list item which we will be packing all of the models into for returning
  if(dump_models){
    working_model_list <- list()
  }


  # if specific models were specified, then we can filter to those
  if(!anyNA(models)){
    dist_list <- dist_list[names(dist_list) %in% models]
  }

  # iterate through each distribution, creating the model if possible and storing results
  iter <- 1
  for(i in dist_list){

    current_dist <- names(dist_list)[iter]
    custom_indicator <- !(current_dist %in% names(flexsurv::flexsurv.dists))
    current_source <- ifelse(custom_indicator, 'survshotgun','flexsurv')

    # skip distribution if in our skip list. This is done inside the loop to allow for broader matching of model names
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
    current_iAUC <- NA
    current_Cindex <- NA
    current_Unos.C <- NA
    current_brier <- NA
    current_mae <- NA
    current_iae <- NA
    current_ise <- NA
    current_ibs <- NA

    # attempt to perform regression
    tryCatch({

      # additional level of obfuscation here to allow for us to continue on warnings
      withCallingHandlers({

        # if working with a native distribution, we'll not specify the <dfns> argument
        if(!custom_indicator){
          if(data_req){
            flexsurv::flexsurvreg(formula, dist=current_dist, data=data, method=opt_method) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvreg(formula, dist=current_dist, method=opt_method) %>% suppressMessages() -> current_model
          }
        }else{ # for custom distributions, we specify the DFNs
          if(data_req){
            flexsurv::flexsurvreg(formula, dist=i, data=data, dfns=list(d=i$d, p=i$p), method=opt_method) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvreg(formula, dist=i, dfns=list(d=i$d, p=i$p), method=opt_method) %>% suppressMessages() -> current_model
          }
        }

        # if model succeeds, collect information
        current_aic <- stats::AIC(current_model)
        current_bic <- stats::BIC(current_model)
        current_ll <- current_model$loglik

        dist_success<-T

        # If we're interested in dumping the models, we pass the current model into the working model list
        if(dump_models){
          working_model_list[[current_dist]] <- current_model
        }


        if(detailed){
        ### should add the case to pase Surv functions of the form Surv(time1, time2, status), which would use length(formula[[2]])
          if(data_req){
            time_portion <-   dplyr::pull(data[c(as.character(formula[[2]][[2]]))])
            status_portion <- dplyr::pull(data[c(as.character(formula[[2]][[3]]))])
            Surv_object <- survival::Surv(time_portion, status_portion)
          }else{
            Surv_object <- survival::Surv(formula[[2]][[2]], formula[[2]][[3]])
          }
          fitstats <- get_fit_stats(Surv_object, model = current_model, ibs)

          current_iAUC <- fitstats$iAUC.Full
          current_Cindex <- fitstats$C.Index
          current_Unos.C <- fitstats$C.Index.Uno
          current_brier <- fitstats$Brier.Median
          current_mae <- fitstats$MAE
          current_iae <- fitstats$IAE.Full
          current_ise <- fitstats$ISE.Full
          current_ibs <- ifelse(ibs, fitstats$IBS.Full, NA)
        }
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
    if(detailed){
      dist_summary %>%
        dplyr::add_row(
          dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll,
          iAUC=current_iAUC, Cindex=current_Cindex, Unos.C=current_Unos.C, Brier.Median=current_brier, MAE=current_mae, IAE=current_iae, ISE=current_ise, IBS=current_ibs
        ) -> dist_summary
    }else{
      dist_summary %>%
        dplyr::add_row(
          dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
        ) -> dist_summary
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
    # if('wy' %in% spline){requireNamespace('splines2')}

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
      current_iAUC <- NA
      current_Cindex <- NA
      current_Unos.C <- NA
      current_brier <- NA
      current_mae <- NA
      current_iae <- NA
      current_ise <- NA
      current_ibs <- NA

      tryCatch({

        # additional level of obfuscation here to allow for us to continue on warnings
        withCallingHandlers({

          # cycle through our spline options
          if(data_req){
            flexsurv::flexsurvspline(formula, data=data, k=kvec[s], scale=svec[s], spline=mvec[s], method=opt_method) %>% suppressMessages() -> current_model
          }else{
            flexsurv::flexsurvspline(formula, k=kvec[s], scale=svec[s], spline=mvec[s], method=opt_method) %>% suppressMessages() -> current_model
          }

          # if model succeeds, collect information
          current_aic <- stats::AIC(current_model)
          current_bic <- stats::BIC(current_model)
          current_ll <- current_model$loglik

          dist_success<-T

          # If we're interested in dumping the models, we pass the current model into the working model list
          if(dump_models){
            working_model_list[[current_dist]] <- current_model
          }

          if(detailed){
            ### should add the case to pass Surv functions of the form Surv(time1, time2, status), which would use length(formula[[2]])
            if(data_req){
              time_portion <-   dplyr::pull(data[c(as.character(formula[[2]][[2]]))])
              status_portion <- dplyr::pull(data[c(as.character(formula[[2]][[3]]))])
              Surv_object <- survival::Surv(time_portion, status_portion)
            }else{
              Surv_object <- survival::Surv(formula[[2]][[2]], formula[[2]][[3]])
            }
            fitstats <- get_fit_stats(Surv_object, model = current_model, ibs)

            current_iAUC <- fitstats$iAUC.Full
            current_Cindex <- fitstats$C.Index
            current_Unos.C <- fitstats$C.Index.Uno
            current_brier <- fitstats$Brier.Median
            current_mae <- fitstats$MAE
            current_iae <- fitstats$IAE.Full
            current_ise <- fitstats$ISE.Full
            current_ibs <- ifelse(ibs, fitstats$IBS.Full, NA)
          }
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
        if(detailed){
          dist_summary %>%
            dplyr::add_row(
              dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll,
              iAUC=current_iAUC, Cindex=current_Cindex, Unos.C=current_Unos.C, Brier.Median=current_brier, MAE=current_mae, IAE=current_iae, ISE=current_ise, IBS=current_ibs
            ) -> dist_summary
        }else{
          dist_summary %>%
            dplyr::add_row(
              dist_name = current_dist, dist_source=current_source,  dist_ran=dist_success, aic=current_aic, bic=current_bic, loglik=current_ll
            ) -> dist_summary
        }
      }
    }

    # close the spline tracker and print message
    if(progress){
      tictoc::toc(quiet =T)$callback_msg %>% message()
    }

  }

  # garbage clean
  gc(verbose = F)

  message('------------------------------------------------------------------------')
  message('Final Summary!')

  ##### We want to aggregate the bests fit via statistics
  dplyr::mutate(
    dist_summary,
    best_aic = (aic== min(aic,na.rm=T)),  # lower = better
    best_bic = (bic== min(bic,na.rm=T)),  # lower = better
    best_loglik = (loglik== max(loglik,na.rm=T)) # greater = better
  ) -> dist_summary

  if(detailed){

    dplyr::mutate(
      dist_summary,
      best_iauc = (iAUC== max(iAUC,na.rm=T)), # greater = better
      best_cin = (Cindex== max(Cindex,na.rm=T)), # greater = better
      best_uno = (Unos.C== max(Unos.C,na.rm=T)), # greater = better
      best_bri = (Brier.Median== min(Brier.Median, na.rm=T)), # lower = better
      best_mae = (MAE== min(MAE,na.rm=T)), # lower = better
      best_iae = (IAE== min(IAE,na.rm=T)), # lower = better
      best_ise = (ISE== min(ISE,na.rm=T)) # lower = better
    ) -> dist_summary

    if(ibs){
      dplyr::mutate(
        dist_summary,
        best_ibs = (IBS== min(IBS,na.rm=T)), # lower = better
      ) -> dist_summary
    }
  }

  message(paste('Model with best AIC:', paste(dplyr::filter(dist_summary, best_aic==TRUE)$dist_name, collapse=', ')))
  message(paste('Model with best BIC:', paste(dplyr::filter(dist_summary, best_bic==TRUE)$dist_name, collapse=', ')))

  if(detailed){
    if(length(dplyr::filter(dist_summary, best_iauc=TRUE)$dist_name) >5){
      message(paste('Many Models tied for best iAUC'))
    }else{
      message(paste('Model with best iAUC:', paste(dplyr::filter(dist_summary, best_iauc==TRUE)$dist_name, collapse=', ')))
    }

    if(length(dplyr::filter(dist_summary, best_cin=TRUE)$dist_name) >5){
      message(paste('Many Models tied for best C-index'))
    }else{
      message(paste('Model with best C-index:', paste(dplyr::filter(dist_summary, best_cin==TRUE)$dist_name, collapse=', ')))
    }
  }


  if(ibs){
    message(paste('Model with best IBS:', paste(dplyr::filter(dist_summary, best_ibs==TRUE)$dist_name, collapse=', ')))
  }else{
     dplyr::select(dist_summary, -IBS) -> dist_summary
  }

  for(q in colnames(dist_summary)){
    if(grepl('best',q)){
      dist_summary %>% dplyr::select(-as.character(q)) -> dist_summary
    }
  }

  # end overall timer (aka fire the shotgun)
  tictoc::toc(quiet=T)$callback_msg %>% message()

  # output is a data-frame, one row for each attempted model
  if(detailed){
    if(ibs){
       dplyr::arrange(dist_summary, aic, bic, dplyr::desc(IBS), iAUC, Cindex) -> out
    }else{
       dplyr::arrange(dist_summary,aic, bic, iAUC, Cindex) -> out
    }
  }else{
    dplyr::arrange(dist_summary, aic, bic, dplyr::desc(loglik)) -> out
  }

  # always print the summary information
  print(out)

  # Return the summary, and a list of models if the option is provided
  if(dump_models){
    invisible(list(
      summary = out,
      models = working_model_list
    ))
  }else{
    invisible(out)
  }
}

