source('R/flexsurvdist.R')

#' Shotgun
surv_shotgun <- function(formula, data=NA, skip=c('default'), dump_models =F, progress=T, warn=F){

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
  tictoc::tic("Kapow")

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
      tictoc::tic(current_dist)
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

    tictoc::toc(quiet =T)$callback_msg %>% message()
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
  tictoc::toc(quiet=T)$callback_msg %>% message()
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
