#' Collect fit statistics for a parametric survival model
#'
#' @param Surv_object Survival Object of form `survival::Surv`. Should be from the left side of the survival model function.
#' @param model Model object. Currently formatted to work with flexsurvreg objects, and may work for other types of models.
#' @param ibs Logical. If True, calculates the integrated Brier Score, which is a helpful fit statistic but is considerably slower to calculate than all other statistics.
#'
#' @returns List of fit statistics for the model.
#'
#' @examples
#' require(survival)
#' require(flexsurv)
#' flexsurvreg(Surv(time,status) ~age +sex, data=cancer, dist= 'weibull') -> model
#' get_fit_stats(Surv(cancer$time, cancer$status), model = model, ibs = FALSE)
#'
#' @export
get_fit_stats <- function(Surv_object, model, ibs=FALSE){
  # Rely on these three packages for statistics
  requireNamespace('survival')
  requireNamespace('survAUC')
  requireNamespace('SurvMetrics')

  # quick sub-function for deconstructing {predict.flexsurvreg(type='survival')}
  debulk_survprob <- function(list_of_tibbles){
    temp <- rbind(as.data.frame(list_of_tibbles$.pred[[1]])$.pred_survival)
    for(i in 2:nrow(list_of_tibbles)){
      as.data.frame(list_of_tibbles$.pred[[i]])$.pred_survival -> newtemp
      temp %>% rbind(newtemp) -> temp
    }
    rownames(temp) <- NULL
    return(as.matrix(temp))
  }

  # Four time-flavors: observed times, IQR, 10-90, 1:max (aka full)
  quants <- stats::quantile(Surv_object, c(.1, .9, .25, .75))$quantile %>% as.numeric()

  times       <- Surv_object[,1] %>% as.numeric() %>% unique() %>% sort()
  max_time    <- max(times)
  med_time    <- stats::median(times)
  avg_time    <- mean(times)

  # alternative time frames
  times_iqr   <- seq(quants[1], quants[2], 1)
  times_wide  <- seq(quants[3], quants[4], 1)
  times_full  <- seq(1, max_time, 1)

  # basic predictions for the base data
  preds      <- tryCatch(as.vector(unlist(stats::predict(model))), error=function(e){
    as.vector(unlist(stats::predict(model, type='rmst')$.pred_rmst)) # use RMST if there are errors in predictions
  })

  # survival rate predictions at {times}, used in IAEISE, IBS
  survprob   <- debulk_survprob(stats::predict(model, type='survival', times= times))
  survprob2  <- debulk_survprob(stats::predict(model, type='survival', times= times_iqr))
  survprob3  <- debulk_survprob(stats::predict(model, type='survival', times= times_wide))
  survprob4  <- debulk_survprob(stats::predict(model, type='survival', times= times_full))

  ### original survival package methods
  # harrel_c   <- survival::concordance(Surv_object ~ preds)$concordance
  # uno_c      <- survival::concordance(Surv_object ~ preds, timewt = 'n/G2')$concordance

  ### survAUC functions
  uno_c_auc  <- tryCatch(survAUC::UnoC(Surv_object, Surv_object, lpnew = -preds)    , error=function(e){NA})
  # iauc       <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = times)$iauc, error=function(e){NA})
  iauc2      <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = times_iqr)$iauc, error=function(e){NA})   # IQR version
  iauc3      <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = times_wide)$iauc, error=function(e){NA})  # Wider Version (10/90)
  iauc4      <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = times_full)$iauc, error=function(e){NA})  # Full version
  iauc5      <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = c(med_time))$iauc, error=function(e){NA}) # AUC at median
  iauc6      <- tryCatch(survAUC::AUC.uno(Surv_object, Surv_object, lpnew = -preds, times = c(avg_time))$iauc, error=function(e){NA}) # AUC at mean

  ### metrics
  c_index    <- tryCatch(SurvMetrics::Cindex(Surv_object, predicted = preds), error=function(e){NA}) %>% as.numeric()
  mae        <- tryCatch(SurvMetrics::MAE(Surv_object, pre_time= preds), error=function(e){NA}) %>% as.numeric()
  # i_stats    <- tryCatch(SurvMetrics::IAEISE(Surv_object, sp_matrix = survprob,  IRange=times), error=function(e){NA}) %>% as.numeric()  ### vastly underestimates I stats?
  i_stats2   <- tryCatch(SurvMetrics::IAEISE(Surv_object, sp_matrix = survprob2, IRange=times_iqr), error=function(e){NA})%>% as.numeric()
  i_stats3   <- tryCatch(SurvMetrics::IAEISE(Surv_object, sp_matrix = survprob3, IRange=times_wide), error=function(e){NA})%>% as.numeric()
  i_stats4   <- tryCatch(SurvMetrics::IAEISE(Surv_object, sp_matrix = survprob4, IRange=times_full), error=function(e){NA})%>% as.numeric()
  brier_med  <- tryCatch(SurvMetrics::Brier(Surv_object, pre_sp = as.vector(unlist(stats::predict(model, type='survival', times=med_time)$.pred_survival)), t_star = med_time), error=function(e){NA}) %>% as.numeric()
  brier_avg  <- tryCatch(SurvMetrics::Brier(Surv_object, pre_sp = as.vector(unlist(stats::predict(model, type='survival', times=avg_time)$.pred_survival)), t_star = avg_time), error=function(e){NA}) %>% as.numeric()

  # output
  out <- list(
    # 'Harrel.C.Index'      = harrel_c,
    # 'Uno.C.Index'         = uno_c,
    'C.Index.Uno'         = uno_c_auc,
    # 'iAUC'                = iauc,
    'iAUC.IQR'            = iauc2,
    'iAUC.Q10.Q90'        = iauc3,
    'iAUC.Full'           = iauc4,
    'AUC.Median'          = iauc5,
    'AUC.Mean'            = iauc6,
    'C.Index'             = c_index,
    'MAE'                 = mae,
    # 'IAE'                 = i_stats[1],
    # 'ISE'                 = i_stats[2],
    'IAE.IQR'             = i_stats2[1],
    'ISE.IQR'             = i_stats2[2],
    'IAE.Q10.Q90'         = i_stats3[1],
    'ISE.Q10.Q90'         = i_stats3[2],
    'IAE.Full'            = i_stats4[1],
    'ISE.Full'            = i_stats4[2],
    'Brier.Median'        = brier_med,
    'Brier.Mean'          = brier_avg
  )

  if(ibs){
    # ibs1       <- SurvMetrics::IBS(Surv_object, sp_matrix = survprob , IBSrange = times)%>% as.numeric()  ### Vastly underestimates
    # ibs2       <- SurvMetrics::IBS(Surv_object, sp_matrix = survprob2, IBSrange = times_iqr)%>% as.numeric()
    # ibs3       <- SurvMetrics::IBS(Surv_object, sp_matrix = survprob3, IBSrange = times_wide)%>% as.numeric()
    ibs4       <- SurvMetrics::IBS(Surv_object, sp_matrix = survprob4, IBSrange = times_full)%>% as.numeric()
    out <- c(
      out,
      # 'IBS'                 = ibs1,
      # 'IBS.IQR'             = ibs2,
      # 'IBS.Q10.Q90'         = ibs3,
      'IBS.Full'            = ibs4
    )
  }

  return(out)
}
