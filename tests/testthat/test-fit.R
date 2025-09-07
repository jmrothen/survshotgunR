test_that("Fit stats seem to work",{
  expect_equal(
    {
      shotgun_dist_list()$lindley -> testdist
      flexsurv::flexsurvreg(survival::Surv(time,status) ~age +sex, data=survival::cancer, dist= testdist, dfns = list(d=testdist$d, p=testdist$p)) -> testmodel
      get_fit_stats(survival::Surv(survival::cancer$time, survival::cancer$status), model = testmodel, ibs = F) -> statlist
      for(i in 1:length(statlist)){
        statlist[[i]] <- round(statlist[[i]],7)
      }
      statlist
    },{
      list(
        C.Index.Uno= round(0.5848048,7),
        iAUC.IQR = round(0.6407606,7),
        iAUC.Q10.Q90 =  round(0.6174132,7),
        iAUC.Full = round(0.6446272,7),
        AUC.Median = round(0.6025951,7),
        AUC.Mean = round(0.6303409,7),
        C.Index = round(0.602809,7),
        MAE = round(193.4215,7),
        IAE.IQR = round(19.1671,7),
        ISE.IQR = round(0.7022,7),
        IAE.Q10.Q90 = round(54.9943,7),
        ISE.Q10.Q90 = round(6.9539,7),
        IAE.Full = round(14.6904,7),
        ISE.Full = round(0.3388,7),
        Brier.Median = round(0.393042,7),
        Brier.Mean = round(0.237512,7)
        #,IBS.Full=round(0.41309,7)
      )
    }
  )
})



# if(F){
#   ### MODEL PREP
#   shotgun_dist_list()$lindley -> testdist
#   shotgun_dist_list()$fatigue -> testdist2
#   shotgun_dist_list()$gompertz -> testdist3
#
#   flexsurvreg(survival::Surv(time,status) ~age +sex, data=survival::cancer, dist= testdist, dfns = list(d=testdist$d, p=testdist$p)) -> pp4
#   flexsurvreg(survival::Surv(time,status) ~age +sex, data=survival::cancer, dist= testdist2, dfns = list(d=testdist2$d, p=testdist2$p)) -> pp5
#   flexsurvreg(survival::Surv(time,status) ~age +sex, data=survival::cancer, dist= testdist3, dfns = list(d=testdist3$d, p=testdist3$p)) -> pp6
#   flexsurvreg(survival::Surv(time,status) ~1, data=survival::cancer, dist= testdist2, dfns = list(d=testdist2$d, p=testdist2$p)) -> pp_ref
#
#
#
#   # examples
#
#   get_fit_stats(Surv(cancer$time, cancer$status), model = pp4, ibs = F) -> cool
#   get_fit_stats(Surv(cancer$time, cancer$status), model = pp4, ibs = T) # much slower, but works
#   get_fit_stats(Surv(cancer$time, cancer$status), model = pp5, ibs = F) -> cool2
#
#   get_fit_stats(Surv(cancer$time, cancer$status), model = pp6, ibs = F) -> cool3
#
#   get_fit_stats(Surv(cancer$time, cancer$status), model = pp_ref, ibs = F) -> cool3
#
# }
