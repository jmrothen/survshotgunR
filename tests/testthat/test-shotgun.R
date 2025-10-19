## raw testing chunk from legacy testing
if(F){

  ###
  # coxph(survival::Surv(time, status) ~ 1, data = survival::aml)
  # coxph(survival::Surv(time, status) ~ factor(x), data = survival::aml)

  ### Verify flexsurvreg and survreg work the same for us (same LL, AIC, etc)
  # survreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')
  # flexsurvreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')

  # simple model tests
  surv_shotgun(survival::Surv(time, status) ~ 1, data = survival::aml,    dump_models = T, warn = F, spline=c('rp','wy'), detailed = T, ibs=T) -> test_models
  surv_shotgun(survival::Surv(time, status) ~ 1, data = survival::cancer, dump_models = T, warn = F) -> test_models2

  # single variable models
  surv_shotgun(survival::Surv(time, status) ~ x,   data = survival::aml,    dump_models = T, warn = F) -> test_models3
  surv_shotgun(survival::Surv(time, status) ~ sex, data = survival::cancer, dump_models = T, warn = F) -> test_models4

  # more variables
  surv_shotgun(survival::Surv(time, status) ~ age + sex, data = survival::cancer, dump_models = T, warn = F, spline=c('rp','wy'), ibs=F, max_knots=5) -> test_models5
}
