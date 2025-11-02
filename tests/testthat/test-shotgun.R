## raw testing chunk from legacy testing
if(F){

  ###
  # coxph(survival::Surv(time, status) ~ 1, data = survival::aml)
  # coxph(survival::Surv(time, status) ~ factor(x), data = survival::aml)

  ### Verify flexsurvreg and survreg work the same for us (same LL, AIC, etc)
  # survreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')
  # flexsurvreg(survival::Surv(time, status) ~ 1, data = survival::aml, dist = 'exp')

  # simple model tests
  fssg(survival::Surv(time, status) ~ 1, data = survival::aml,    dump_models = T, warn = F, spline=c('rp','wy'), detailed = T, ibs=T) -> test_models
  fssg(survival::Surv(time, status) ~ 1, data = survival::cancer, dump_models = T, warn = F) -> test_models2

  # specifying a model
  fssg(survival::Surv(time, status) ~ 1, data = survival::aml, detailed = T, ibs=T, model='gamgomp')

  fssg(survival::Surv(time, status) ~ 1, data = survival::aml, detailed = T, ibs=T, model=c('gamgomp','lindley'))

  # single variable models
  fssg(survival::Surv(time, status) ~ x,   data = survival::aml,    dump_models = T, warn = F) -> test_models3
  fssg(survival::Surv(time, status) ~ sex, data = survival::cancer, dump_models = T, warn = F) -> test_models4

  # more variables
  fssg(survival::Surv(time, status) ~ age + sex, data = survival::cancer, dump_models = T, warn = F, spline=NA, ibs=F, max_knots=5, skip='NA') -> test_models5
}
