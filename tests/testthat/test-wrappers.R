test_that("d functions work", {
  ### Erlang
  # Wrapped R function
  expect_equal(
    round(derlang(1,1,1),7),
    round(0.3678794,7)
  )
  # Direct cpp version
  expect_equal(
    round(derlang_c(1,1,1),7),
    round(0.3678794,7)
  )
  # Vectorized function
  expect_equal(
    round(derlang(c(1,3,4,5),c(1,2),1),7),
    round(c(0.36787944, 0.14936121, 0.01831564, 0.03368973),7)
  )


  ### Gamma Gompertz
  # Wrapped R function
  expect_equal(
    round(dgamgomp(1,1,1,1),7),
    round(0.3678794,7)
  )
  # Direct cpp version
  expect_equal(
    round(dgamgomp_c(1,1,1,1),7),
    round(0.3678794,7)
  )
  # Vectorized
  expect_equal(
    round(dgamgomp(c(1,2,3,4),1,c(2,3),1),7),
    round(c(2.706706e-01, 7.436257e-03, 4.957504e-03, 1.843264e-05),7)
  )


  ### Log Cauchy
  # Wrapped
  expect_equal(
    round(dlogcauchy(1,1,1),7),
    round(0.1591549,7)
  )
  # CPP
  expect_equal(
    round(dlogcauchy_c(1,1,1),7),
    round(0.1591549,7)
  )
  # Vectorized
  expect_equal(
    round(dlogcauchy(c(1,2,3,4),c(1,2),1),7),
    round(c(0.1591549, 0.0587751, 0.1050814, 0.0578058),7)
  )


  ### Hypertabastic
  # Wrapped
  expect_equal(
    round(dhypertab(1,1,1),7),
    round(0.1701686,7)
  )
  # CPP
  expect_equal(
    round(dhypertab_c(1,1,1),7),
    round(0.1701686,7)
  )
  # Vectorized
  expect_equal(
    round(dhypertab(c(1,2,3,4),c(2,3),1),7),
    round(c(0.5440161, 0.2104008, 0.0692722, 0.0007311),7)
  )


  ### Inverse Lindley
  # Wrapped
  expect_equal(
    round(dinvlind(1,1),7),
    round(0.3678794,7)
  )
  # CPP
  expect_equal(
    round(dinvlind_c(1,1),7),
    round(0.3678794,7)
  )
  # Vectorized
  expect_equal(
    round(dinvlind(c(1,2,3,4),c(1,2)),7),
    round(c(0.3678794, 0.1839397, 0.0530764, 0.0631803),7)
  )
})




test_that("p functions work", {
  ### Erlang
  # Wrapped R function
  expect_equal(
    round(perlang(1,1,1),7),
    round(0.6321206,7)
  )
  # Direct cpp version
  expect_equal(
    round(perlang_c(1,1,1),7),
    round(0.6321206,7)
  )
  # Vectorized function
  expect_equal(
    round(perlang(c(1,3,4,5),c(1,2),1),7),
    round(c(0.6321206, 0.8008517, 0.9816844, 0.9595723),7)
  )


  ### Gamma Gompertz
  # Wrapped R function
  expect_equal(
    round(pgamgomp(1,1,1,1),7),
    round(0.6321206,7)
  )
  # Direct cpp version
  expect_equal(
    round(pgamgomp_c(1,1,1,1),7),
    round(0.6321206,7)
  )
  # Vectorized
  expect_equal(
    round(pgamgomp(c(1,2,3,4),1,c(2,3),1),7),
    round(c(0.8646647, 0.9975212, 0.9975212, 0.9999939),7)
  )


  ### Log Cauchy
  # Wrapped
  expect_equal(
    round(plogcauchy(1,2,1),7),
    round(0.1475836,7)
  )
  # CPP
  expect_equal(
    round(plogcauchy_c(1,2,1),7),
    round(0.1475836,7)
  )
  # Vectorized
  expect_equal(
    round(plogcauchy(c(1,2,3,4),c(1,2),1),7),
    round(c(0.2500000, 0.2079062, 0.5312881, 0.3247907),7)
  )


  ### Hypertabastic
  # Wrapped
  expect_equal(
    round(phypertab(1,1,1),7),
    round(0.0470717,7)
  )
  # CPP
  expect_equal(
    round(phypertab_c(1,1,1),7),
    round(0.0470717,7)
  )
  # Vectorized
  expect_equal(
    round(phypertab(c(1,2,3,4),c(2,3),1),7),
    round(c(0.1683770, 0.9205259, 0.9644561, 0.9997552),7)
  )


  ### Inverse Lindley
  # Wrapped
  expect_equal(
    round(pinvlind(1,1),7),
    round(0.5518192,7)
  )
  # CPP
  expect_equal(
    round(pinvlind_c(1,1),7),
    round(0.5518192,7)
  )
  # Vectorized
  expect_equal(
    round(pinvlind(c(1,2,3,4),c(1,2)),7),
    round(c(0.5518192, 0.4905059, 0.8359532, 0.7076191),7)
  )
})




test_that("Log and Tail works", {
  # ds
  expect_equal(
    log(derlang(1,1,1)),
    derlang(1,1,1,T)
  )
  expect_equal(
    log(dgamgomp(1,1,1,1)),
    dgamgomp(1,1,1,1,T)
  )
  expect_equal(
    log(dlogcauchy(1,1,1)),
    dlogcauchy(1,1,1,T)
  )
  expect_equal(
    log(dhypertab(1,1,1)),
    dhypertab(1,1,1,T)
  )
  expect_equal(
    log(dinvlind(1,1)),
    dinvlind(1,1,T)
  )


  # ps
  expect_equal(
    log(perlang(1,1,1)),
    perlang(1,1,1,log.p=T)
  )
  expect_equal(
    log(pgamgomp(1,1,1,1)),
    pgamgomp(1,1,1,1,log.p=T)
  )
  expect_equal(
    log(plogcauchy(1,1,1)),
    plogcauchy(1,1,1,log.p=T)
  )
  expect_equal(
    log(phypertab(1,1,1)),
    phypertab(1,1,1,log.p=T)
  )
  expect_equal(
    log(pinvlind(1,1)),
    pinvlind(1,1,log.p=T)
  )


  # tails
  expect_equal(
    1-perlang(1,1,1),
    perlang(1,1,1,F)
  )
  expect_equal(
    1-pgamgomp(1,1,1,1),
    pgamgomp(1,1,1,1,F)
  )
  expect_equal(
    1-plogcauchy(1,1,1),
    plogcauchy(1,1,1,F)
  )
  expect_equal(
    1-phypertab(1,1,1),
    phypertab(1,1,1,F)
  )
  expect_equal(
    1-pinvlind(1,1),
    pinvlind(1,1,F)
  )


  # tail and log
  expect_equal(
    log(1-perlang(1,1,1)),
    perlang(1,1,1,F,T)
  )
  expect_equal(
    log(1-pgamgomp(1,1,1,1)),
    pgamgomp(1,1,1,1,F,T)
  )
  expect_equal(
    log(1-plogcauchy(1,1,1)),
    plogcauchy(1,1,1,F,T)
  )
  expect_equal(
    log(1-phypertab(1,1,1)),
    phypertab(1,1,1,F,T)
  )
  expect_equal(
    log(1-pinvlind(1,1)),
    pinvlind(1,1,F,T)
  )
})
