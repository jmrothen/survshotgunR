test_that("extra distribution functions work",{
  ### Survival Functions
  # Survivify with one of our functions
  expect_equal(
    round(pgamgomp(1,1,1,1, lower.tail = F),7),
    {
      sgamgomp <- survivify(pgamgomp)
      round(sgamgomp(1,1,1,1),7)
    }
  )

  # Survivify with a base p-function
  expect_equal(
    round(pnorm(1,1,1, lower.tail = F),7),
    {
      snorm <- survivify(pnorm)
      round(snorm(1,1,1),7)
    }
  )

  # Survivify with imported p function
  expect_equal(
    round(VGAM::pfrechet(2,1,1,1, lower.tail = F),7),
    {
      sfrechet <- survivify(VGAM::pfrechet)
      round(sfrechet(2,1,1,1),7)
    }
  )

  ### Hazard Functions
  # Hazardify with one of our functions
  expect_equal(
    round(dgamgomp(1,1,1,1) /  (1-pgamgomp(1,1,1,1)),7),
    {
      hgamgomp <- hazardify(dgamgomp, pgamgomp)
      round(hgamgomp(1,1,1,1),7)
    }
  )

  # Hazardify with a base p-function
  expect_equal(
    round(dnorm(1,1,1) / (1-pnorm(1,1,1)),7),
    {
      hnorm <- hazardify(dnorm, pnorm)
      round(hnorm(1,1,1),7)
    }
  )

  # Hazardify with imported p function
  expect_equal(
    round(VGAM::dfrechet(2,1,1,1) / (1-VGAM::pfrechet(2,1,1,1)),7),
    {
      hfrechet <- hazardify(VGAM:::dfrechet, VGAM::pfrechet)
      round(hfrechet(2,1,1,1),7)
    }
  )

  ### Cumulative Hazard functions
  # Cum-Hazardify with one of our functions
  expect_equal(
    round(-pgamgomp(1,1,1,1,F,T),7),
    {
      Hgamgomp <- cumhazardify(pgamgomp)
      round(Hgamgomp(1,1,1,1),7)
    }
  )

  # Cum-Hazardify with a base p-function
  expect_equal(
    round(-pnorm(1,1,1,F,T), 7),
    {
      Hnorm <- cumhazardify(pnorm)
      round(Hnorm(1,1,1),7)
    }
  )

  # Cum-Hazardify with imported p function
  expect_equal(
    round(-VGAM::pfrechet(2,1,1,1,F,T),7),
    {
      Hfrechet <- cumhazardify(VGAM::pfrechet)
      round(Hfrechet(2,1,1,1),7)
    }
  )

})


# simple checks about the distribution function
test_that("Distribution list runs",{
  expect_type(
    fssg_dist_list(),
    'list'
  )
  expect_gt(
    length(fssg_dist_list()),
    50
  )
})


# tests for accessing attributes of the elements within fssg_dist_list
test_that("Distribution list items works", {
  expect_equal(
    fssg_dist_list()[[1]]$name, 'genf'
  )
  expect_equal(
    fssg_dist_list()$genf$location, 'mu'
  )
  expect_equal(
    round(fssg_dist_list()$lindley$inits(c(1,2,3,4,5)),7), round(0.3333333,7)
  )
  expect_equal(
    fssg_dist_list()$lindley$d(1,2), VGAM::dlind(1,2)
  )
  expect_equal(
    fssg_dist_list()$lindley$h(1,2) ,
    {
      hazardify(VGAM::dlind, VGAM::plind) -> hlind
      hlind(1,2)
    }
  )
  expect_in(
    c('lindley', 'genf','weibull','gamgomp'),
    names(fssg_dist_list())
  )
})
