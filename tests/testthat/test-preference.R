context('preference function calls')

#############################
### Sample Size Functions ###
#############################

test_that("selection_sample_size function works", {
  resp <- selection_sample_size(power=0.8, phi=c(0.5, 0.5), sigma2=c(1,1), 
    delta_pi=1, delta_nu=0.5,xi=c(0.3,0.7),nstrata=2)
  expect_equal(resp, 566)
  expect_is(resp, 'numeric')
})


test_that("preference_sample_size function works", {
  resp <- preference_sample_size(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), 
    delta_pi=1, delta_nu=0.5,xi=c(0.3,0.7),nstrata=2)
  expect_equal(resp, 130)
  expect_is(resp, 'numeric')
})


test_that("treatment_sample_size function works", {
  resp <- treatment_sample_size(power=0.8, sigma2=c(1, 1), delta_tau=1.5, 
    xi=c(0.3,0.7),nstrata=2)
  expect_equal(resp, 28)
  expect_is(resp, 'numeric')
})

test_that("overall_sample_size function works", {
  resp <- overall_sample_size(power=0.8, phi=c(0.5,0.4), sigma2=c(1, 1), 
    delta_pi=1, delta_nu=0.5, delta_tau=1.5, xi=c(0.3,0.7),nstrata=2)
  expect_equal(resp$treatment, 28)
  expect_equal(resp$selection, 596)
  expect_equal(resp$preference, 138)
  expect_is(resp, 'list')
})

#######################
### Power Functions ###
#######################

test_that("selection_power function works", {
  resp <- selection_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, 
    delta_nu=0.5, xi=c(0.5,0.5), nstrata=2)
  expect_equal(round(resp,6), 0.508326)
  expect_is(resp, 'numeric')
})

test_that("preference_power function works", {
  resp <- preference_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, 
    delta_nu=0.5, xi=c(0.5,0.5), nstrata=2)
  expect_equal(round(resp,6), 0.984880)
  expect_is(resp, 'numeric')
})

test_that("treatment_power function works", {
  resp <- treatment_power(N=300, sigma2=c(1,1), delta_tau=0.5, xi=c(0.5,0.5), 
    nstrata=2)
  expect_equal(round(resp,6), 0.864747)
  expect_is(resp, 'numeric')
})

test_that("overall_power function works", {
  resp <- overall_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, 
                    delta_nu=0.5, delta_tau=0.5, xi=c(0.5,0.5), nstrata=2)
  expect_equal(round(resp$trt_pwr,6), 0.864747)
  expect_equal(round(resp$pref_pwr,6), 0.984880)
  expect_equal(round(resp$sel_pwr,6), 0.508326)
  expect_is(resp, 'list')
})

##########################
### Analysis Functions ###
##########################

test_that("analyze_raw_data function works", {
  x1 <- c(10,8,6,10,5)
  s11 <- c(1,1,2,2,2)
  x2 <- c(8,7,6,10,12,11,6,8)
  s22 <- c(1,1,1,1,2,2,2,2)
  y1 <- c(10,5,7,9,12,6)
  s1 <- c(1,1,1,2,2,2)
  y2 <- c(8,9,10,7,8,11)
  s2 <- c(1,1,1,2,2,2)
  resp <- analyze_raw_data(x1, x2, y1, y2, s11=s11, s22=s22, s1=s1, s2=s2, 
    xi=c(0.5,0.5), nstrata=2)
  expect_equal(round(resp$pref_test,6), -0.292573)
  expect_equal(round(resp$sel_test,6), 0.333207)
  expect_equal(round(resp$treat_test,6), -0.453945)
  expect_is(resp, 'data.frame')
})

test_that("analyze_summary_data function works", {
  x1mean <- 5
  x1var <- 1
  m1 <- 15
  x2mean <- 7
  x2var <- 1.1
  m2 <- 35
  y1mean <- 6
  y1var <- 1
  n1 <- 25
  y2mean <- 8
  y2var <- 1.2
  n2 <- 25
  resp <- analyze_summary_data(x1mean, x2var, m1, x2mean, x2var, m2, y1mean, 
    y2var, n1, y2mean, y2var, n2)
  expect_equal(round(resp$pref_test,6), -4.461299)
  expect_equal(round(resp$sel_test,6), 1.544837)
  expect_equal(round(resp$treat_test,6), -6.454972)
  expect_is(resp, 'data.frame')
})

#######################
### Other Functions ###
#######################

test_that("treatment_effect_size function works", {
  resp <- treatment_effect_size(N=300, power=0.9, sigma2=c(1,0.8), 
    xi=c(0.3,0.7), nstrata=2)
  expect_equal(round(resp,6), 0.490887)
  expect_is(resp, 'numeric')
})

test_that("optimal_proportion function works", {
  resp <- optimal_proportion(w_sel=0.2, w_pref=0.4, w_treat=0.4, sigma2=1, 
    phi=0.5, delta_pi=1, delta_nu=0.5)
  expect_equal(round(resp,6), 0.451026)
  expect_is(resp, 'numeric')
})

test_that("effects_from_means function works", {
  resp <- effects_from_means(mu1=1, mu2=2, mu11=1.5, mu22=2.5, phi=0.5)
  expect_equal(resp$treatment, -1)
  expect_equal(resp$selection, 0)
  expect_equal(resp$preference, 1)
  expect_is(resp, 'list')
})

