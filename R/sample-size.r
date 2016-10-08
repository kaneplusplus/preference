#' Stratified Selection Effect Sample Size
#'
#' Calculates the sample size required to detect a given selection effect 
#' in a stratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_selection<-function(beta, phi, sigma2, delta_pi, delta_nu, 
                          alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                          nstrata=2) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
              *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_nu+delta_pi)^2
              +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_nu^2)*sum_total
  return(N)
}

#' Stratified Preference Effect Sample Size
#'
#' Calculates the sample size required to detect a given preference effect 
#' in a stratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_preference(zbeta=1.282, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_pi=1, 
#'  delta_nu=0)
#' @export
strat_preference<-function(beta, phi, sigma2, delta_pi, delta_nu, 
                           alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                           nstrata=2) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
              *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_pi+delta_nu)^2
              +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_pi^2)*sum_total
  return(N)
}

#' Stratified Treatment Effect Sample Size
#'
#' Calculates the sample size required to detect a given treatment effect 
#' in a stratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_treatment<-function(beta, phi, sigma2, delta_tau, alpha=0.05,
                          theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) xi[x]*sigma2[x])
  sum_total=sum(terms)
  N=4*(zalpha+zbeta)^2/((1-theta)*delta_tau^2)*sum_total
  return(N)
}

#' Stratified Overall Sample Size
#'
#' Calculates the sample size required to detect a given set of effects 
#' in a stratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_overall<-function(beta, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                           alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                           nstrata=2) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  pref=strat_preference(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, 
                        theta, xi, nstrata)
  sel=strat_selection(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, 
                      theta, xi, nstrata)
  treat=strat_treatment(zbeta, phi, sigma2, delta_tau, zalpha, 
                        theta, xi, nstrata)
  return(max(pref,sel,treat))
}

#' Unstratified Selection Effect Sample Size
#'
#' Calculates the sample size required to detect a given selection effect 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi proportion of patients preferring treatment 1.
#' @param sigma2 variance estimate.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization (default=0.5).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_selection<-function(beta, phi, sigma2, delta_pi, delta_nu, 
                            alpha=0.05, theta=0.5) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  longterm=sigma2+phi*(1-phi)*((2*phi-1)*delta_nu+delta_pi)^2
          +2*(theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  N=(zalpha+zbeta)^2/(4*theta*delta_nu^2*phi^2*(1-phi)^2)*longterm
  return(N)
}

#' Unstratified Preference Effect Sample Size
#'
#' Calculates the sample size required to detect a given preference effect 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param beta desired study power
#' @param phi proportion of patients preferring treatment 1.
#' @param sigma2 variance estimate.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization (default=0.5).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_preference<-function(beta, phi, sigma2, delta_pi, delta_nu, 
                             alpha=0.05, theta=0.5) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  longterm=sigma2+phi*(1-phi)*((2*phi-1)*delta_pi+delta_nu)^2
          +2*(theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  N=(zalpha+zbeta)^2/(4*theta*delta_pi^2*phi^2*(1-phi)^2)*longterm
  return(N)
}

#' Unstratified Treatment Effect Sample Size
#'
#' Calculates the sample size required to detect a given treatment effect 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi proportion of patients preferring treatment 1.
#' @param sigma2 variance estimate.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization (default=0.5).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_treatment<-function(beta, sigma2=1, delta_tau, alpha=0.05, 
                            theta=0.5){
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  N=(4*sigma2*(zbeta+zalpha)^2)/((1-theta)*delta_tau^2)
  return(N)
}

#' Unstratified Overall Sample Size
#'
#' Calculates the sample size required to detect a set of effects 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param beta desired study power.
#' @param phi proportion of patients preferring treatment 1.
#' @param sigma2 variance estimate.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization (default=0.5).
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_overall<-function(beta, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                            alpha=0.05, theta=0.5) {
  zbeta<-qnorm(beta)
  zalpha<-qnorm(1-(alpha/2))
  sel=unstrat_selection(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, theta)
  pref=unstrat_preference(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, 
                          theta)
  treat=unstrat_treatment(zbeta, phi, sigma2, delta_tau, zalpha, theta)
  return(max(sel, pref, treat))
}

#' Unstratified Optimized Theta
#'
#' Calculates the optimal proportion of patients assigned to the choice arm
#' in a two-stage randomized trial
#'
#' @param w_sel weight assigned to the estimation of the selection effect. Sum 
#'              of thre weights should be 1.
#' @param w_pref weight assigned to the estimation of the preference effect. Sum 
#'               of thre weights should be 1.
#' @param w_treat weight assigned to estimation of the treatment effect. Sum 
#'                of thre weights should be 1.
#' @param sigma2 variance estimate.
#' @param phi proportion of patients preferring treatment 1.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect. 
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
theta_optim<-function(w_sel,w_pref,w_treat,sigma2,phi,delta_pi,delta_nu) {
  if (w_sel+w_pref+w_treat!=1) stop('weights do not sum to 1')
  # Based on Equation 16 in Walter paper
  num<-w_sel+w_pref+phi*(1-phi)*((w_sel*((2*phi-1)*delta_nu+delta_pi)^2
                                  +w_pref*((2*phi-1)*delta_pi+delta_nu)^2)/sigma2)
  denom<-16*w_treat*phi^2*(1-phi)^2+2*(w_sel+w_pref)*(phi^2+(1-phi)^2)
  
  return(uniroot(f,c(0,1),value=(num/denom))$root)
}

# Function used in theta optimization function
f<-function(theta,value) {
  (theta/(1-theta))^2-value
}  


###################
# Extra functions #
###################

# Find sigma^2 for unstratified case
sigma2_mixtures<-function(sigma2,means,prop){
  # Overall mixture
  val=prop[1]*sigma2[1]+prop[2]*sigma2[2]+(prop[1]*means[1]^2+
      prop[2]*means[2]^2-(prop[1]*means[1]+prop[2]*means[2])^2)
  return(val)
}

# Find sigma^2 for each stratum based on means
sigma2_stratum<-function(sigma,mu,tau,nu,pi,prop,theta=c(0.5,0.5)){
  mean_choice=sapply(1:2, function(x) mu[x]+tau[x]+nu[x]+pi[x])
  mean_random=sapply(1:2, function(x) mu[x]+tau[x])
  # Compute mean, variance for each stratum
  means=sapply(1:2, function(x) theta[x]*mean_choice[x]+
                 (1-theta[x])*mean_random[x])
  sigma2_derived=sapply(1:2, function(x) theta[x]*sigma[x]+(1-theta[x])
                        *sigma[x]+(theta[x]*mean_choice[x]^2+(1-theta[x])
                        *mean_random[x]^2-(theta[x]*mean_choice[x]
                        +(1-theta[x])*mean_random[x])^2))
  return(sigma2_derived)
}

# Find means from effects
means_stratum<-function(sigma,mu,tau,nu,pi,prop,theta=c(0.5,0.5)){
  mean_choice=sapply(1:2, function(x) mu[x]+tau[x]+nu[x]+pi[x])
  mean_random=sapply(1:2, function(x) mu[x]+tau[x])
  # Compute mean value for each stratum
  means=sapply(1:2, function(x) theta[x]*mean_choice[x]+
              (1-theta[x])*mean_random[x])
  return(means)
}


