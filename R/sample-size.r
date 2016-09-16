#############################################
# Functions with stratum-specific variances #
#############################################

# TODO: Fix the order of the parameters (ones with defaults go last.

#' Stratified Selection Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha Fill this in.
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_pi Fill this in.
#' @param delta_nu Fill this in.
#' @param xi Fill this in.
#' @param nstrata Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_selection<-function(zalpha=qnorm(0.975), zbeta, theta=0.5, phi, sigma2, delta_pi, delta_nu,  xi=c(0.5,0.5), nstrata=2) {
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))*(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_nu+delta_pi)^2+2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_nu^2)*sum_total
  return(N)
}

#' Stratified Preference Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha  Fill this in.
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_pi Fill this in.
#' @param delta_nu Fill this in.
#' @param xi Fill this in.
#' @param nstrata Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_preference<-function(zalpha=qnorm(0.975), zbeta, theta=0.5, phi, 
                           sigma2, delta_pi, delta_nu, xi=c(0.5,0.5), 
                           nstrata=2) {
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))*(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_pi+delta_nu)^2+2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_pi^2)*sum_total
  return(N)
}

#' Stratified Preference Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha  Fill this in.
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_tau Fill this in.
#' @param xi Fill this in.
#' @param nstrata Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
strat_treatment<-function(zalpha=qnorm(0.975), zbeta, theta=0.5, phi, 
                          sigma2, delta_tau, xi=c(0.5,0.5), nstrata=2) {
  terms=sapply(1:nstrata, function(x) xi[x]*sigma2[x])
  sum_total=sum(terms)
  N=4*(zalpha+zbeta)^2/((1-theta)*delta_tau^2)*sum_total
  return(N)
}

#' Stratified Preference Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha Fill this in. 
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_pi Fill this in.
#' @param delta_nu Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_selection<-function(zalpha=qnorm(0.975), zbeta, theta=0.5, phi, 
                            sigma2=1, delta_pi, delta_nu) {
  longterm=sigma2+phi*(1-phi)*((2*phi-1)*delta_nu+delta_pi)^2+2*(theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  N=(zalpha+zbeta)^2/(4*theta*delta_nu^2*phi^2*(1-phi)^2)*longterm
  return(N)
}

#' Unstratified Preference Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha Fill this in. 
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_pi Fill this in.
#' @param delta_nu Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_preference<-function(zalpha=qnorm(0.975), zbeta, theta=0.5, phi, 
                             sigma2=1, delta_pi, delta_nu) {
  longterm=sigma2+phi*(1-phi)*((2*phi-1)*delta_pi+delta_nu)^2+2*(theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  N=(zalpha+zbeta)^2/(4*theta*delta_pi^2*phi^2*(1-phi)^2)*longterm
  return(N)
}

#' Unstratified Treatment Effect
#'
#' Put a longer description of what the function does here.
#'
#' @param zalpha Fill this in. 
#' @param zbeta Fill this in.
#' @param theta Fill this in.
#' @param phi Fill this in.
#' @param sigma2 Fill this in.
#' @param delta_tau Fill this in.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
unstrat_treatment<-function(zalpha=qnorm(0.0975), zbeta, theta=0.5, sigma2=1, 
                            delta_tau){
  N=(4*sigma2*(zbeta+zalpha)^2)/((1-theta)*delta_tau^2)
  return(N)
}

#######################################################
# Function to compute variance of mixture normal dist #
#######################################################

# Find sigma^2 for unstratified case
sigma2_mixtures<-function(sigma2,means,prop){
  # Overall mixture
  val=prop[1]*sigma2[1]+prop[2]*sigma2[2]+(prop[1]*means[1]^2+prop[2]*means[2]^2-(prop[1]*means[1]+prop[2]*means[2])^2)
  return(val)
}

# Find sigma^2 for each stratum based on means
sigma2_stratum<-function(sigma,mu,tau,nu,pi,prop,theta=c(0.5,0.5)){
  mean_choice=sapply(1:2, function(x) mu[x]+tau[x]+nu[x]+pi[x])
  mean_random=sapply(1:2, function(x) mu[x]+tau[x])
  # Compute mean, variance for each stratum
  means=sapply(1:2, function(x) theta[x]*mean_choice[x]+(1-theta[x])*mean_random[x])
  sigma2_derived=sapply(1:2, function(x) theta[x]*sigma[x]+(1-theta[x])*sigma[x]+(theta[x]*mean_choice[x]^2+(1-theta[x])*mean_random[x]^2-(theta[x]*mean_choice[x]+(1-theta[x])*mean_random[x])^2))
  return(sigma2_derived)
}

# Find means from effects
means_stratum<-function(sigma,mu,tau,nu,pi,prop,theta=c(0.5,0.5)){
  mean_choice=sapply(1:2, function(x) mu[x]+tau[x]+nu[x]+pi[x])
  mean_random=sapply(1:2, function(x) mu[x]+tau[x])
  # Compute mean value for each stratum
  means=sapply(1:2, function(x) theta[x]*mean_choice[x]+(1-theta[x])*mean_random[x])
  return(means)
}

