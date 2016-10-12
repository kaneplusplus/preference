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
  # Error messages
  if(!is.numeric(beta) || beta <= 0 || beta >= 1)
    stop('Power must be numeric in [0,1]')
  if (length(phi) != nstrata) 
    stop('Length vector does not match number of strata')
  if(any(!is.numeric(phi)) || any(phi <= 0) || any(phi >= 1))
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2) != nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(!is.numeric(sigma2) || any(sigma2 <= 0)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) || !is.numeric(delta_nu))
    stop('Effect size must be numeric value')
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop('Type I error rate must be numeric in [0,1]')
  if(!is.numeric(theta) || theta <= 0 || theta >= 1)
    stop('Theta must be numeric in [0,1]')
  if(any(!is.numeric(xi) || any(xi <= 0) || any(xi >= 1)))
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi) != nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi) != 1) 
    stop('Stratum proportions do not sum to 1')
  if(!is.numeric(nstrata) || nstrata <=0)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata))
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata))
    stop('Number of strata must be numeric greater than 0')
  
  #Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) | !is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata))
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate sample size
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
  # Error messages
  if(beta<0 | beta>1 | !is.numeric(beta)) 
    stop('Power must be numeric in [0,1]')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) | !is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate sample size
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
  if(w_sel<0 | w_sel>1 | w_pref<0 | w_pref>1 | w_treat<0 | w_treat>1 | 
     any(!is.numeric(c(w_sel,w_pref,w_treat))))
    stop('Weights must be numeric value in [0,1]')
  if (w_sel+w_pref+w_treat!=1) 
    stop('weights do not sum to 1')
  if(sigma2<=0 | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(phi<0 | phi>1 | !is.numeric(phi)) 
    stop('Preference rate must be numeric value in [0,1]')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu))
    stop('Effect size must be numeric value')
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

#' Calculate Effect Sizes from Means
#'
#' Calculates the preference, selection and treatment effects given the means
#' of each treatment group in the choice and random arms for the 2-stage 
#' randomized study.
#'
#' @param mu1 mean response of the patients receiving treatment 1 in the 
#'            random arm. For unstratified design, should be numeric value.
#'            For the stratified design, should be vector of length equal to
#'            number of strata with each entry corresponding to stratum-
#'            specific mean.
#' @param mu2 mean response of the patients receiving treatment 2 in the 
#'            random arm. For unstratified design, should be numeric value.
#'            For the stratified design, should be vector of length equal to
#'            number of strata with each entry corresponding to stratum-
#'            specific mean.
#' @param mu11 mean response of the patients choosing treatment 1 in the choice
#'             arm. For unstratified design, should be numeric value.
#'            For the stratified design, should be vector of length equal to
#'            number of strata with each entry corresponding to stratum-
#'            specific mean.
#' @param mu22 mean response of the patients choosing treatment 2 in the choice
#'             arm. For unstratified design, should be numeric value.
#'            For the stratified design, should be vector of length equal to
#'            number of strata with each entry corresponding to stratum-
#'            specific mean.
#' @param phi proportion of patients preferring treatment 1. For unstratified 
#'            design, should be numeric value. For the stratified design, 
#'            should be vector of length equal to number of strata with each 
#'            entry corresponding to stratum-specific preference rate.
#' @param nstrata number of strata. Default is 1 (unstratified design).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. Should only be specified for stratified
#'          design.
#' @examples
#' # Put example code here.
#' rnorm(10)
#' @export
calc_effects<-function(mu1,mu2,mu11,mu22,phi,nstrata=1,xi=NULL) {
  # Error messages
  if(nstrata<=0 | !is.numeric(nstrata))
    stop('Number of strata must be numeric greater than 0')
  if (nstrata>1 & is.null(xi))
    stop('Must define xi for stratified design')
  if (length(phi)!=nstrata | length(mu1)!=nstrata | length(mu2)!=nstrata
      | length(mu11)!=nstrata |length(mu22)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(!is.numeric(mu1) | !is.numeric(mu2) | !is.numeric(mu11) | !is.numeric(mu22))
    stop('Mean must be numeric value')
  if((any(xi<0) | any(xi>1) | any(!is.numeric(xi))) & !is.null(xi))
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata & nstrata!=1) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1 & !is.null(xi)) 
    stop('Stratum proportions do not sum to 1')
  
  
  # Calculate unobserved means
  mu12<-(mu1-phi*mu11)/(1-phi)
  mu21<-(mu2-(1-phi)*mu22)/phi
  
  # Calculate effect sizes
  delta_tau<-mu1-mu2
  delta_nu<-(mu11+mu21-mu12-mu22)/2
  delta_pi<-(mu11-mu21-mu12+mu22)/2
  
  if (nstrata==1) { 
    # Unstratified case
    effects<-list("delta_tau"=delta_tau,"delta_nu"=delta_nu,
                  "delta_pi"=delta_pi)
  } else {
    # Stratified case
    effects<-list("delta_tau"=sum(sapply(1:nstrata, function(x) 
                              phi[x]*delta_tau[x])),
                  "delta_nu"=sum(sapply(1:nstrata, function(x) 
                              phi[x]*delta_nu[x])),
                  "delta_pi"=sum(sapply(1:nstrata,function(x) 
                              phi[x]*delta_pi[x])))
  }

  return(effects)
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


