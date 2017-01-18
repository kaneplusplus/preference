##############################
#### STRATIFIED FUNCTIONS ####
##############################

#' Stratified Selection Effect Sample Size
#'
#' Calculates the sample size required to detect a given selection effect 
#' in a stratified two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_selection(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_pi=1, 
#'  delta_nu=0.5)
#' @export
strat_selection<-function(power, phi, sigma2, delta_pi, delta_nu, 
                          alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                          nstrata=2) {
  # Error messages
  if(!is.numeric(power) || power <= 0 || power >= 1)
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
  zbeta<-qnorm(power)
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
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_preference(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_pi=1, 
#'  delta_nu=0.5)
#' @export
strat_preference<-function(power, phi, sigma2, delta_pi, delta_nu, 
                           alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                           nstrata=2) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
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
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_treatment(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_tau=1)
#' @export
strat_treatment<-function(power, phi, sigma2, delta_tau, alpha=0.05,
                          theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) xi[x]*sigma2[x])
  sum_total=sum(terms)
  N=4*(zalpha+zbeta)^2/((1-theta)*delta_tau^2)*sum_total
  return(N)
}

#' Stratified Overall Sample Size
#'
#' Calculates the sample size required to detect a given set of effects 
#' in a stratified two-stage randomized clinical trial. Returns the largest
#' of the required sample sizes for a given set of treatment, selection, and 
#' preference effects.
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_overall(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_pi=1, 
#'  delta_nu=0.5, delta_tau=0.5)
#' @export
strat_overall<-function(power, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                           alpha=0.05, theta=0.5, xi=c(0.5,0.5), 
                           nstrata=2) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  pref=strat_preference(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, 
                        theta, xi, nstrata)
  sel=strat_selection(zbeta, phi, sigma2, delta_pi, delta_nu, zalpha, 
                      theta, xi, nstrata)
  treat=strat_treatment(zbeta, phi, sigma2, delta_tau, zalpha, 
                        theta, xi, nstrata)
  return(max(pref,sel,treat))
}

#' Stratified Treatment Effect Power Calculation
#'
#' Calculates the study power to detect the treatment effect given a particular 
#' sample size in an stratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_trt_pwr(N=300, sigma2=c(1,1), delta_tau=0.5)
#' @export
strat_trt_pwr<-function(N, sigma2, delta_tau, alpha=0.05, 
                        theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  power=pnorm(sqrt((((1-theta)*delta_tau^2*N)/4)*
              sum(sapply(1:nstrata, function(i) xi[i]*sigma2[i])))-zalpha)
  
  return(power)
}

#' Stratified Preference Effect Power Calculation
#'
#' Calculates the study power to detect the preference effect given a particular 
#' sample size in an stratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_pref_pwr(N=300, phi=c(0.5,0.6), sigma2=c(1,1), delta_pi=0.5, 
#'                delta_nu=0.5)
#' @export
strat_pref_pwr<-function(N, phi, sigma2, delta_pi, delta_nu, alpha=0.05, 
                         theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
               *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_pi+delta_nu)^2
                 +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  power=pnorm(sqrt((4*theta*delta_pi^2*N)/(sum_total))-zalpha)
  
  return(power)
}

#' Stratified Selection Effect Power Calculation
#'
#' Calculates the study power to detect the selection effect given a particular 
#' sample size in an stratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_sel_pwr(N=300, phi=c(0.5,0.6), sigma2=c(1,1), delta_pi=0.5, 
#'                delta_nu=0.5)
#' @export
strat_sel_pwr<-function(N, phi, sigma2, delta_pi, delta_nu, alpha=0.05, 
                         theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
               *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_nu+delta_pi)^2
                 +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  power=pnorm(sqrt((4*theta*delta_nu^2*N)/(sum_total))-zalpha)
  
  return(power)
}

#' Stratified Power Calculation from Sample Size
#'
#' Calculates the study power to detect a set of effects given a particular 
#' sample size in an stratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi vector of the proportion of patients preferring treatment 1 within
#'            each stratum. Length of vector should equal number of strata 
#'            in the study. Vector elements should be numeric values between 0
#'            and 1.
#' @param sigma2 vector of within-stratum variances. Length of vector should 
#'               equal number of strata in the study. Vector elements should be
#'               positive numeric values.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_power(N=300, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5,
#' delta_tau=0.5)
#' @export
strat_power<-function(N, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                        alpha=0.05, theta=0.5, xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  
  # Calculate study power
  trt_pwr<-strat_trt_pwr(N=N,sigma2=sigma2,delta_tau=delta_tau,alpha=alpha,
                         theta=theta,xi=xi,nstrata=nstrata)
  pref_pwr<-strat_pref_pwr(N=N,phi=phi,sigma2=sigma2,delta_pi=delta_pi,
                           delta_nu=delta_nu,alpha=alpha,theta=theta,xi=xi,
                           nstrata=nstrata)
  sel_pwr<-strat_sel_pwr(N=N,phi=phi,sigma2=sigma2,delta_pi=delta_pi,
                         delta_nu=delta_nu,alpha=alpha,theta=theta,xi=xi,
                         nstrata=nstrata)
  
  return(data.frame(trt_pwr=trt_pwr,pref_pwr=pref_pwr,sel_pwr=sel_pwr))  
}

#' Stratified Treatment Effect Back Calculation
#' 
#' Calculates the treatment effect that can be detected given a desired study 
#' power and overall study sample size for the stratified two-stage 
#' randomized design
#' 
#' @param N overall study sample size.
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2).
#' @examples
#' strat_trt_effect(N=300,power=0.9,sigma2=c(1,0.8))
#' @export
strat_trt_effect<-function(N, power, sigma2, alpha=0.05, theta=0.5, 
                           xi=c(0.5,0.5), nstrata=2) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
  if(power<0 | power>1 | !is.numeric(power)) 
    stop('Power must be numeric in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
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
  
  # Calculate effect size
  zbeta=qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  effect=sqrt(((4*(zbeta+zalpha)^2)/((1-theta)*N))*
                sum(sapply(1:nstrata, function(i) xi[i]*sigma2[i])))
  
  return(effect)
}


#' Stratified Analysis Function
#' 
#' Computes the test statistic and p-value for the preference, selection, and 
#' treatment effects for the stratified two-stage randomized trial
#' 
#' @param x1 vector of responses for patients choosing treatment 1
#' @param s11 vector of stratum membership for patients choosing treatment 1. 
#'            Should be a vector of the same length as x1 with the number of
#'            unique values equal to the number of strata.
#' @param x2 vector of responses for patients choosing treatment 2
#' @param s22 vector of stratum membership for patients choosing treatment 2. 
#'            Should be a vector of the same length as x2 with the number of
#'            unique values equal to the number of strata.
#' @param y1 vector of responses for patients randomized to treatment 1
#' @param s1 vector of stratum membership for patients randomized to treatment 1. 
#'            Should be a vector of the same length as y1 with the number of
#'            unique values equal to the number of strata.
#' @param y2 vector of responses for patients randomized to treatment 2.
#' @param s2 vector of stratum membership for patients randomized to treatment 2. 
#'            Should be a vector of the same length as y2 with the number of
#'            unique values equal to the number of strata.
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1.
#' @param nstrata number of strata (default=2)
#' @examples
#' x1<-c(10,8,6,10,5)
#' s11<-c(1,1,2,2,2)
#' x2<-c(8,7,6,10,12,11,6,8)
#' s22<-c(1,1,1,1,2,2,2)
#' y1<-c(10,5,7,9,12,6)
#' s1<-c(1,1,1,2,2,2)
#' y2<-c(8,9,10,7,8,11)
#' s2<-c(1,1,1,2,2,2)
#' @export
strat_analysis<-function(x1,s11,x2,s22,y1,s1,y2,s2,xi=c(0.5,0.5),nstrata=2){
  # Error messages
  if(!is.numeric(x1) | !is.numeric(x1) | !is.numeric(y1) | !is.numeric(y2))
    stop("Arguments must be numeric vectors")
  if(length(s11)!=length(x1))
    stop("Length of s11, x1 must match")
  if(length(s22)!=length(x2))
    stop("Length of s22, x2 must match")
  if(length(s1)!=length(y1))
    stop("Length of s1, y1 must match")
  if(length(s2)!=length(y2))
    stop("Length of s2, y2 must match")
  if(length(unique(s11))!=nstrata | length(unique(s22))!=nstrata | 
     length(unique(s11))!=nstrata | length(unique(s11))!=nstrata)
    stop("Number of unique elements in strata membership not equal to nstrata")
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata))
    stop('Number of strata must be numeric greater than 0')
  
  unstrat_stats<-matrix(NA,nrow=nstrata,ncol=3)
  
  # Compute unstratified test statistics
  for(i in 1:nstrata){
    x1i<-x1[as.factor(s11)==levels(as.factor(s11))[i]]
    x2i<-x2[as.factor(s22)==levels(as.factor(s22))[i]]
    y1i<-y1[as.factor(s1)==levels(as.factor(s1))[i]]
    y2i<-y2[as.factor(s2)==levels(as.factor(s2))[i]]
    
    unstrat_stats[i,]<-unlist(unstrat_analysis(x1i,x2i,y1i,y2i)[c(1,3,5)])
  }
  
  # Compute stratified test statistics and p-values
  pref_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,1]))
  sel_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,2]))
  treat_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,3]))
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval<-(1-pnorm(abs(pref_test)))*2 # Preference effect
  sel_pval<-(1-pnorm(abs(sel_test)))*2 # Selection effect
  treat_pval<-(1-pnorm(abs(treat_test)))*2
  
  results<-data.frame(pref_test,pref_pval,sel_test,sel_pval,treat_test,treat_pval)
  
  return(results)
}

################################
#### UNSTRATIFIED FUNCTIONS ####
################################

#' Unstratified Selection Effect Sample Size
#'
#' Calculates the sample size required to detect a given selection effect 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_selection(power=0.8, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5)
#' @export
unstrat_selection<-function(power, phi, sigma2, delta_pi, delta_nu, 
                            alpha=0.05, theta=0.5) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
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
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_preference(power=0.8, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5)
#' @export
unstrat_preference<-function(power, phi, sigma2, delta_pi, delta_nu, 
                             alpha=0.05, theta=0.5) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
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
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_treatment(power=0.8, sigma2=1, delta_tau=0.5)
#' @export
unstrat_treatment<-function(power, sigma2, delta_tau, alpha=0.05, 
                            theta=0.5){
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
    stop('Power must be numeric in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate sample size
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  N=(4*sigma2*(zbeta+zalpha)^2)/((1-theta)*delta_tau^2)
  return(N)
}

#' Unstratified Overall Sample Size
#'
#' Calculates the sample size required to detect a set of effects 
#' in an unstratified two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_overall(power=0.8, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5,
#' delta_tau=0.5)
#' @export
unstrat_overall<-function(power, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                            alpha=0.05, theta=0.5) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power)) 
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
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  sel=unstrat_selection(power, phi, sigma2, delta_pi, delta_nu, alpha, theta)
  pref=unstrat_preference(power, phi, sigma2, delta_pi, delta_nu, alpha, 
                          theta)
  treat=unstrat_treatment(power, sigma2, delta_tau, alpha, theta)
  return(max(sel, pref, treat))
}

#' Unstratified Treatment Effect Power Calculation
#'
#' Calculates the study power to detect the treatment effect given a particular 
#' sample size in an unstratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_trt_pwr(N=300, sigma2=1, delta_tau=0.5)
#' @export
unstrat_trt_pwr<-function(N, sigma2, delta_tau, alpha=0.05, theta=0.5) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau))
    stop('Effect size must be numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  power=pnorm(sqrt(((1-theta)*delta_tau^2*N)/(4*sigma2))-zalpha)

  return(power)
}

#' Unstratified Preference Effect Power Calculation
#'
#' Calculates the study power to detect the preference effect given a particular 
#' sample size in an unstratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_pref_pwr(N=300, phi=0.5, sigma2=1, delta_pi=0.5, delta_nu=0.5)
#' @export
unstrat_pref_pwr<-function(N, phi, sigma2, delta_pi, delta_nu, 
                          alpha=0.05, theta=0.5) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  numerator=4*theta*phi^2*(1-phi)^2*delta_pi^2*N
  denominator=sigma2+phi*(1-phi)*((2*phi-1)*delta_pi+delta_nu)^2+2*
    (theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  power=pnorm(sqrt(numerator/denominator)-zalpha)
  
  return(power)
}

#' Unstratified Selection Effect Power Calculation
#'
#' Calculates the study power to detect the selection effect given a particular 
#' sample size in an unstratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_sel_pwr(N=300, phi=0.5, sigma2=1, delta_pi=0.5, delta_nu=0.5)
#' @export
unstrat_sel_pwr<-function(N, phi, sigma2, delta_pi, delta_nu, 
                           alpha=0.05, theta=0.5) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  numerator=4*theta*phi^2*(1-phi)^2*delta_nu^2*N
  denominator=sigma2+phi*(1-phi)*((2*phi-1)*delta_nu+delta_pi)^2+2*
    (theta/(1-theta))*(phi^2+(1-phi)^2)*sigma2
  power=pnorm(sqrt(numerator/denominator)-zalpha)
  
  return(power)
}

#' Unstratified Power Calculation from Sample Size
#'
#' Calculates the study power to detect a set of effects given a particular 
#' sample size in an unstratified two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_power(N=300, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5,
#' delta_tau=0.5)
#' @export
unstrat_power<-function(N, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                          alpha=0.05, theta=0.5) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
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
  
  # Calculate study power
  trt_pwr<-unstrat_trt_pwr(N,sigma2,delta_tau,alpha,theta)
  pref_pwr<-unstrat_pref_pwr(N,phi,sigma2,delta_pi,delta_nu,alpha,theta)
  sel_pwr<-unstrat_sel_pwr(N,phi,sigma2,delta_pi,delta_nu,alpha,theta)

  return(data.frame(trt_pwr=trt_pwr,pref_pwr=pref_pwr,sel_pwr=sel_pwr))  
}

#' Unstratified Treatment Effect Back Calculation
#' 
#' Calculates the treatment effect that can be detected given a desired study 
#' power and overall study sample size for the unstratified two-stage 
#' randomized design
#' 
#' @param N overall study sample size.
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @examples
#' unstrat_trt_effect(N=300,power=0.9,sigma2=1)
#' @export
unstrat_trt_effect<-function(N, power, sigma2, alpha=0.05, theta=0.5) {
  # Error messages
  if(N<0 | !is.numeric(N)) 
    stop('N must be a positive numeric value')
  if(power<0 | power>1 | !is.numeric(power)) 
    stop('Power must be numeric in [0,1]')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(alpha<0 | alpha>1 | !is.numeric(alpha))
    stop('Type I error rate must be numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta)) 
    stop('Theta must be numeric in [0,1]')
  
  # Calculate effect size
  zbeta=qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  effect=sqrt((4*sigma2*(zbeta+zalpha)^2)/((1-theta)*N))
  
  return(effect)
}

#' Unstratified Optimized Theta
#'
#' Calculates the optimal proportion of patients assigned to the choice arm
#' in a two-stage randomized trial
#'
#' @param w_sel weight assigned to the estimation of the selection effect. Each 
#'              weight should be a numeric value between 0 and 1 and sum of 
#'              three weights should be 1.
#' @param w_pref weight assigned to the estimation of the preference effect. Each 
#'              weight should be a numeric value between 0 and 1 and sum of 
#'              three weights should be 1.
#' @param w_treat weight assigned to estimation of the treatment effect. Each 
#'              weight should be a numeric value between 0 and 1 and sum of 
#'              three weights should be 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect. 
#' @examples
#' theta_optim(w_sel=0.2, w_pref=0.4, w_treat=0.4, sigma1=1, phi=0.5,
#' delta_pi=1, delta_nu=0.5)
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

#' Unstratified Analysis Function
#'
#' Computes the test statistic and p-value for the preference, selection, and 
#' treatment effects for the unstratified two-stage randomized trial
#'
#' @param x1 vector of responses for patients choosing treatment 1
#' @param x2 vector of responses for patients choosing treatment 2
#' @param y1 vector of responses for patients randomized to treatment 1
#' @param y2 vector of responses for patients randomized to treatment 2
#' @examples
#' x1<-c(10,9,5,12,14)
#' x2<-c(16,12,14,10,8,9,11)
#' y1<-c(10,6,5,14,8,10)
#' y2<-c(12,10,15,11,11,13)
#' unstrat_analysis(x1,x2,y1,y2)
#' @export
unstrat_analysis<-function(x1,x2,y1,y2) {
  # Error messages
  if(!is.numeric(x1) | !is.numeric(x1) | !is.numeric(y1) | !is.numeric(y2))
    stop("Arguments must be numeric vectors")
  
  # Define sample sizes
  m1<-length(x1)
  m2<-length(x2)
  n1<-length(y1)
  n2<-length(y2)
  m<-m1+m2
  n<-n1+n2
  N<-m+n
  
  # Calculate z statistic
  z1<-sum(x1)-m1*mean(y1)
  z2<-sum(x2)-m2*mean(y2)
  
  # Calculate variances (formulas from Rucker paper)
  var1<-m1*var(x1)+(1+((m-1)/m)*m1)*m1*(var(y1)/n1)+
    (m1*m2/m)*(mean(x1)-mean(y1))^2
  var2<-m2*var(x2)+(1+((m-1)/m)*m2)*m2*(var(y2)/n2)+
    (m1*m2/m)*(mean(x2)-mean(y2))^2
  cov<--(m1*m2/m)*(mean(x1)-mean(y1))*(mean(x2)-mean(y2))
  
  # Compute test statistics (from Rucker paper)
  pref_test<-(z1+z2)/sqrt(var1+var2+2*cov) # Preference effect
  sel_test<-(z1-z2)/sqrt(var1+var2-2*cov) # Selection effect
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval<-(1-pnorm(abs(pref_test)))*2 # Preference effect
  sel_pval<-(1-pnorm(abs(sel_test)))*2 # Selection effect
  
  # Compute treatment effect t-test from random arm
  treat_test<-t.test(y1,y2)$statistic
  treat_pval<-t.test(y1,y2)$p.value
  
  results<-data.frame(pref_test,pref_pval,sel_test,sel_pval,treat_test,treat_pval)
  
  return(results)
}

#########################
#### OTHER FUNCTIONS ####
#########################

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
#'            entry corresponding to stratum-specific preference rate. All 
#'            elements should be numeric values between 0 and 1.
#' @param nstrata number of strata. Default is 1 (unstratified design).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. Should only be specified for stratified
#'          design.
#' @examples
#' calc_effects(mu1=1, mu2=2, mu11=1.5, mu22=2.5, phi=0.5)
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


