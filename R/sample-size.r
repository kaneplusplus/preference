############################
### SAMPLE SIZE FORMULAS ###
############################

#' Selection Effect Sample Size
#'
#' Calculates the sample size required to detect a given selection effect 
#' in a two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1. 
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If 
#'               study is stratified, should be vector of within-stratum 
#'               variances with length equal to the number of strata in the 
#'               study. 
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' selection_sample_size(power=0.8, phi=0.6, sigma2=1, delta_pi=1, delta_nu=0.5)
#' # Stratified
#' selection_sample_size(power=0.8, phi=c(0.5, 0.5), sigma2=c(1,1), delta_pi=1, 
#'  delta_nu=0.5,xi=c(0.3,0.7),nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @importFrom stats qnorm
#' @export
selection_sample_size<-function(power, phi, sigma2, delta_pi, delta_nu, 
                alpha=0.05, theta=0.5, xi=1, 
                nstrata=1) {
  # Error messages
  if(!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Power must be single numeric value in [0,1]')
  if (length(phi) != nstrata) 
    stop('Length vector does not match number of strata')
  if(any(!is.numeric(phi)) || any(phi <= 0) || any(phi >= 1))
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2) != nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(!is.numeric(sigma2) || any(sigma2 <= 0)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) || !is.numeric(delta_nu) ||
     length(delta_pi)!=1 || length(delta_nu)!=1)
    stop('Effect size must be single numeric value')
  if(!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be single numeric value in [0,1]')
  if(!is.numeric(theta) || theta <= 0 || theta >= 1 || length(theta)!=1)
    stop('Theta must be single numeric value in [0,1]')
  if(any(!is.numeric(xi) || any(xi <= 0) || any(xi > 1)))
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi) != nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi) != 1) 
    stop('Stratum proportions do not sum to 1')
  if(!is.numeric(nstrata) || nstrata <=0 || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
               *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_nu+delta_pi)^2
                 +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_nu^2)*sum_total
  return(ceiling(N))
}

#' Preference Effect Sample Size
#'
#' Calculates the sample size required to detect a given preference effect 
#' in a two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' preference_sample_size(power=0.8, phi=0.6, sigma2=1, delta_pi=1, delta_nu=0.5)
#' # Stratified
#' preference_sample_size(power=0.8, phi=c(0.5, 0.5), sigma2=c(1, 1), delta_pi=1, 
#'  delta_nu=0.5,xi=c(0.3,0.7),nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
preference_sample_size<-function(power, phi, sigma2, delta_pi, delta_nu, 
                 alpha=0.05, theta=0.5, xi=1, nstrata=1) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power) || length(power)!=1) 
    stop('Power must be single numeric value in [0,1]')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) ||
     length(delta_pi)!=1 || length(delta_nu)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) | length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) (xi[x]/(phi[x]^2*(1-phi[x])^2))
               *(sigma2[x]+phi[x]*(1-phi[x])*((2*phi[x]-1)*delta_pi+delta_nu)^2
                 +2*(theta/(1-theta))*sigma2[x]*(phi[x]^2+(1-phi[x])^2)))
  sum_total=sum(terms)
  N=(zalpha+zbeta)^2/(4*theta*delta_pi^2)*sum_total
  return(ceiling(N))
}

#' Treatment Effect Sample Size
#'
#' Calculates the sample size required to detect a given treatment effect 
#' in a two-stage randomized clinical trial
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' treatment_sample_size(power=0.8, sigma2=1, delta_tau=1.5)
#' # Stratified
#' treatment_sample_size(power=0.8, sigma2=c(1, 1), delta_tau=1.5, 
#'                       xi=c(0.3,0.7),
#' nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." 
#' \emph{Stat Methods Med Res}.  
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
treatment_sample_size<-function(power, sigma2, delta_tau, alpha=0.05,theta=0.5, 
                                xi=1, nstrata=1) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power) || length(power)!=1) 
    stop('Power must be single numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau) || length(delta_tau)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  #Calculate sample size
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  terms=sapply(1:nstrata, function(x) xi[x]*sigma2[x])
  sum_total=sum(terms)
  N=4*(zalpha+zbeta)^2/((1-theta)*delta_tau^2)*sum_total
  return(ceiling(N))
}

#' Overall Sample Size
#'
#' Calculates the sample size required to detect a given set of effects 
#' in a two-stage randomized clinical trial. Returns the largest
#' of the required sample sizes for a given set of treatment, selection, and 
#' preference effects.
#'
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
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
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' overall_sample_size(power=0.8, phi=0.5, sigma2=1, delta_pi=1, delta_nu=0.5, 
#' delta_tau=1.5)
#' # Stratified
#' overall_sample_size(power=0.8, phi=c(0.5,0.4), sigma2=c(1, 1), delta_pi=1, 
#' delta_nu=0.5, delta_tau=1.5, xi=c(0.3,0.7),nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
overall_sample_size<-function(power, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                    alpha=0.05, theta=0.5, xi=1, nstrata=1) {
  # Error messages
  if(power<0 | power>1 | !is.numeric(power) || length(power)!=1) 
    stop('Power must be single numeric value in [0,1]')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) | !is.numeric(delta_tau) 
     || length(delta_pi)!=1 || length(delta_nu)!=1 || length(delta_tau)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate sample size
  zbeta<-qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  pref=preference_sample_size(power, phi, sigma2, delta_pi, delta_nu, alpha, 
                              theta, xi, nstrata)
  sel=selection_sample_size(power, phi, sigma2, delta_pi, delta_nu, alpha, theta, 
                            xi, nstrata)
  treat=treatment_sample_size(power, sigma2, delta_tau, alpha, theta, xi, 
                              nstrata)
  
  ss<-list("treatment"=treat,"selection"=sel,"preference"=pref)
  return(ss)
}

###################################
### POWER CALCULATION FUNCTIONS ###
###################################

#' Treatment Effect Power Calculation
#'
#' Calculates the study power to detect the treatment effect given a particular 
#' sample size in a two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param delta_tau overall study treatment effect.
#' @param alpha desired type I error rate.. 
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' treatment_power(N=300, sigma2=1, delta_tau=0.5)
#' # Stratified
#' treatment_power(N=300, sigma2=c(1,1), delta_tau=0.5, xi=c(0.5,0.5), nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
treatment_power<-function(N, sigma2, delta_tau, alpha=0.05, theta=0.5, xi=1, 
                  nstrata=1) {
  # Error messages
  if(N<0 | !is.numeric(N) | length(N)!=1) 
    stop('N must be a single positive numeric value')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_tau) || length(delta_tau)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate study power
  zalpha<-qnorm(1-(alpha/2))
  power=pnorm(sqrt((((1-theta)*delta_tau^2*N)/(4*sum(sapply(1:nstrata,function(i)
    xi[i]*sigma2[i])))))-zalpha)
  
  return(power)
}

#' Preference Effect Power Calculation
#'
#' Calculates the study power to detect the preference effect given a particular 
#' sample size in a two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' preference_power(N=300, phi=0.6, sigma2=1, delta_pi=1, delta_nu=0.5)
#' # Stratified
#' preference_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, 
#' delta_nu=0.5, xi=c(0.5,0.5), nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
preference_power<-function(N, phi, sigma2, delta_pi, delta_nu, alpha=0.05, 
                   theta=0.5, xi=1, nstrata=1) {
  # Error messages
  if(N<0 | !is.numeric(N) | length(N)!=1) 
    stop('N must be a single positive numeric value')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) ||
     length(delta_pi)!=1 || length(delta_nu)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
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

#' Selection Effect Power Calculation
#'
#' Calculates the study power to detect the selection effect given a particular 
#' sample size in a two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' selection_power(N=300, phi=0.6, sigma2=1, delta_pi=1, delta_nu=0.5)
#' # Stratified
#' selection_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, 
#' delta_nu=0.5, xi=c(0.5,0.5), nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
selection_power<-function(N, phi, sigma2, delta_pi, delta_nu, alpha=0.05, 
                  theta=0.5, xi=1, nstrata=1) {
  # Error messages
  if(N<0 | !is.numeric(N) | length(N)!=1) 
    stop('N must be a single positive numeric value')
  if (length(phi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) ||
     length(delta_pi)!=1 || length(delta_nu)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
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

#' Power Calculation from Sample Size
#'
#' Calculates the study power to detect a set of effects given a particular 
#' sample size in a two-stage randomized clinical trial
#'
#' @param N overall study sample size.
#' @param phi the proportion of patients preferring treatment 1. Should be
#'            numeric value between 0 and 1. If study is stratified, should be
#'            vector with length equal to the number of strata in the study.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
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
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' # Unstratified
#' overall_power(N=300, phi=0.6, sigma2=1, delta_pi=1, delta_nu=0.5, 
#' delta_tau=1.5)
#' # Stratified
#' overall_power(N=300, phi=c(0.6,0.5), sigma2=c(1,1), delta_pi=1, delta_nu=0.5, 
#' delta_tau=0.5, xi=c(0.5,0.5), nstrata=2)
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and 
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
overall_power<-function(N, phi, sigma2, delta_pi, delta_nu, delta_tau, 
                      alpha=0.05, theta=0.5, xi=1, nstrata=1) {
  # Error messages
  if(N<0 | !is.numeric(N) | length(N)!=1) 
    stop('N must be a single positive numeric value')
  if (length(phi)!=nstrata) 
    stop('Length vector does not match number of strata')
  if(any(phi<0) | any(phi>1) | any(!is.numeric(phi))) 
    stop('Preference rate must be numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) | !is.numeric(delta_tau) 
     || length(delta_pi)!=1 || length(delta_nu)!=1 || length(delta_tau)!=1)
    stop('Effect size must be single numeric value')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be single numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  
  # Calculate study power
  trt_pwr<-treatment_power(N=N,sigma2=sigma2,delta_tau=delta_tau,alpha=alpha,
                   theta=theta,xi=xi,nstrata=nstrata)
  pref_pwr<-preference_power(N=N,phi=phi,sigma2=sigma2,delta_pi=delta_pi,
                     delta_nu=delta_nu,alpha=alpha,theta=theta,xi=xi,
                     nstrata=nstrata)
  sel_pwr<-selection_power(N=N,phi=phi,sigma2=sigma2,delta_pi=delta_pi,
                   delta_nu=delta_nu,alpha=alpha,theta=theta,xi=xi,
                   nstrata=nstrata)
  
  return(list("trt_pwr"=trt_pwr,"pref_pwr"=pref_pwr,"sel_pwr"=sel_pwr))  
}

##########################
### ANALYSIS FUNCTIONS ###
##########################

#' Analysis Function: Raw Data
#' 
#' Computes the test statistic and p-value for the preference, selection, and 
#' treatment effects for the two-stage randomized trial using
#' provided raw data
#' 
#' @param x1 vector of responses for patients choosing treatment 1
#' @param x2 vector of responses for patients choosing treatment 2
#' @param y1 vector of responses for patients randomized to treatment 1
#' @param y2 vector of responses for patients randomized to treatment 2.
#' @param s11 (if study is stratified) vector of stratum membership for 
#'            patients choosing treatment 1. Should be a vector of the same 
#'            length as x1 with the number of unique values equal to the 
#'            number of strata.
#' @param s22 (if study is stratified) vector of stratum membership for 
#'            patients choosing treatment 2. Should be a vector of the same 
#'            length as x2 with the number of unique values equal to the number
#'            of strata.
#' @param s1  (if study is stratified) vector of stratum membership for 
#'            patients randomized to treatment 1. Should be a vector of the 
#'            same length as y1 with the number of unique values equal to the 
#'            number of strata.
#' @param s2 (if study is stratified) vector of stratum membership for patients
#'           randomized to treatment 2. Should be a vector of the same length 
#'           as y2 with the number of unique values equal to the number of 
#'           strata.
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' #Unstratified
#' x1<-c(10,8,6,10,5)
#' x2<-c(8,7,6,10,12,11,6,8)
#' y1<-c(10,5,7,9,12,6)
#' y2<-c(8,9,10,7,8,11)
#' analyze_raw_data(x1,x2,y1,y2)
#' #Stratified
#' x1<-c(10,8,6,10,5)
#' s11<-c(1,1,2,2,2)
#' x2<-c(8,7,6,10,12,11,6,8)
#' s22<-c(1,1,1,1,2,2,2,2)
#' y1<-c(10,5,7,9,12,6)
#' s1<-c(1,1,1,2,2,2)
#' y2<-c(8,9,10,7,8,11)
#' s2<-c(1,1,1,2,2,2)
#' analyze_raw_data(x1,x2,y1,y2,s11=s11,s22=s22,s1=s1,s2=s2,xi=c(0.5,0.5),
#' nstrata=2)
#' @references Rucker G (1989). "A two-stage trial design for testing treatment, 
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
analyze_raw_data<-function(x1,x2,y1,y2,s11=1,s22=1,s1=1,s2=1,xi=1,nstrata=1){
  # Check stratum assignments
  if(nstrata==1){
    s11=rep(1,length(x1))
    s22=rep(1,length(x2))
    s1=rep(1,length(y1))
    s2=rep(1,length(y2))
  }
  
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
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  unstrat_stats<-matrix(NA,nrow=nstrata,ncol=6)
  
  # Compute unstratified test statistics
  for(i in 1:nstrata){
    x1i<-x1[as.factor(s11)==levels(as.factor(s11))[i]]
    x2i<-x2[as.factor(s22)==levels(as.factor(s22))[i]]
    y1i<-y1[as.factor(s1)==levels(as.factor(s1))[i]]
    y2i<-y2[as.factor(s2)==levels(as.factor(s2))[i]]
    
    unstrat_stats[i,]<-unlist(unstrat_analyze_raw_data(x1i,x2i,y1i,y2i))
  }
  
  # Compute stratified test statistics and p-values
  pref_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,1]))
  sel_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,3]))
  treat_test<-sum(sapply(1:nstrata, function(i) xi[i]*unstrat_stats[i,5]))
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval<-pnorm(abs(pref_test), lower.tail = FALSE)*2 # Preference effect
  sel_pval<-pnorm(abs(sel_test), lower.tail = FALSE)*2 # Selection effect
  treat_pval<-pnorm(abs(treat_test), lower.tail=FALSE)*2
  
  results<-data.frame(pref_test, pref_pval, sel_test, sel_pval, 
                      treat_test, treat_pval)
  
  return(results)
}

#' Analysis Function: Summary Data
#' 
#' Computes the test statistic and p-value for the preference, selection, and 
#' treatment effects for the two-stage randomized trial using provided summary 
#' data
#' 
#' @param x1mean mean of responses for patients choosing treatment 1. If study
#'               is stratified, should be vector with length equal to the
#'               number of strata.
#' @param x1var variance of responses for patients choosing treatment 1. If 
#'              study is stratified, should be vector with length equal to the
#'              number of strata.
#' @param m1 number of patients choosing treatment 1. If study
#'               is stratified, should be vector with length equal to the
#'               number of strata.
#' @param x2mean mean of responses for patients choosing treatment 2. If study
#'               is stratified, should be vector with length equal to the
#'               number of strata.
#' @param x2var variance of responses for patients choosing treatment 2. If 
#'              study is stratified, should be vector with length equal to the
#'              number of strata.
#' @param m2 number of patients choosing treatment 2. If study
#'               is stratified, should be vector with length equal to the
#'               number of strata.
#' @param y1mean mean of responses for patients randomized to treatment 1. If 
#'               study is stratified, should be vector with length equal to the
#'               number of strata.
#' @param y1var variance of responses for patients randomized to treatment 1. If 
#'              study is stratified, should be vector with length equal to the
#'              number of strata.
#' @param n1 number of patients randomized to treatment 1. If study is 
#'           stratified, should be vector with length equal to the number of 
#'           strata.
#' @param y2mean mean of responses for patients randomized to treatment 2. If 
#'               study is stratified, should be vector with length equal to the
#'               number of strata.
#' @param y2var variance of responses for patients randomized to treatment 2. If 
#'              study is stratified, should be vector with length equal to the
#'              number of strata.
#' @param n2 number of patients randomized to treatment 2. If study is 
#'           stratified, should be vector with length equal to the number of 
#'           strata.
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' x1mean<-5
#' x1var<-1
#' m1<-15
#' x2mean<-7
#' x2var<-1.1
#' m2<-35
#' y1mean<-6
#' y1var<-1
#' n1<-25
#' y2mean<-8
#' y2var<-1.2
#' n2<-25
#' analyze_summary_data(x1mean,x2var,m1,x2mean,x2var,m2,y1mean,y2var,n1,y2mean,
#' y2var,n2)
#' @references Rucker G (1989). "A two-stage trial design for testing treatment,
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
analyze_summary_data<-function(x1mean,x1var,m1,x2mean,x2var,m2,y1mean,y1var,
                           n1,y2mean,y2var,n2,xi=1,nstrata=1){
  # Error messages
  if(!is.numeric(x1mean) | !is.numeric(x1var) | 
     !is.numeric(x2mean) | !is.numeric(x2var) |
     !is.numeric(y1mean) | !is.numeric(y1var) |
     !is.numeric(y2mean) | !is.numeric(y2var) |
     !is.numeric(m1) | !is.numeric(m2) | !is.numeric(n1) | !is.numeric(n2))
    stop("Arguments must be numeric vectors")
  if(length(x1mean)!=nstrata | length(x1var)!=nstrata | length(m1)!=nstrata)
    stop("Length of vector must match number of strata")
  if(length(x2mean)!=nstrata | length(x2var)!=nstrata | length(m2)!=nstrata)
    stop("Length of vector must match number of strata")
  if(length(y1mean)!=nstrata | length(y1var)!=nstrata | length(n1)!=nstrata)
    stop("Length of vector must match number of strata")
  if(length(y2mean)!=nstrata | length(y2var)!=nstrata | length(n2)!=nstrata)
    stop("Length of vector must match number of strata")
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Compute unstratified test statistics
  unstrat_stats<-sapply(1:nstrata, function(i) 
    unstrat_analyze_summary_data(x1mean[i],x1var[i],m1[i],x2mean[i],x2var[i],
                             m2[i],y1mean[i],y1var[i],n1[i],y2mean[i],y2var[i],
                             n2[i]))
  
  # Compute stratified test statistics and p-values
  pref_test<-sum(sapply(1:nstrata, function(i) xi[i]*unlist(unstrat_stats[1,i])))
  sel_test<-sum(sapply(1:nstrata, function(i) xi[i]*unlist(unstrat_stats[3,i])))
  treat_test<-sum(sapply(1:nstrata,function(i) xi[i]*unlist(unstrat_stats[5,i])))
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval<-pnorm(abs(pref_test), lower.tail = FALSE)*2 # Preference effect
  sel_pval<-pnorm(abs(sel_test), lower.tail = FALSE)*2 # Selection effect
  treat_pval<-pnorm(abs(treat_test), lower.tail = FALSE)*2
  
  results<-data.frame(pref_test, pref_pval, sel_test, sel_pval, treat_test, 
                      treat_pval)
  
  return(results)
}

######################
### MISC FUNCTIONS ###
######################

#' Treatment Effect Back Calculation
#' 
#' Calculates the treatment effect that can be detected given a desired study 
#' power and overall study sample size for the two-stage randomized design
#' 
#' @param N overall study sample size.
#' @param power desired study power. Should be numeric value between 0 and 1.
#' @param sigma2 variance estimate. Should be positive numeric values. If study
#'               is stratified, should be vector of within-stratum variances 
#'               with length equal to the number of strata in the study.
#' @param alpha desired type I error rate.
#' @param theta proportion of patients assigned to choice arm in the initial
#'              randomization. Should be numeric value between
#'              0 and 1 (default=0.5).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. All vector elements should be numeric 
#'          values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @examples
#' treatment_effect_size(N=300,power=0.9,sigma2=c(1,0.8), xi=c(0.3,0.7), 
#' nstrata=2)
#' @export
treatment_effect_size<-function(N, power, sigma2, alpha=0.05, theta=0.5, xi=1, 
                     nstrata=1) {
  # Error messages
  if(N<0 | !is.numeric(N) | length(N)!=1) 
    stop('N must be a single positive numeric value')
  if(power<0 | power>1 | !is.numeric(power) || length(power)!=1) 
    stop('Power must be single numeric value in [0,1]')
  if(length(sigma2)!=nstrata)
    stop('Length of variance vector does not match number of strata')
  if(any(sigma2<=0) | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(alpha<0 | alpha>1 | !is.numeric(alpha) || length(alpha)!=1)
    stop('Type I error rate must be alpha numeric in [0,1]')
  if(theta<0 | theta>1 | !is.numeric(theta) || length(theta)!=1) 
    stop('Theta must be single numeric in [0,1]')
  if(any(xi<0) | any(xi>1) | any(!is.numeric(xi))) 
    stop('Proportion of patients in strata must be numeric value in [0,1]')
  if (length(xi)!=nstrata) 
    stop('Length of vector does not match number of strata')
  if (sum(xi)!=1) 
    stop('Stratum proportions do not sum to 1')
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
    stop('Number of strata must be numeric greater than 0')
  
  # Calculate effect size
  zbeta=qnorm(power)
  zalpha<-qnorm(1-(alpha/2))
  effect=sqrt(((4*(zbeta+zalpha)^2)/((1-theta)*N))*
                sum(sapply(1:nstrata, function(i) xi[i]*sigma2[i])))
  
  return(effect)
}

#' Unstratified Optimized Theta
#'
#' Calculates the optimal proportion of patients assigned to the choice arm
#' in an unstratified two-stage randomized trial
#'
#' @param w_sel weight assigned to the estimation of the selection effect. Each 
#'              weight should be a numeric value between 0 and 1 and sum of 
#'              three weights should be 1.
#' @param w_pref weight assigned to the estimation of the preference effect. 
#'               Each weight should be a numeric value between 0 and 1 and 
#'               sum of three weights should be 1.
#' @param w_treat weight assigned to estimation of the treatment effect. Each 
#'              weight should be a numeric value between 0 and 1 and sum of 
#'              three weights should be 1.
#' @param sigma2 variance estimate. Should be a positive numeric value.
#' @param phi proportion of patients preferring treatment 1. Should be numeric
#'            value between 0 and 1.
#' @param delta_pi overall study preference effect.
#' @param delta_nu overall study selection effect. 
#' @examples
#' optimal_proportion(w_sel=0.2, w_pref=0.4, w_treat=0.4, sigma2=1, phi=0.5,
#' delta_pi=1, delta_nu=0.5)
#' @references Walter et. al. (2011). "Optimal allocation of participants for
#' the estimation of selection, preference and treatment effects in the 
#' two-stage randomised trial design." \emph{Stat Med}, 
#' \strong{31}(13):1307-1322.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22362374}{PubMed})
#' @importFrom stats uniroot
#' @export
optimal_proportion<-function(w_sel,w_pref,w_treat,sigma2,phi,delta_pi,delta_nu) {
  if(w_sel<0 | w_sel>1 | w_pref<0 | w_pref>1 | w_treat<0 | w_treat>1 | 
     any(!is.numeric(c(w_sel,w_pref,w_treat))) | length(w_sel)!=1 | 
     length(w_pref)!=1 | length(w_treat)!=1)
    stop('Weights must be single numeric value in [0,1]')
  if (w_sel+w_pref+w_treat!=1) 
    stop('weights do not sum to 1')
  if(sigma2<=0 | any(!is.numeric(sigma2)))
    stop('Variance estimate must be numeric value greater than 0')
  if(phi<0 | phi>1 | !is.numeric(phi)) 
    stop('Preference rate must be numeric value in [0,1]')
  if(!is.numeric(delta_pi) | !is.numeric(delta_nu) ||
     length(delta_pi)!=1 || length(delta_nu)!=1)
    stop('Effect size must be single numeric value')
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
#'            entry corresponding to stratum-specific preference rate. All 
#'            elements should be numeric values between 0 and 1.
#' @param nstrata number of strata. Default is 1 (unstratified design).
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'          Length of vector should equal the number of strata in the study and 
#'          sum of vector should be 1. Should only be specified for stratified
#'          design.
#' @examples
#' effects_from_means(mu1=1, mu2=2, mu11=1.5, mu22=2.5, phi=0.5)
#' @references Rucker G (1989). "A two-stage trial design for testing treatment, 
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @export
effects_from_means<-function(mu1,mu2,mu11,mu22,phi,nstrata=1,xi=NULL) {
  # Error messages
  if(nstrata<=0 | !is.numeric(nstrata) || length(nstrata)!=1)
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
    effects<-list("treatment"=delta_tau,"selection"=delta_nu,
                  "preference"=delta_pi)
  } else {
    # Stratified case
    effects<-list("treatment"=sum(sapply(1:nstrata, function(x) 
                              phi[x]*delta_tau[x])),
                  "selection"=sum(sapply(1:nstrata, function(x) 
                              phi[x]*delta_nu[x])),
                  "preference"=sum(sapply(1:nstrata,function(x) 
                              phi[x]*delta_pi[x])))
  }

  return(effects)
}

######################################
### Extra (non-exported) functions ###
######################################

### Find sigma^2 for unstratified case
sigma2_mixtures<-function(sigma2,means,prop){
  # Overall mixture
  val=prop[1]*sigma2[1]+prop[2]*sigma2[2]+(prop[1]*means[1]^2+
      prop[2]*means[2]^2-(prop[1]*means[1]+prop[2]*means[2])^2)
  return(val)
}

### Find sigma^2 for each stratum based on means
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

### Find means from effects
means_stratum<-function(sigma,mu,tau,nu,pi,prop,theta=c(0.5,0.5)){
  mean_choice=sapply(1:2, function(x) mu[x]+tau[x]+nu[x]+pi[x])
  mean_random=sapply(1:2, function(x) mu[x]+tau[x])
  # Compute mean value for each stratum
  means=sapply(1:2, function(x) theta[x]*mean_choice[x]+
              (1-theta[x])*mean_random[x])
  return(means)
}

#### Analaysis Function (Raw Data)

#' @importFrom stats pnorm t.test var
unstrat_analyze_raw_data<-function(x1,x2,y1,y2) {
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
  pref_pval<-pnorm(abs(pref_test), lower.tail = FALSE)*2 # Preference effect
  sel_pval<-pnorm(abs(sel_test), lower.tail = FALSE)*2 # Selection effect
  
  # Compute treatment effect t-test from random arm
  treat_test<-t.test(y1,y2)$statistic
  treat_pval<-t.test(y1,y2)$p.value
  
  results<-data.frame(pref_test, pref_pval, sel_test, sel_pval, treat_test, 
                      treat_pval)
  
  return(results)
}


### Analysis Function (Summary Data)
unstrat_analyze_summary_data<-function(x1mean, x1var, m1, x2mean, x2var, m2, 
                                       y1mean, y1var,n1, y2mean, y2var, n2) {
  # Error messages
  if(!is.numeric(x1mean) | !is.numeric(x1var) | 
     !is.numeric(x2mean) | !is.numeric(x2var) |
     !is.numeric(y1mean) | !is.numeric(y1var) |
     !is.numeric(y2mean) | !is.numeric(y2var) |
     !is.numeric(m1) | !is.numeric(m2) | !is.numeric(n1) | !is.numeric(n2))
    stop("Arguments must be numeric vectors")
  
  # Define sample sizes
  m<-m1+m2
  n<-n1+n2
  N<-m+n
  
  # Calculate z statistic
  z1<-m1*x1mean-m1*y1mean
  z2<-m2*x2mean-m2*y2mean
  
  # Calculate variances (formulas from Rucker paper)
  var1<-m1*x1var+(1+((m-1)/m)*m1)*m1*(y1var/n1)+
    (m1*m2/m)*(x1mean-y1mean)^2
  var2<-m2*x2var+(1+((m-1)/m)*m2)*m2*(y2var/n2)+
    (m1*m2/m)*(x2mean-y2mean)^2
  cov<--(m1*m2/m)*(x1mean-y1mean)*(x2mean-y2mean)
  
  # Compute test statistics (from Rucker paper)
  pref_test<-(z1+z2)/sqrt(var1+var2+2*cov) # Preference effect
  sel_test<-(z1-z2)/sqrt(var1+var2-2*cov) # Selection effect
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval<-pnorm(abs(pref_test), lower.tail = FALSE)*2 # Preference effect
  sel_pval<-pnorm(abs(sel_test), lower.tail = FALSE)*2 # Selection effect
  
  # Compute treatment effect t-test from random arm
  treat_test<-t.test2(y1mean,y2mean,y1var,y2var,n1,n2)$t
  treat_pval<-t.test2(y1mean,y2mean,y1var,y2var,n1,n2)$p.value
  
  results<-data.frame(pref_test, pref_pval, sel_test, sel_pval, treat_test, 
                      treat_pval)
  
  return(results)
}


### T-test from summary data (null hypothesis of no difference, no assumption 
# of equal variances)
# m1,m1: sample means
# s1,s2: sample variances
# n1, n2: sample sizes

#' @importFrom stats pt
t.test2 <- function(m1,m2,s1,s2,n1,n2)
{
  se <- sqrt( (s1/n1) + (s2/n2) )
  # Welch-satterthwaite df
  df <- ( (s1/n1 + s2/n2)^2 )/( (s1/n1)^2/(n1-1) + (s2/n2)^2/(n2-1) )
  
  t <- (m1-m2)/se 
  dat <- data.frame(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Mean.Diff", "Std.Err", "t", "p.value")
  return(dat) 
}

