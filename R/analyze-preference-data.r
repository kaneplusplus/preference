##########################
### ANALYSIS FUNCTIONS ###
##########################

#' Analysis Function: Raw Data
#' 
#' Computes the test statistic and p-value for the preference, selection, and 
#' treatment effects for the two-stage randomized trial using
#' provided raw data.
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
#'
#' #Unstratified
#'
#' x1 <- c(10, 8, 6, 10, 5)
#' x2 <- c(8, 7, 6, 10, 12, 11, 6, 8)
#' y1 <- c(10, 5, 7, 9, 12, 6)
#' y2 <- c(8, 9, 10, 7, 8, 11)
#' analyze_raw_data(x1, x2, y1, y2)
#'
#' #Stratified
#'
#' x1 <- c(10, 8, 6, 10, 5)
#' s11 <- c(1, 1, 2, 2, 2)
#' x2 <- c(8, 7, 6, 10, 12, 11, 6, 8)
#' s22 <- c(1, 1, 1, 1, 2, 2, 2, 2)
#' y1 <- c(10, 5, 7, 9, 12, 6)
#' s1 <- c(1, 1, 1, 2, 2, 2)
#' y2 <- c(8, 9, 10, 7, 8, 11)
#' s2 <- c(1, 1, 1, 2, 2, 2)
#' analyze_raw_data(x1, x2, y1, y2, s11=s11, s22=s22, s1=s1, s2=s2,
#'                  xi=c(0.5,0.5), nstrata=2)
#' @references Rucker G (1989). "A two-stage trial design for testing treatment,
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." 
#' \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
analyze_raw_data <- function(x1, x2, y1, y2, s11=1, s22=1, s1=1, s2=1, xi=1, 
                             nstrata=1) {

  # Check stratum assignments
  if(nstrata == 1) {
    s11 <- rep(1,length(x1))
    s22 <- rep(1,length(x2))
    s1 <- rep(1,length(y1))
    s2 <- rep(1,length(y2))
  }
  
  # Error messages
  if(!is.numeric(x1) | !is.numeric(x1) | !is.numeric(y1) | !is.numeric(y2)) {
    stop("Arguments must be numeric vectors")
  }

  if(length(s11) != length(x1)) {
    stop("Length of s11, x1 must match")
  }

  if(length(s22) != length(x2)) {
    stop("Length of s22, x2 must match")
  }

  if(length(s1) != length(y1)) {
    stop("Length of s1, y1 must match")
  }

  if(length(s2) != length(y2)) {
    stop("Length of s2, y2 must match")
  }

  if(length(unique(s11)) != nstrata | length(unique(s22)) != nstrata | 
     length(unique(s11)) != nstrata | length(unique(s11)) != nstrata) {
    stop("Number of unique elements in strata membership not equal to nstrata")
  }

  if (length(xi) != nstrata) {
    stop('Length of vector does not match number of strata')
  }

  if (sum(xi) != 1) {
    stop('Stratum proportions do not sum to 1')
  }

  if(nstrata <=0 | !is.numeric(nstrata) || length(nstrata)!=1) {
    stop('Number of strata must be numeric greater than 0')
  }
  
  unstrat_stats <- matrix(NA,nrow=nstrata,ncol=6)
  
  # Compute unstratified test statistics
  for(i in 1:nstrata){
    x1i <- x1[as.factor(s11)==levels(as.factor(s11))[i]]
    x2i <- x2[as.factor(s22)==levels(as.factor(s22))[i]]
    y1i <- y1[as.factor(s1)==levels(as.factor(s1))[i]]
    y2i <- y2[as.factor(s2)==levels(as.factor(s2))[i]]
    
    unstrat_stats[i,] <- unlist(unstrat_analyze_raw_data(x1i, x2i, y1i, y2i))
  }
  
  # Compute stratified test statistics and p-values
  pref_test <- sum(vapply(seq_len(nstrata), 
                   function(i) xi[i] * unstrat_stats[i, 1], 0.0))
  sel_test <- sum(vapply(seq_len(nstrata), 
                  function(i) xi[i]*unstrat_stats[i, 3], 0.0))
  treat_test <- sum(vapply(seq_len(nstrata), 
                    function(i) xi[i]*unstrat_stats[i, 5], 0.0))
  
  # Compute p-values (Assume test stats approximately normally distributed)

  # preference effect
  pref_pval <- 2 * pnorm( abs( pref_test / sum(xi^2) ), lower.tail = FALSE ) 

  # selection effect
  sel_pval <- 2 * pnorm( abs( sel_test /sum(xi^2) ), lower.tail = FALSE )

  # treatment effect
  treat_pval <- 2 * pnorm( abs( treat_test / sum(xi^2) ), lower.tail=FALSE )
  
  data.frame(pref_test, pref_pval, sel_test, sel_pval, treat_test, treat_pval)
}

#' Fit the preference data collected from a clinical trial
#'
#' TODO fill this in.
#'
#' @param outcome (numeric) individual trial outcomes.
#' @param random (logical) was this individual part of the random arm?
#' @param treatment (character, factor, or integer) which treatment an 
#' individual received
#' @param strata (optional integer) which strata the individual belongs to.
#' @examples
#' ##Unstratified
#' outcome <- c(10, 8, 6, 10, 5, 8, 7, 6, 10, 12, 11, 6, 8, 10, 5, 7, 9, 12, 6,
#'              8, 9, 10, 7, 8,11)
#' random <- c(rep(FALSE, 13), rep(TRUE, 12))
#' treatment <- c(rep(1, 5), rep(2, 8), rep(1, 6), rep(2, 6))
#' fit_preference(outcome, random, treatment)
#' 
#' ##Stratified
#' #x1 <- c(10,8,6,10,5)
#' #s11 <- c(1,1,2,2,2)
#' #x2 <- c(8,7,6,10,12,11,6,8)
#' #s22 <- c(1,1,1,1,2,2,2,2)
#' #y1 <- c(10,5,7,9,12,6)
#' #s1 <- c(1,1,1,2,2,2)
#' #y2 <- c(8,9,10,7,8,11)
#' #s2 <- c(1,1,1,2,2,2)
#' #analyze_raw_data(x1, x2, y1, y2, s11=s11, s22=s22, s1=s1, s2=s2,
#' #                 xi=c(0.5,0.5), nstrata=2)
#' @export
fit_preference <- function(outcome, random, treatment, strata) {
  if (missing(strata)) {
    strata <- rep(1, length(outcome))
  }
  
  if (length(unique(treatment)) != 2) {
    stop("You may only have two treatments.")
  }

  pd <- data.frame(outcome=outcome, random=random, treatment=treatment,
                   strata=strata)
 
  # Compute unstratified test statistics
  strat_split <- split(seq_len(length(strata)), strata)
  
  treatments <- sort(unique(treatment)) 
  unstrat_stats <- NULL
  for (ss in strat_split) { 
    pds <- pd[ss , ]
    unstrat_stats <- rbind(unstrat_stats, 
      # unstrat_analyze_raw_data NEEDS TO CHANGE to include other values.
      unstrat_analyze_raw_data(
        pds$outcome[pds$random == FALSE & pds$treatment == treatments[1]],
        pds$outcome[pds$random == FALSE & pds$treatment == treatments[2]],
        pds$outcome[pds$random == TRUE & pds$treatment == treatments[1]],
        pds$outcome[pds$random == TRUE & pds$treatment == treatments[2]]))
  }
 
  #calculate xi
  xi <- table(strata) / length(outcome)
  nstrata <- length(unique(strata))
 
  # Compute stratified test statistics and p-values
  pref_test <- sum(vapply(seq_len(nstrata), 
                          function(i) xi[i] * unstrat_stats[i,1], 0.0))
  sel_test <- sum(vapply(seq_len(nstrata), 
                         function(i) xi[i] * unstrat_stats[i,3], 0.0))
  treat_test <- sum(vapply(seq_len(nstrata), 
                           function(i) xi[i] * unstrat_stats[i,5], 0.0))
 
  # Compute p-values (Assume test stats approximately normally distributed)
  
  # preference effect
  pref_pval <- 2 * pnorm(abs(pref_test), lower.tail = FALSE)

  # selection effect
  sel_pval <- 2 * pnorm(abs(sel_test), lower.tail = FALSE)

  # treatment effect
  treat_pval <- 2 * pnorm(abs(treat_test), lower.tail=FALSE)
 
  #Might need for unstrat_analyze_raw_data to also return the variance values
  data.frame(pref_test, pref_pval, sel_test, sel_pval, treat_test, treat_pval)
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
  
  t <- (m1 - m2)/se 
  
  dat <- data.frame(m1-m2, se, t, 2*pt(-abs(t),df))    
  
  names(dat) <- c("Mean.Diff", "Std.Err", "t", "p.value")
  
  dat
}



### Analysis Function (Summary Data)
#' @importFrom stats qnorm
unstrat_analyze_summary_data <- function(x1mean, x1var, m1, x2mean, x2var, m2, 
                                         y1mean, y1var, n1, y2mean, y2var, n2, 
                                         alpha) {
  # Define sample sizes
  m <- m1 + m2
  n <- n1 + n2
  N <- m + n
  
  # Calculate z values as defined by Rucker
  z1 <- m1*x1mean - m1*y1mean
  z2 <- m2*x2mean - m2*y2mean
  
  # Calculate variance components (formulas from Rucker paper)
  var1 <- m1 * x1var + 
    (1 + ((m - 1)/m)*m1)*m1*(y1var/n1) + (m1*m2/m)*(x1mean - y1mean)^2

  var2 <- m2 * x2var + 
    (1 + ((m - 1)/m)*m2)*m2*(y2var/n2) + (m1*m2/m)*(x2mean - y2mean)^2
  
  cov <- -(m1*m2/m)*(x1mean - y1mean)*(x2mean - y2mean)
  
  #Call the test2 function
  treat_out <- t.test2(y1mean, y2mean, y1var, y2var, n1, n2)
  
  # Compute effect estimates 
  pref_effect <- 0.5*(m/(m1*m2))*(z1+z2)
  sel_effect <- 0.5*(m/(m1*m2))*(z1-z2)
  treat_effect <- treat_out$Mean.Diff
  
  #Compute SE estimates 
  pref_SE <- sqrt(var1 + var2 + 2*cov)*0.5*(m/(m1*m2))
  sel_SE <- sqrt(var1 + var2 - 2*cov)*0.5*(m/(m1*m2))
  treat_SE <- treat_out$Std.Err
  
  # Compute test statistics 
  pref_test <- pref_effect/pref_SE 
  sel_test <- sel_effect/sel_SE 
  treat_test <- treat_out$t
  
  # Compute p-values (Assume test stats approximately normally distributed)
  pref_pval <- pnorm(abs(pref_test), lower.tail = FALSE)*2 # Preference effect
  sel_pval <- pnorm(abs(sel_test), lower.tail = FALSE)*2 # Selection effect
  treat_pval <- treat_out$p.value
  
  #Compute approximate (1-alpha)% confidence intervals
  zalpha <- qnorm(1 - (alpha/2))
  pref_LB <- pref_effect - zalpha*pref_SE
  pref_UB <- pref_effect + zalpha*pref_SE
  sel_LB <- sel_effect - zalpha*sel_SE
  sel_UB <- sel_effect + zalpha*sel_SE
  treat_LB <- treat_effect - zalpha*treat_SE
  treat_UB <- treat_effect + zalpha*treat_SE
  

  data.frame(pref_effect, pref_SE, pref_test, pref_pval,pref_LB, pref_UB,
             sel_effect, sel_SE, sel_test, sel_pval, sel_LB, sel_UB,
             treat_effect, treat_SE,  treat_test, treat_pval, treat_LB, 
             treat_UB)
}


#' Fit Preference Model from Summary Data
#' 
#' Computes the test statistic and p-value for the preference, 
#' selection, and treatment effects for the two-stage randomized trial using 
#' provided summary data.
#' 
#' @param x1mean mean of responses for patients choosing treatment 1. If study
#'   is stratified, should be vector with length equal to the
#'   number of strata.
#' @param x1var variance of responses for patients choosing treatment 1. If 
#'   study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param m1 number of patients choosing treatment 1. If study
#'   is stratified, should be vector with length equal to the
#'   number of strata.
#' @param x2mean mean of responses for patients choosing treatment 2. If study
#'   is stratified, should be vector with length equal to the
#'   number of strata.
#' @param x2var variance of responses for patients choosing treatment 2. If 
#'   study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param m2 number of patients choosing treatment 2. If study
#'   is stratified, should be vector with length equal to the
#'   number of strata.
#' @param y1mean mean of responses for patients randomized to treatment 1. If 
#'   study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param y1var variance of responses for patients randomized to treatment 1. 
#'   If study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param n1 number of patients randomized to treatment 1. If study is 
#'   stratified, should be vector with length equal to the number of 
#'   strata.
#' @param y2mean mean of responses for patients randomized to treatment 2. If 
#'   study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param y2var variance of responses for patients randomized to treatment 2. 
#'   If study is stratified, should be vector with length equal to the
#'   number of strata.
#' @param n2 number of patients randomized to treatment 2. If study is 
#'   stratified, should be vector with length equal to the number of 
#'   strata.
#' @param xi a numeric vector of the proportion of patients in each stratum. 
#'   Length of vector should equal the number of strata in the study and 
#'   sum of vector should be 1. All vector elements should be numeric 
#'   values between 0 and 1. Default is 1 (i.e. unstratified design).
#' @param nstrata number of strata. Default is 1 (i.e. unstratified design).
#' @param alpha Type I error rate, used to determine confidence interval level 
#'   for the effect estimates. Default is 0.05 (i.e. 95\% confidence interval)
#' @examples
#' x1mean <- 5
#' x1var <- 1
#' m1 <- 15
#' x2mean <- 7
#' x2var <- 1.1
#' m2 <- 35
#' y1mean <- 6
#' y1var <- 1
#' n1 <- 25
#' y2mean <- 8
#' y2var <- 1.2
#' n2 <- 25
#' fit_preference_summary(x1mean, x2var, m1, x2mean, x2var, m2, y1mean, y2var,
#'                n1, y2mean, y2var, n2)
#' @references Rucker G (1989). "A two-stage trial design for testing treatment,
#' self-selection and treatment preference effects." \emph{Stat Med}, 
#' \strong{8}(4):477-485. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/2727471}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a 
#' Stratified Doubly Randomized Preference Design." 
#' \emph{Stat Methods Med Res}. 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
fit_preference_summary <- function(x1mean, x1var, m1, x2mean, x2var, m2, y1mean,
                                   y1var, n1, y2mean, y2var, n2, xi=1, 
                                   nstrata=1, alpha=0.05) {
  
  # Compute unstratified test statistics
  unstrat_stats <- vapply(seq_len(nstrata), 
    function(i) {
      unstrat_analyze_summary_data(x1mean[i], x1var[i], m1[i], x2mean[i], 
                                   x2var[i], m2[i], y1mean[i], y1var[i], n1[i],
                                   y2mean[i], y2var[i], n2[i], alpha)
    }, data.frame(pref_effect=NA, pref_SE=NA, pref_test = NA, pref_pval = NA, 
                  pref_LB=NA, pref_UB=NA, sel_effect=NA, sel_SE=NA, 
                  sel_test = NA , sel_pval = NA, sel_LB=NA, sel_UB=NA, 
                  treat_effect=NA, treat_SE=NA, treat_test = NA, 
                  treat_pval = NA, treat_LB=NA, treat_UB=NA))
  
  #Calculate the overall effect estimate
  overall_pref_effect <- sum(
    vapply(seq_len(nstrata), 
           function(i) xi[i] * unlist(unstrat_stats[1, i]), 0.0))
  
  overall_sel_effect <- sum(
    vapply(seq_len(nstrata), 
           function(i) xi[i] * unlist(unstrat_stats[7, i]), 0.0))
  
  overall_treat_effect <- sum(
    vapply(seq_len(nstrata),
           function(i) xi[i] * unlist(unstrat_stats[13, i]), 0.0))
  
  #Calculate the overall SE
  overall_pref_SE <- sqrt(sum(
    vapply(seq_len(nstrata), 
           function(i) xi[i]^2 * unlist(unstrat_stats[2, i])^2, 0.0)))
  
  overall_sel_SE <- sqrt(sum(
    vapply(seq_len(nstrata), 
           function(i) xi[i]^2 * unlist(unstrat_stats[8, i])^2, 0.0)))
  
  overall_treat_SE <- sqrt(sum(
    vapply(seq_len(nstrata),
           function(i) xi[i]^2 * unlist(unstrat_stats[14, i])^2, 0.0)))
  
  #Calculate overall test statistic
  overall_pref_test <- overall_pref_effect/overall_pref_SE
  
  overall_sel_test <- overall_sel_effect/overall_sel_SE
  
  overall_treat_test <- overall_treat_effect/overall_treat_SE
  
  # Compute p-values (Assume test stats approximately normally distributed)
  
  # preference effect
  overall_pref_pval <- 2 * pnorm(abs(overall_pref_test), lower.tail = FALSE)
  
  # selection effect
  overall_sel_pval <- 2 * pnorm(abs(overall_sel_test), lower.tail = FALSE) 
  
  # treatment effect
  overall_treat_pval <- 2 * pnorm(abs(overall_treat_test), lower.tail = FALSE)
  
  #Compute the upper and lower bounds of the confidence interval

  zalpha <- qnorm(1-(alpha/2))
  
  overall_pref_LB <- overall_pref_effect - zalpha*overall_pref_SE
  
  overall_pref_UB <- overall_pref_effect + zalpha*overall_pref_SE
  
  overall_sel_LB <- overall_sel_effect - zalpha*overall_sel_SE
  
  overall_sel_UB <- overall_sel_effect + zalpha*overall_sel_SE
  
  overall_treat_LB <- overall_treat_effect - zalpha*overall_treat_SE
  
  overall_treat_UB <- overall_treat_effect + zalpha*overall_treat_SE

  overall_stats<-data.frame(
    overall_pref_effect = overall_pref_effect, 
    overall_pref_SE = overall_pref_SE, 
    overall_pref_test = overall_pref_test,
    overall_pref_pval = overall_pref_pval, 
    overall_pref_LB = overall_pref_LB, 
    overall_pref_UB = overall_pref_UB,
    overall_sel_effect = overall_sel_effect, 
    overall_sel_SE = overall_sel_SE, 
    overall_sel_test = overall_sel_test,
    overall_sel_pval = overall_sel_pval,  
    overall_sel_LB = overall_sel_LB, 
    overall_sel_UB = overall_sel_UB,
    overall_treat_effect = overall_treat_effect, 
    overall_treat_SE = overall_treat_SE, 
    overall_treat_test = overall_treat_test,
    overall_treat_pval = overall_treat_pval, 
    overall_treat_LB = overall_treat_LB, 
    overall_treat_UB = overall_treat_UB)
  
  list(unstratified_statistics=unstrat_stats, overall_statistics=overall_stats) 
}


#' @title Fit preference data collected from a clinical trial
#'
#' @description The variables in the formula should reference columns in the
#' data parameter and should have the following characteristics.
#' \itemize{
#'   \item{outcome: }{Numeric values giving the outcome of interest.}
#'   \item{treatment: }{Character, categorical, or integer values denoting the 
#'                      treatment received by an individual.}
#'   \item{arm: }{Character or categorical variable denoting the arm 
#'                individuals belong to. Note that these values should be 
#'                either "random" or "choice".}
#'   \item{strata: }{An optional integer value denoting which strata 
#'         individuals belong to.}
#' }
#'
#' @param form a formula of the form outcome ~ treatment:arm \{| strata\}. See
#' Details for more explanation.
#' @param data a data.frame containing variables specified in the formula
#' @importFrom stats terms
#' @examples
#' 
#' #Unstratified
#' 
#' outcome <- c(10, 8, 6, 10, 5, 8, 7, 6, 10, 12, 11, 6, 8, 10, 5, 7, 9, 12, 6,
#'              8, 9, 10, 7, 8, 11)
#' arm <- c(rep("choice", 13), rep("random", 12))
#' treatment <- c(rep(1, 5), rep(2, 8), rep(1, 6), rep(2, 6))
#' d <- data.frame(outcome=outcome, treatment=treatment, arm=arm)
#' preference(outcome ~ treatment:arm, d)
#' @export
preference <- function(form, data) {
  outcome_var <- as.character(terms(form)[[2]])
  rhs <- as.character(terms(form)[[3]])
  if (rhs[1] == "|") {
    # It's stratified.
    strat_var <- rhs[3]
    interact <- unlist(strsplit(rhs[2], ":"))
    treatment_var <- interact[1]
    arm_var <- interact[2]
  } else if (rhs[1] == ":") {
    treatment_var <- rhs[2]
    arm_var <- rhs[3]
    strat_var <- NULL
  } else {
    stop(paste0("Formula argument should be of the form outcome ~ treatment:",
                "arm {|strata}."))
  }
  if (is.null(arm_var)) {
    fit_preference(data[, outcome_var], 
                        data[, arm_var] == "random",
                        data[, treatment_var],
                        data[, arm_var])
  } else {
    fit_preference(data[, outcome_var], 
                        data[, arm_var] == "random",
                        data[, treatment_var])
  }
}

