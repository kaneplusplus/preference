
# Internal function for creating a single preference trial object.
preference.trial.single <- function(pref_ss, pref_effect, selection_ss, 
  selection_effect, treatment_ss, treatment_effect, sigma2, 
  pref_prop, choice_prop=0.05, stratum_prop=1, alpha=0.05) {
  if (!is.numeric(pref_ss) || pref_ss < 1) {
    stop("The pref_ss parameter should be a numeric value of at least 1.")
  }
  if (!is.numeric(pref_effect) || length(pref_effect) != 1) {
    stop("The pref_effect parameter should be a single numeric value.")
  }
  if (!is.numeric(selection_ss) || selection_ss < 1) {
    stop("The selection_ss parameter should be a numeric value of at least 1.")
  }
  if (!is.numeric(selection_effect) || length(selection_effect) != 1) {
    stop("The selection_effect parameter should be a single numeric value.")
  }
  if (!is.numeric(treatment_ss) || treatment_ss < 1) {
    stop("The treatment_ss parameter should be a numeric value of at least 1.")
  }
  if (!is.numeric(treatment_effect) || length(treatment_effect) != 1) {
    stop("The treatment_effect parameter should be a single numeric value.")
  }
  if (!is.list(pref_prop)) pref_prop <- list(pref_prop)
  if (!is.list(choice_prop)) choice_prop <- list(choice_prop)
  if (!is.list(stratum_prop)) stratum_prop <- list(stratum_prop)
  if (!is.list(sigma2)) sigma2<- list(sigma2)
   
  if (length(stratum_prop[[1]]) != length(choice_prop[[1]]) ||
      length(choice_prop[[1]]) != length(pref_prop[[1]])) {
    stop(paste("The stratum_prop, choice_prop, and pref_pop parameters",
               "must have the same length."))
  }
  if (sum(stratum_prop[[1]]) != 1) {
    stop("Elements of the stratum_prop parameter must sum to 1.")
  }
  if (any(!is.numeric(stratum_prop[[1]])) || any(stratum_prop[[1]] < 0) ||
      any(stratum_prop[[1]] > 1)) {
    stop("The stratum_prop parameter must be a list numeric values in [0, 1]")
  }
  if (any(!is.numeric(choice_prop[[1]])) || any(choice_prop[[1]] < 0) ||
      any(choice_prop[[1]] > 1)) {
    stop("The choice_prop parameter must be a list of numeric values in [0, 1]")
  }
  if (any(!is.numeric(pref_prop[[1]])) || any(pref_prop[[1]] < 0) ||
      any(pref_prop[[1]] > 1)) {
    stop("The pref_prop parameter must be a list of numeric values in [0, 1]")
  }
  if (!is.numeric(sigma2[[1]]) || 
      length(sigma2[[1]]) != length(stratum_prop[[1]]) ||
      any(sigma2[[1]] < 0)) {
    stop(paste("The sigma2 parameter must be a list of numeric values",
               "greater than zero and the length is equal to",
               "the stratum_prop, choice_prop, and pref_prop parameters."))
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || any(alpha < 0 | alpha > 1)) {
    stop("The alpha parameter must be a single numeric value in [0, 1].")
  }
  ret <- data.frame(pref_ss=pref_ss, pref_effect=pref_effect, 
    selection_ss=selection_ss, selection_effect=selection_effect, 
    treatment_ss=treatment_ss, treatment_effect=treatment_effect, 
    alpha=alpha)
  ret$pref_prop=list(pref_prop)
  ret$choice_prop=list(choice_prop)
  ret$stratum_prop=list(stratum_prop)
  ret$sigma2=list(sigma2)
  class(ret) <- c("preference.trial", class(ret))
  ret
}

# Circular indexing to recycle values.
cind <- function(i, vec_len) {
  (i-1) %% vec_len + 1
}

#' Create a Preference Trial
#' 
#' @param pref_ss the sample size of the preference arm.
#' @param pref_effect the effect size of the preference arm (delta_pi). 
#' @param selection_ss the sample size of the selection arm.
#' @param selection_effect the effect size of selection arm (delta_nu).
#' @param treatment_ss the sample size of the treatment arm .
#' @param treatment_effect the sample size of the treatment arm (delta_tau)
#' @param sigma2 the variance estimate of ???. This value should be 
#' positive numeric values. If study is stratified, should be vector of 
#' within-stratum variances with length equal to the number of strata in the
#' study.
#' @param pref_prop the proportion of patients preferring treatment 1. This
#' value should be between 0 and 1 (phi).
#' @param choice_prop the proportion of patients assigned to choice arm in 
#' the initial randomization. Should be numeric value between
#' 0 and 1 (default=0.5) (theta).
#' @param stratum_prop xi a numeric vector of the proportion of patients in 
#' each stratum. Length of vector should equal the number of strata in the 
#' study and sum of vector should be 1. All vector elements should be numeric
#' values between 0 and 1. Default is 1 (i.e. unstratified design) (xi).
#' @param alpha the desired type I error rate.
#' @examples
#'
#' # Unstratified single trial.
#' preference.trial(pref_ss=100, pref_effect=1, selection_ss=100, 
#'   selection_effect=1, treatment_ss=100, treatment_effect=1,
#'   sigma2=1, pref_prop=0.6)
#' 
#' # Stratified single trial.
#' preference.trial(pref_ss=100, pref_effect=1, selection_ss=100,
#'   selection_effect=1, treatment_ss=100, treatment_effect=1,
#'   sigma2=list(c(1, 0.8)), pref_prop=list(c(0.6, 0.3)),
#'   choice_prop=list(c(0.5, 0.5)), stratum_prop=list(c(0.3, 0.7)))
#' 
#' @references Turner RM, et al. (2014). "Sample Size and Power When Designing
#'  a Randomized Trial for the Estimation of Treatment, Selection, and
#'  Preference Effects." \emph{Medical Decision Making}, \strong{34}:711-719.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24695962}{PubMed})
#' @references Cameron B, Esserman D (2016). "Sample Size and Power for a
#' Stratified Doubly Randomized Preference Design." \emph{Stat Methods Med Res}.
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27872194}{PubMed})
#' @export
preference.trial <- function(pref_ss, pref_effect, selection_ss, 
  selection_effect, treatment_ss, treatment_effect, sigma2, 
  pref_prop, choice_prop=0.5, stratum_prop=1, alpha=0.05) {

  # Evaluate the arguments once from the match.call return.
  args <- lapply(as.list(match.call())[-1], eval)
  arg_lens <- vapply(args, length, 0L)

  if (all(arg_lens == 1)) {
    preference.trial.single(pref_ss, pref_effect, selection_ss,
      selection_effect, treatment_ss, treatment_effect, sigma2,
      pref_prop, choice_prop, stratum_prop, alpha)
  } else {
    # Use circular to create multiple trials.
    max_arg_len <- max(arg_lens)
    if (any(max_arg_len %% arg_lens != 0)) {
      warning(paste("One argument is not a sub-multiple or multiple",
                    "of the longest argument."))
    }

    # Create a set of preference trials.
    exp_arg_df <- as.data.frame(
      Map(function(x) x[cind(1:max_arg_len, length(x))], args))
    expanded_args <- 
      Map(function(x) as.list(exp_arg_df[x,]), 1:nrow(exp_arg_df))
    Reduce(rbind,
      Map(function(x) do.call(preference.trial.single, x), expanded_args))
  }
}

#' Design Preference Trials with Power Constraint(s)
#' @examples
#' pt_from_power(power=seq(.1, 0.8, by=0.1), pref_effect=1, selection_effect=1, 
#'   treatment_effect=1, sigma2=1, pref_prop=0.6)
#' @export
pt_from_power <- function(power, pref_effect, selection_effect, 
  treatment_effect, sigma2, pref_prop, choice_prop=0.5, stratum_prop,
  alpha=0.05) {

  # Check the power parameter. Other parameters will be checked later.
  if(!is.numeric(power) || power <= 0 || power >= 1) {
    stop('Power must be numeric in [0,1]')
  }


  # Use preference.trial to create the data frame. Use power as a
  # place holder. Fill it in after the object is created.
  args <- lapply(as.list(match.call())[-1], eval)
  args$pref_ss <- args$selection_ss <- args$treatment_ss <- rep(1,length(power))
  args <- args[names(args) != "power"]
  ret <- do.call(preference.trial, args)
  for (i in 1:nrow(ret)) {
    sss <- overall_sample_size(power[cind(i, length(power))], 
                               ret$pref_prop[[i]],
                               ret$sigma2[[i]], ret$pref_effect[i], 
                               ret$selection_effect[i], ret$treatment_effect[i],
                               ret$alpha[i], ret$choice_prop[[i]], 
                               ret$stratum_prop[[i]], 
                               length(ret$stratum_prop[[i]]))
    ret$treatment_ss[i] <- sss$treatment[1]
    ret$pref_ss[i] <- sss$preference[1]
    ret$selection_ss[i] <- sss$selection[1]
  }
  ret
}

