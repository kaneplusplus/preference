
# Circular indexing to recycle values.
cind <- function(i, vec) {
  (i-1) %% length(vec) + 1
}

# TODO: if theta is TRUE *and* strata = 1, use the optimal_proportion function
# to find theta. Otherwise, theta=0.5.

#setGeneric("preference.trial",
#  function(sample_size, power, sigma2, phi, delta_nu, delta_pi, delta_tau, 
#           alpha=0.05, theta=0.5, xi=1, nstrata=1L, k=1L) {
#    standardGeneric("preference.trial")
#  })

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
#' # Unstratified
#' preference.trial(pref_ss=100, pref_effect=1, selection_ss=100, 
#'   selection_effect=1, treatment_ss=100, treatment_effect=1,
#'   sigma2=1, pref_prop=0.6)
#' 
#' # Stratified
#' preference.trial(pref_ss=100, pref_effect=1, selection_ss=100, 
#'   selection_effect=1, treatment_ss=100, treatment_effect=1,
#'   sigma2=c(1, 0.8), pref_prop=c(0.6, 0.3), choice_prop=c(0.5, 0.5),
#'   stratum_prop=c(0.3, 0.7))
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
  if (length(stratum_prop) != length(choice_prop) ||
      length(choice_prop) != length(pref_prop)) {
    stop(paste("The stratum_prop, choice_prop, and pref_pop parameters",
               "must have the same length."))
  }
  if (sum(stratum_prop) != 1) {
    stop("Elements of the stratum_prop parameter must sum to 1.")
  }
  if (any(!is.numeric(stratum_prop)) || any(stratum_prop < 0) ||
      any(stratum_prop > 1)) {
    stop("The stratum_prop parameter must be a numeric value in [0, 1]")
  }
  if (any(!is.numeric(choice_prop)) || any(choice_prop < 0) ||
      any(choice_prop > 1)) {
    stop("The choice_prop parameter must be a numeric value in [0, 1]")
  }
  if (any(!is.numeric(pref_prop)) || any(pref_prop < 0) ||
      any(pref_prop > 1)) {
    stop("The pref_prop parameter must be a numeric value in [0, 1]")
  }
  if (!is.numeric(sigma2) || length(sigma2) != length(stratum_prop) ||
      any(sigma2 < 0)) {
    stop(paste("The sigma2 parameter must be a numeric vector where each",
               "value is greater than zero and the length is equal to",
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

#setMethod("preference.trial",
#  signature(sample_size="integer", power="missing", sigma2="numeric", 
#            phi="numeric", delta_nu="numeric", delta_pi="numeric", 
#            delta_tau="numeric", alpha="numeric", theta="numeric",
#            xi="numeric", nstrata="integer", k="integer"),
#  function(sample_size, sigma2, phi, delta_nu, delta_pi, delta_tau, 
#           alpha, theta, xi, nstrata, k) {
#    ret <- NULL
#    max_len <- max(c(length(sample_size), length(sigma2), length(phi), 
#                     length(delta_nu), length(delta_pi), length(delta_tau), 
#                     length(alpha), length(theta), length(xi), length(nstrata), 
#                     length(k)))
#    for (i in 1:length(max_len)) {
#      ret <- rbind(ret, 
#        data.frame(sample_size=sample_size[cind(i, sample_size)],
#                   sigma2=sigma2[cind(i, sigma2)], 
#                   phi=phi[cind(i, phi)],
#                   delta_nu=delta_nu[cind(i, delta_nu)],
#                   delta_pi=delta_pi[cind(i, delta_pi)],
#                   delta_tau=delta_tau[cind(i, delta_tau)],
#                   alpha=alpha[cind(i, alpha)], 
#                   theta=theta[cind(i, theta)],
#                   xi=xi[cind(i, xi)],
#                   nstrata=nstrata[cind(i, nstrata)],
#                   k=k[cind(i, k)]))
#    }
#    class(ret) <- c("preference.trial", class(ret))
#    ret
#  })

#setMethod("preference.trial",
#  signature(sample_size="missing", power="numeric", sigma2="numeric", 
#            phi="numeric", delta_nu="numeric", delta_pi="numeric", 
#            delta_tau="numeric", alpha="numeric", theta="numeric", 
#            xi="numeric", nstrata="integer", k="integer"),
#  function(power, sigma2, phi, delta_nu, delta_pi, delta_tau, 
#           alpha, theta, xi, nstrata, k) {
#    # We'll set sample size to power to create the ret variable. Then
#    # we'll fill it in with the actual sample size.
#    sample_size <- power
#    ret <- NULL
#    max_len <- max(length(power), length(sigma2), length(phi), 
#                   length(delta_nu), length(delta_pi), length(delta_tau), 
#                   length(alpha), length(theta), length(xi), 
#                   length(nstrata), length(k))
#    for (i in 1:length(max_len)) {
#      ret <- rbind(ret, 
#        data.frame(sample_size=power[cind(i, power)], 
#                   sigma2=sigma2[cind(i, sigma2)], 
#                   phi=phi[cind(i, phi)], 
#                   delta_nu=delta_nu[cind(i, delta_nu)],
#                   delta_pi=delta_pi[cind(i, delta_pi)],
#                   delta_tau=delta_tau[cind(i, delta_tau)], 
#                   alpha=alpha[cind(i, alpha)], 
#                   theta=theta[cind(i, theta)], 
#                   xi=xi[cind(i, xi)],
#                   nstrata=nstrata[cind(i, nstrata)],
#                   k=k[cind(i, k)]))
#    }
#    for (i in 1:nrow(ret)) {
#      ret$sample_size[i] <- max(unlist(
#        overall_sample_size(power[cind(i, seq_len(max_len))], ret$phi[i], 
#                            ret$sigma2[i], ret$delta_pi[i], ret$delta_nu[i], 
#                            ret$delta_tau[i], ret$alpha[i], ret$theta[i], 
#                            ret$xi[i], ret$nstrata[i])))
#    }
#    class(ret) <- c("preference.trial", class(ret))
#    ret
#  })
