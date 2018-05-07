
#' @title Plot the effect sizes of a preference trial
#' 
#' @description The pt_plot() function visualizes the change in the
#' preference effect, the selection effect, or both as a function of the
#' total sample size of the trial. If the preference effect varies but the
#' selection effect does not, then it plots the preference effect by the
#' total sample size. Similarly if the selection effect varies but not the
#' preference effect then selection effect vs total sample size is shown.
#' When both preference and selection effect vary then the selection effect
#' is shown conditioned on the given preference effects.
#' 
#' It is assumed that the set of trial provided as a parameter are related
#' and are comparable. For example, the function does not check to if
#' the strata are the same for all trials. If some other visualization is
#' required then the user is reminded that a preference.trial object is
#' a data frame and can be visualized in the usual way.
#'
#' @param pt an object of class preference.trial.
#' @examples
#'
#' # Plot trials with fixed power and varying preference effect.
#' trials <- pt_from_power(power = 0.8, pref_effect = seq(0.5, 2, by = 0.1), 
#'                         selection_effect = 1, treatment_effect = 1, 
#'                         sigma2 = 1, pref_prop = 0.6)
#' pt_plot(trials)
#'  
#' # Plot trials with fixed power and varying selection effect.
#' trials <- pt_from_power(power = 0.8, pref_effect = 1,
#'                         selection_effect = seq(0.5, 2, by = 0.1), 
#'                         treatment_effect = 1, sigma2 = 1, pref_prop = 0.6)
#' pt_plot(trials)
#'
#' # Plot trials with fixed power and varying preference and 
#' # selection effects.
#'
#' # the selection effects of interest
#' selection_effects <- rep(seq(0.5, 2, by = 0.1), 4)
#'
#' # the preference effects to condition on
#' pref_effects <- rep(seq(0.4, 1, by = 0.2), 
#'                     each = length(selection_effects)/4)
#'
#' trials <- pt_from_power(power = 0.8, pref_effect = pref_effects,
#'                         selection_effect = selection_effects,
#'                         treatment_effect = 1, sigma2 = 1, pref_prop = 0.6)
#' pt_plot(trials)
#'  
#' 
#' @importFrom ggplot2 ggplot aes_string geom_line ylab xlab aes
#' @importFrom tidyr gather
#' @export
pt_plot <- function(pt) {
  if (!inherits(pt, "preference.trial")) {
    stop(paste0("pt_plot doesn't know how to plot an object of type ",
                class(pt), "."))
  }
  `Preference Effect` <- `Sample Size` <- Power <- Type <- NULL
  selection_effect <- NULL
  ret <- NULL
  pt$sample_size <- apply(pt, 1, 
    function(x) {
      max(as.data.frame(x)[,c("pref_ss", "selection_ss", "treatment_ss")])
    })
  if (length(unique(pt[, "pref_effect"])) > 1 && 
      length(unique(pt[, "selection_effect"])) == 1) {
    pt$sample_size <- apply(pt, 1, 
      function(x) {
        max(as.data.frame(x)[,c("pref_ss", "selection_ss", "treatment_ss")])
      })
    ret <- ggplot(data=pt, 
      aes_string(x="pref_effect", y="sample_size")) +
      geom_line() + xlab("Preference Effect") + 
      ylab("Sample Size")
  } else if (length(unique(pt[, "pref_effect"])) == 1 &&
           length(unique(pt[, "selection_effect"])) > 1) {
    ret <- ggplot(data=pt, 
      aes_string(x="selection_effect", y="sample_size")) +
      geom_line() + xlab("Selection Effect") + 
      ylab("Selection Sample Size")
  } else if (length(unique(pt[, 'pref_effect'])) > 1 &&
           length(unique(pt[, 'selection_effect'])) > 1) {
    pt$`Preference Effect` <- factor(pt[, "pref_effect"])
    ret <- ggplot(data=pt, aes(x=selection_effect, y=sample_size, 
      group=`Preference Effect`, color=`Preference Effect`)) + 
      geom_line() + xlab("Selection Effect") + ylab("Sample Size")
  } else if ( all(pt$treatment_power == pt$selection_power) &&
              all(pt$treatment_power == pt$pref_power) &&
              all(pt$selection_power == pt$pref_power) &&
              length(unique(pt$pref_ss)) > 1 &&
              length(unique(pt$selection_ss)) > 1 &&
              length(unique(pt$treatment_ss)) > 1 ) {
    x <- pt[,c("treatment_power", "pref_ss", "selection_ss", "treatment_ss")]
    names(x) <- c("Power", "Preference", "Selection", "Treatment")
    x <- gather(x, key="Type", value=`Sample Size`, 2:4)
    ret <- ggplot(data=x, 
      aes(x=Power, y=`Sample Size`, group=Type, col=Type)) + 
        geom_line()
  } else if ( all(pt$pref_ss == pt$selection_ss) &&
              all(pt$pref_ss== pt$treatment_ss) &&
              all(pt$selection_ss== pt$treatment_ss) &&
              (length(unique(pt$pref_power)) > 1 || 
              length(unique(pt$selection_power)) > 1 ||
              length(unique(pt$treatment_power)) > 1 )) {
    x <- pt[,c("pref_ss", "pref_power", "selection_power", "treatment_power")]
    names(x) <- c("Sample Size", "Preference", "Selection", "Treatment")
    x <- gather(x, key="Type", value=`Power`, 2:4)
    ret <- ggplot(data=x, 
      aes(x=`Sample Size`, y=Power, group=Type, col=Type)) + 
        geom_line()
  } else {
    stop(paste0("Don't know how to visualize the set of trials.  ",
                "Consider building the visualization yourself."))
  }
  ret
}
