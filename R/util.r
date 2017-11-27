
plot.selection_designs <- function(x, y, ...) {
}

plot.preference_designs <- function(x, y, ...) {
}

create_strat_selection_plot = function(params, df) {
  ret = NULL
  if (length(unique(df$preference_effect)) > 1 && 
      length(unique(df$selection_effect)) == 1) {
    ret = ggplot2::ggplot(data=df, 
      ggplot2::aes_string(x="preference_effect", y="sample_size")) +
      ggplot2::geom_line() + ggplot2::xlab("Preference Effect") + 
      ggplot2::ylab("Sample Size")
  }
  else if (length(unique(df$preference_effect)) == 1 &&
           length(unique(df$selection_effect)) > 1) {
    ret = ggplot2::ggplot(data=df, 
      ggplot2::aes_string(x="selection_effect", y="sample_size")) +
      ggplot2::geom_line() + ggplot2::xlab("Selection Effect") + 
      ggplot2::ylab("Sample Size")
  }
  else if (length(unique(df$preference_effect)) > 1 &&
           length(unique(df$selection_effect)) > 1) {
    df$`Preference Effect` <- factor(df$preference_effect)
    ret <- ggplot(data=df, aes(x=selection_effect, y=sample_size, 
      group=`Preference Effect`, color=`Preference Effect`)) + 
      geom_line() + xlab("Selection Effect") + ylab("Sample Size")

    # For Selection Effect
    # df$`Selection Effect` <- factor(df$selection_effect)
    # ret <- ggplot(data=df, aes(x=preference_effect, y=sample_size, 
    #   group=`Selection Effect`, color=`Selection Effect`)) + 
    #   geom_line() + xlab("Preference Effect") + ylab("Sample Size")
  }
  ret
}

create_strat_selection_df = function(params) {
  df = data.frame(
    list(delta_pi=rep(params$delta_pi, each=length(params$delta_nu)),
         delta_nu=rep(params$delta_nu, times=length(params$delta_pi))))
  ss = rep(NA, nrow(df))
  for (i in 1:nrow(df)) {
    ss[i] = 
      ceiling(selection_sample_size(params$power, params$phi, params$sigma2, 
                    df$delta_pi[i], df$delta_nu[i], params$alpha,
                    params$theta, params$xi, params$num_strata))
  }
  df$sample_size = ss
  names(df) = c("preference_effect", "selection_effect", "sample_size")
  df     
}

