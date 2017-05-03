
create_strat_selection_plot = function(params, df) {
  ret = NULL
  if (length(unique(df$preference_effect)) > 1 && 
      length(unique(df$selection_effect)) == 1) {
    ret = ggplot(data=df, aes(x=preference_effect, y=sample_size)) +
      geom_line() + xlab("Preference Effect") + ylab("Sample Size")
  }
  else if (length(unique(df$preference_effect)) == 1 &&
           length(unique(df$selection_effect)) > 1) {
    ret = ggplot(data=df, aes(x=selection_effect, y=sample_size)) +
      geom_line() + xlab("Selection Effect") + ylab("Sample Size")
  }
  else if (length(unique(df$preference_effect)) > 1 &&
           length(unique(df$selection_effect)) > 1) {
    uss = length(unique(df$sample_size))
    ret=ggplot(data=df, aes(x=preference_effect, y=selection_effect, 
      fill=factor(sample_size)))+
      geom_tile() + xlab("Preference Effect") + ylab("Selection Effect") +
      scale_fill_manual(values=rev(heat.colors(uss)),
        guide=guide_legend(title="Sample Size"))
      #scale_fill_continuous(guide=guide_legend(title = "Sample Size")) 
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
      ceiling(n_sel(params$power, params$phi, params$sigma2, 
                              df$delta_pi[i], df$delta_nu[i], params$alpha,
                              params$theta, params$xi, params$num_strata))
  }
  df$sample_size = ss
  names(df) = c("preference_effect", "selection_effect", "sample_size")
  df     
}

