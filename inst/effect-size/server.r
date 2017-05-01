library(rbokeh)
library(preference)
library(ggplot2)
library(foreach)
library(tidyverse)
library(trelliscopejs)

parse_sequence_text = function(x) {
  suppressWarnings({ret = as.numeric(x)})
  if (is.na(ret)) {
    # Try to parse it as a sequence.

    to_by = regexpr("\\d+\\.?\\d*\\s?to\\s?\\d+\\.?\\d*\\s+by", x)
    to = regexpr("\\d+\\.?\\d*\\s?to\\s?\\d+\\.?\\d*", x)
    if (to_by != -1) {
      # It is a "to-by" statement?
      vals = as.numeric(unlist(strsplit(x, "to|by")))
      ret = seq(from=vals[1], to=vals[2], by=vals[3])
    } else if (to != -1) {
      # It is a "to" statement?    
      vals = as.numeric(unlist(strsplit(x, "to")))
      ret = seq(from=vals[1], to=vals[2], by=1)
    } else {
      # We can't parse it.
      ret = NA
    }
  }
  ret
}

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {

  get_inputs = reactive({
    mu1 = parse_sequence_text(input$mu1)
    mu2 = parse_sequence_text(input$mu2)
    mu11 = parse_sequence_text(input$mu11)
    mu22 = parse_sequence_text(input$mu22)
    phi = parse_sequence_text(input$phi)
    vary_param = input$vary_param
    list(mu1=mu1, mu2=mu2, mu11=mu11, mu22=mu22, phi=phi, vary_param=vary_param)
  })

  get_effects = reactive({ 
    
    params = get_inputs()
    ret = foreach(m1=params$mu1, .combine=rbind) %:% 
      foreach(m2=params$mu2, .combine=rbind) %:%
      foreach(m11=params$mu11, .combine=rbind) %:% 
      foreach(m22=params$mu22, .combine=rbind) %:% 
      foreach(p=params$phi, .combine=rbind) %do% {
        effects = calc_effects(m1, m2, m11, m22, p)
        c(m1, m2, m11, m22, p, effects$delta_tau, effects$delta_nu, 
        effects$delta_pi)
      }
    ret = as.data.frame(ret)
    colnames(ret) = c("mu1", "mu2", "mu11", "mu22", "phi", "delta_tau",
      "delta_nu", "delta_pi")
    ret
  })

  output$sample_size = renderDataTable({
    get_effects()
  })

  output$effect_viz = renderTrelliscope({
    my_viz = function(xs, symbol, removes) {
      xs %>% gather(Effect, Value, delta_tau:delta_pi) %>% 
        ggplot(aes_string(x=symbol, y="Value")) + 
          geom_line() + facet_grid(Effect ~ .)

    }
    vary_param = get_inputs()$vary_param
    x = get_effects()
    
    group_by_names = colnames(x)[(1:5)[-which(colnames(x) == vary_param)]]
    group_by_symbols = lapply(group_by_names, as.symbol)
    xn = x %>% group_by_(.dots=group_by_symbols) %>% nest %>% 
      mutate(effect_viz= map(data, 
        function(x) my_viz(x, vary_param, group_by_names))) %>%
      trelliscope(name="Effect Sizes", nrow=1, ncol=2, panel_col="effect_viz")
  })

})

