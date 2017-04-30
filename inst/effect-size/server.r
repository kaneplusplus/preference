library(rbokeh)
library(preference)
library(ggplot2)

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
    delta_pi = parse_sequence_text(input$delta_pi)
    delta_nu = parse_sequence_text(input$delta_nu)
    # Turn the text inputs into vectors of numeric values.
    phi = as.numeric(gsub(" ", "", unlist(strsplit(input$phi, ","))))
    sigma2 = as.numeric(gsub(" ", "", unlist(strsplit(input$sigma2, ","))))
    xi = as.numeric(gsub(" ", "", unlist(strsplit(input$xi, ","))))
    nr = length(delta_pi) * length(delta_nu)

    list(power=input$power, phi=phi, sigma2=sigma2, delta_pi=delta_pi,
         delta_nu=delta_nu, alpha=input$alpha, theta=input$theta,
         xi=xi, num_strata=input$num_strata)
  })

  get_strat_selection = reactive({ 
    params = get_inputs()
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
    df
  })

  output$sample_size = renderDataTable({
    x = get_strat_selection()
    names(x) = c("Preference Effect", "Selection Effect", "Sample Size")
    x
  })

  output$line_graph = renderPlot({
    params = get_inputs()
    df = get_strat_selection()
    ret = NULL
    if (length(unique(df$delta_pi)) > 1 && length(unique(df$delta_nu)) == 1) {
      ret = ggplot(data=df, aes(x=delta_pi, y=sample_size)) +
        geom_line() + xlab("Preference Effect") + ylab("Sample Size")
    }
    else if (length(unique(df$delta_pi)) == 1 && 
             length(unique(df$delta_nu)) > 1) {
      ret = ggplot(data=df, aes(x=delta_nu, y=sample_size)) +
        geom_line() + xlab("Selection Effect") + ylab("Sample Size")
    }
    else if (length(unique(df$delta_pi)) > 1 && 
             length(unique(df$delta_nu)) > 1) {
      uss = length(unique(df$sample_size))
      ret = ggplot(data=df, aes(x=delta_pi, y=delta_nu, fill=factor(sample_size))) +
        geom_tile() + xlab("Preference Effect") + ylab("Selection Effect") +
        scale_fill_manual(values=rev(heat.colors(uss)), 
          guide=guide_legend(title="Sample Size"))
        #scale_fill_continuous(guide=guide_legend(title = "Sample Size")) 
    }
    ret
  })

})

