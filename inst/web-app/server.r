library(shiny)
library(rbokeh)
source('sample-size.r')

parse_sequence_text = function(x) {
  ret = as.numeric(x)
  if (is.na(ret)) {
    # Try to parse it as a sequence.

    to_by = regexpr("\\d+\\.?\\d*\\s?to\\s?\\d+\\.?\\d*\\s+by", x)
    to = regexpr("\\d+\\.?\\d*\\s?to\\s?\\d+\\.?\\d*\\s+", x)
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
    list(zbeta=input$zbeta, phi=phi, sigma2=sigma2, delta_pi=delta_pi,
         delta_nu=delta_nu, zalpha=input$zalpha, theta=input$theta,
         xi=xi)
  })

  get_strat_selection = reactive({ 
    params = get_inputs()
    strat_selection(params$zbeta, params$phi, params$sigma2, params$delta_pi, 
                    params$delta_nu, params$zalpha, params$theta, params$xi)
  })

  output$sample_size = renderTable({
    get_strat_selection()
  })

  output$line_graph = renderRbokeh({
    params = get_inputs()
    ss = get_strat_selection()
    figure() %>% ly_lines(lowess(cars), color = "red", width = 2)
  })

})

