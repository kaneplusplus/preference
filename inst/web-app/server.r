library(shiny)
source('sample-size.r')

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  output$sample_size = renderTable({
    # Turn the text inputs into vectors of numeric values.
    phi = as.numeric(gsub(" ", "", unlist(strsplit(input$phi, ","))))
    sigma2 = as.numeric(gsub(" ", "", unlist(strsplit(input$sigma2, ","))))
    xi = as.numeric(gsub(" ", "", unlist(strsplit(input$xi, ","))))
    strat_selection(input$zbeta, phi, sigma2, input$delta_pi, input$delta_nu,
                    input$zalpha, input$theta, xi)
  })
})

