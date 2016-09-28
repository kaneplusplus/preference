library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("Old Faithful Geyser Data"),

   # Need to add tooltips for descriptions of inputs.
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("num_strata",
                     "Number of Strata:",
                     min=1,
                     max=10,
                     value=2),
         sliderInput("zbeta",
                     "Z-Power:",
                     min=0,
                     max=3,
                     value=1.282),
         # slider input, 1 per strata, between 0 and 1.
         # numeric input 1 per strata, greater than 0.
         # single value, numeric
         # single value, numeric
         sliderInput("zalpha",
                     "Z Type-1 Error:",
                     min=0,
                     max=3,
                     value=1.96),
         sliderInput("theta",
                     "Theta:",
                     min=0,
                     max=1,
                     value=0.5),
         # slider input, 1 per strata, must sum to 1. Text box for now.
         # numeric input for treatment effect.
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
})

# Run the application 
shinyApp(ui = ui, server = server)


