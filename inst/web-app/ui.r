library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Stratified Sample-Size Calculator"),

  # Need to add tooltips for descriptions of inputs.
   
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    sliderInput("num_strata", "Number of Strata:", min=1, max=10, value=2),
    sliderInput("zbeta", "Z-Power:", min=0, max=3, value=1.282),
    textInput("phi", "Preference Rate:", value="0.5, 0.5"),
    textInput("sigma2", "Within-Stratum Variances:", value="1, 1"),
    numericInput("delta_pi", "Preference Effect:", value=1, min=0, step=0.01),
    numericInput("delta_nu", "Selection Effect:", value=1, min=0, step=0.01),
    sliderInput("zalpha", "Z Type-1 Error:", min=0, max=3, value=1.96),
    sliderInput("theta", "Theta:", min=0, max=1, value=0.5),
    textInput("xi", "Proporion of Patients in Each Stratum:", value="0.5, 0.5")
  ),
      
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("sample_size")
  )

))

