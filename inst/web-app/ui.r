library(shiny)
library(shinyBS)
library(rbokeh)

# Define UI for application that draws a histogram
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Stratified Sample-Size Calculator"),

  # Need to add tooltips for descriptions of inputs.

  # Pass syntax: 1 to by 0.4
   
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    sliderInput("num_strata", "Number of Strata:", min=1, max=10, value=2),
    sliderInput("beta", "Power:", min=0, max=3, value=1.282),
    textInput("phi", "Preference Rate:", value="0.5, 0.5"),
    textInput("sigma2", "Within-Stratum Variances:", value="1, 1"),
    textInput("delta_pi", "Preference Effect:", value="1"),
    textInput("delta_nu", "Selection Effect:", value="1"),
    sliderInput("alpha", "Type-1 Error Rate:", min=0, max=3, value=1.96),
    sliderInput("theta", "Choice-arm Proportion:", min=0, max=1, value=0.5),
    textInput("xi", "Proporion of Patients in Each Stratum:", value="0.5, 0.5"),
#    bsTooltip("num_strata", "The number of Strata", "top", 
#              options(container="list")),
#    bsTooltip("beta", "The desired study power", "top", options(container="list")),
    bsTooltip("phi", "The proportion of patients preferring treatment 1 within each stratum. There should be one proportion per stratum.", "top", options(container="list")),
    bsTooltip("sigma2", "The within-stratum variances. There should be one variance value per stratum.", "top", options(container="list")),
    bsTooltip("delta_pi", "The overall study preference effect.", "top", options(container="list")),
    bsTooltip("delta_nu", "The overall study selection effect.", "top", options(container="list")),
#    bsTooltip("alpha", "The desired type I error rate.", "top", options(container="list")),
#    bsTooltip("theta", "The proportion of patients assigned to the choice arm in the initial randomization.", "top", options(container="list")),
    bsTooltip("xi", "The proportion of patients in each stratum. There should be one proportion-value per stratum.", "top", options(container="list"))
    # It would be nice to get rid of tooltips when the user doesn't want
    # them.
    #checkboxInput("tooltip_checkbox", label="Include help tooltips", 
    #              value=FALSE),
  ),
      
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("sample_size"),
    rbokehOutput("line_graph")
  )

))

