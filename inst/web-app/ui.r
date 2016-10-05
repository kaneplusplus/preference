library(shiny)
library(shinyBS)

# Define UI for application that draws a histogram
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Stratified Sample-Size Calculator"),

  # Need to add tooltips for descriptions of inputs.

  # Pass syntax: 1 to by 0.4
   
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    sliderInput("num_strata", "Number of Strata:", min=1, max=10, value=2),
    sliderInput("zbeta", "Z-Power:", min=0, max=3, value=1.282),
    textInput("phi", "Preference Rate:", value="0.5, 0.5"),
    textInput("sigma2", "Within-Stratum Variances:", value="1, 1"),
    textInput("delta_pi", "Preference Effect:", value="1"),
    textInput("delta_nu", "Selection Effect:", value="1"),
    sliderInput("zalpha", "Z Type-1 Error:", min=0, max=3, value=1.96),
    sliderInput("theta", "Theta:", min=0, max=1, value=0.5),
    textInput("xi", "Proporion of Patients in Each Stratum:", value="0.5, 0.5"),
    bsTooltip("num_strata", "PUT TEXT HERE.", "right", 
              options(container="list")),
    bsTooltip("zbeta", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("phi", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("sigma2", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("delta_pi", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("delta_nu", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("zalpha", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("theta", "PUT TEXT HERE.", "right", options(container="list")),
    bsTooltip("xi", "PUT TEXT HERE.", "right", options(container="list"))
    # It would be nice to get rid of tooltips when the user doesn't want
    # them.
    #checkboxInput("tooltip_checkbox", label="Include help tooltips", 
    #              value=FALSE),
  ),
      
  # Show a plot of the generated distribution
  mainPanel(
    tableOutput("sample_size"),
    uiOutput("update_tooltip")
  )

))

