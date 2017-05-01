library(shiny)
library(shinyBS)
library(ggplot2)
library(trelliscopejs)

# Define UI for application that draws a histogram
ui <- shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Stratified Selection Sample-Size Calculator"),

  # Need to add tooltips for descriptions of inputs.

  # Pass syntax: 1 to by 0.4
   
  # Sidebar with a slider input for number of bins 
  sidebarPanel(
    textInput("mu1", "Random Arm Treatment 1 Mean Response:", 
      value="2 to 10"),
    textInput("mu2", "Random Arm Treatment 2 Mean Response:", 
      value="2"),
    textInput("mu11", "Preference Arm Treatment 1 Mean Response:", 
      value="2"),
    textInput("mu22", "Preference Arm Treatment 2 Mean Response:", 
      value="2"),
    textInput("phi", "Proportion of Patients Preferring Treatment 1:",
      value="0.4 to 0.6 by 0.1"),
    textInput("vary_param", "Parameter to vary over", value="mu1")
  ),
      
  # Show a plot of the generated distribution
  mainPanel(
    trelliscopeOutput("effect_viz"),
    dataTableOutput("effect_size")
  )

))

