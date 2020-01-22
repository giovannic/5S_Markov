#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#install.packages('shiny', 'MCMCpack', 'ggplot2')
library(shiny)
source('./Markov_PSA.R')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Price Sensitivity Analysis - Adherence Intervention"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("int_costs_yearly",
                        "Yearly cost of the adherence intervention",
                        min = 100,
                        max = 100000,
                        value = 5500)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plane")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$plane <- renderPlot({
        plot_plane(input$int_costs_yearly)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
