#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # h1(fluidRow(img(src = "theStrikeZone.png", height = 150))),
    h1(fluidRow(img(src = "theStrikeZone.logo.png"))),#,img(src = "theStrikeZone.contactMe.png", height = 250))),

    # Application title
    #titlePanel(img(src = "theStrikeZone.png", height = 150)),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("pitcher", "Pitcher", 
                        choices=uniquePitcherNamesList),
            selectInput("pitchType", "Pitch Type", 
                        choices=uniquePitchTypes)
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(plotOutput("pitchLocs")),
            fluidRow(plotOutput("pitchSpeeds"),
            plotOutput("pitchBreaks"))
        )
    )
))
