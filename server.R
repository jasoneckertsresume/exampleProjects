#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

setwd("C:\\Users\\jason\\Desktop\\baseBall")

# Data Import:
# Original Data Location: https://www.kaggle.com/pschale/mlb-pitch-data-20152018#pitches.csv
#
# atBats <- read.csv("atbats.csv")
# ejections <- read.csv("ejections.csv")
# games <- read.csv("games.csv")
# pitches <- read_csv("pitches.csv")
# pitches2 <- read_csv("pitches.csv", col_types = cols(.default = col_character()))
# playerNames <- read.csv("player_names.csv")

# Used to correct an issue with ab_id that was causing truncation and loss of specificity due to
# ab_id being larger than num can handle
#
# a <- pitches2$ab_id
# a <- substr(a, 0, 10)
# pitches$ab_id <- strtoi(a)
# 
# temp <- data.frame(atBats$ab_id, atBats$batter_id, atBats$pitcher_id)
# colnames(temp)<- c("ab_id","batter_id", "pitcher_id")
# pitchesTemp <- left_join(pitches, temp)
# pitches <- pitchesTemp

# Creation of sample dataset 
#
pitchesSample <- pitches [1:1000000, ]
write.csv(pitchesSample, file = "pitchesSample.csv")

# Added for utility so that a sample could be used for testing in-order to reduce latency and the original imports
# do not need to repeated
#
# pitchesSample <- read.csv("pitchesSample.csv")
# pitchesSample1 <- left_join(pitchesSample, playerNames, by = c("pitcher_id" = "id"))
# pitchesSample1$fullName <- paste(pitchesSample1$first_name, pitchesSample1$last_name)
# # 
# uniqueNames <- data.frame(pitchesSample1$first_name, pitchesSample1$last_name)
# colnames(uniqueNames) <- c("first_name", "last_name")
# uniqueNames <- unique(uniqueNames)
# uniquePitcherNamesList <- paste(uniqueNames$first_name, uniqueNames$last_name)
# uniquePitcherNamesList <- sort(uniquePitcherNamesList)
# uniquePitcherNamesList <- prepend(uniquePitcherNamesList, "All")
# # 
# uniquePitchTypes <- unique(pitchesSample1$pitch_type)
# uniquePitchTypes <- prepend(as.character(uniquePitchTypes), "All")

# Creation of the shiny server
shinyServer(function(input, output) {
    output$pitchLocs <- renderPlot({
        
        hmTab <- data.frame(pitchesSample1$px, pitchesSample1$pz, pitchesSample1$fullName, pitchesSample1$pitch_type)
        colnames(hmTab) <- c("px", "pz", "fullName", "pitchType")
        if(input$pitcher != "All")
        { 
            hmTab <- subset(hmTab, fullName == input$pitcher)
        }
        if(input$pitchType != "All")
        { 
            hmTab <- subset(hmTab, pitchType == input$pitchType)
        }
        ?subset
        colnames(hmTab) <- c("px", "pz")
        #hmTab <- as.matrix(hmTab)
         
        ggplot(hmTab, aes(hmTab$px, hmTab$pz)) +
            geom_raster(aes(fill = density))
        
        locs <- ggplot(hmTab, aes(px, pz))
        locs + geom_tile()
        locs + stat_bin2d(aes(fill = stat(density)), binwidth = c(.2,.2)) + scale_fill_gradientn(colours = heat.colors(15))
        #locs + stat_bin2d(aes(fill = stat(density)), binwidth = c(.3333,.3333)) + scale_fill_gradientn(colours = heat.colors(15))
    })
    output$pitchSpeeds <- renderPlot({
        
        hmTab1 <- data.frame(pitchesSample1$start_speed, pitchesSample1$fullName, pitchesSample1$pitch_type)
        colnames(hmTab1) <- c("speed", "fullName", "pitchType")
        if(input$pitcher != "All")
        { 
            hmTab1 <- subset(hmTab1, fullName == input$pitcher)
        }
        if(input$pitchType != "All")
        { 
            hmTab1 <- subset(hmTab1, pitchType == input$pitchType)
        }
        #hmTab <- as.matrix(hmTab)
        
        
        ggplot(data = hmTab1, aes(hmTab1$speed)) + 
            geom_histogram(binwidth = 3)
    })
    output$pitchBreaks <- renderPlot({
        
        hmTab1 <- data.frame(pitchesSample1$break_y, pitchesSample1$fullName, pitchesSample1$pitch_type)
        colnames(hmTab1) <- c("breaks", "fullName", "pitchType")
        if(input$pitcher != "All")
        { 
            hmTab1 <- subset(hmTab1, fullName == input$pitcher)
        }
        if(input$pitchType != "All")
        { 
            hmTab1 <- subset(hmTab1, pitchType == input$pitchType)
        }
        #hmTab <- as.matrix(hmTab)
        
        
        ggplot(data = hmTab1, aes(hmTab1$breaks)) + 
            geom_histogram(binwidth = .1)
    })
})
