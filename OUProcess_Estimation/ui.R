
# SOFTWAREPROJEKT OU-Prozesse 
# Josef Starkmann
# dX_t = \kappa(\alpha-X_t)dt + \sigma dB(t)
# X_0 = a
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

# IMPORT SHINY LIBRARY
library(shiny)
X=# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  # TITLE
  withMathJax(titlePanel("Computation of $$dX_t = \\kappa(\\alpha-X_t)dt + \\sigma dB(t)$$ $$ X_0 = a$$")),
  # SHORT HELPTEXT WITH GENERAL INFORMATION
  helpText("Softwareprojekt OU-Process - Josef Starkmann"),
  
  # CREATE A SIDEBARLAYOUT WITH SIDEBAR PANEL
  sidebarLayout(
    sidebarPanel(
      ###### Parameter ###############
      # ANZAHL Realisierungen pro Datensatz n
      
      # ANZAHL DATENSAETZE  D
      
      # X_0 = a
      
      # \kappa
      
      # \alpha
      
      # \sigma
      
      # \delta, dt
      ##############################
       # CREATE ALL IMPORTANT INPUT FIELDS
       numericInput("D",h4('number of time series\' D'),value = 1000, step =1,width = '50%'),#'
       numericInput("n",h4('steps per time series n'),value = 1000, step =1,width = '50%'),
       numericInput("a",h4(withMathJax("starting value \\( a \\)")),value = 0, step =0.01,width = '50%'),
       numericInput("delta",h4(withMathJax("time lag \\( \\delta \\)")),value = 1e-2, step =1e-5,width = '50%'),
       numericInput("kappa",h4(withMathJax("\\( \\kappa \\)")),value = 1, step =0.01,width = '50%'),
       numericInput("alpha",h4(withMathJax("\\( \\alpha \\)")),value = 1, step =0.01,width = '50%'),
       numericInput("sigma",h4(withMathJax("\\( \\sigma \\)")),value = 1, step =0.01,width = '50%'),
       numericInput("seed",div('seed',style = "color:green"),value = 3141593,step = 1, width = '50%'),
       actionButton("create_data_button","Create Data and Compute Estimates")
       
    ),
    
    # CREATE MAINPANEL WITH PLOTS AND INFORMATION
    mainPanel(
       #numericInput with Number of the Plot to show(can be one of the D Datasets)
       numericInput("plt",h4('index of shown plot'),value = 10,min = 1,max = 20, step =1,width = '50%'),
       
       # Plot of the choosen process, can be changed dynamically
       plotOutput("process_plot"),
       
       # Plot of the estimatetes of alpha, kappa, sigma +mean_error_depending on delta and variance depending on delta
       plotOutput("distPlot"),
       
       # FLUID ROW WITH TWO COLUMNS OF WIDTH 6 EACH ONE ROW HAS WIDTH 12
       fluidRow(
         # FIRST LEFT ROW WITH THREE ELEMENTS A SLDER FOR DELTA, a correlating numeric input with the same value and a button 
         column(6,
                sliderInput("delta_new_slide",h4(withMathJax("Increase \\( \\delta \\) with factor")),min = 1, max = 100, value =1),
                numericInput("delta_new_in", label = "",value = 1),
                uiOutput('text_delta'),
                actionButton("recompute","Recompute Estimates")
                ),
         # SECOND OUTPUT IN THE COLUMN WITH THE GENERAL INFOS     
         column(6,uiOutput('infos'),uiOutput('infos2'))
          
             
       ),
       # FURTHER PLOTS THAT OCCUR AFTER WE CLICK RECOMPUTE ESTIMATES
       plotOutput("further_plots")
     )
  )
))
