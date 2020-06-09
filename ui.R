#ライブラリ読み込み
library(shiny)

library(DT)



shinyUI(fluidPage(

  sidebarLayout(
    sidebarPanel(

      #ファイル選択
      fileInput("file", "Choose CSV File",
                accept = c("csv")
                ),
      
      #変数の選択
      htmlOutput("colname"),
      
      
      #分布の選択
      htmlOutput("distr.sel")

    ),
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Histgram", 
                           plotOutput("gg.hist")
                           ),
                  
                  tabPanel("Table", 
                           DT::dataTableOutput("result")
                           ),
                  
                  tabPanel("Result",
                           plotOutput("result.plot"),
                           verbatimTextOutput("summary"),
                           verbatimTextOutput("gofstat")
                           )
      )
      
      
      
      
      
      
    )
  )
))
