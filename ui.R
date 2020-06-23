#ライブラリ読み込み
library(shiny)

library(DT)

#分布関数のデータ読み込み
dist <- read.csv("dist.csv")

#使用する手法
use.dist <- dist[, "distr"] %>% as.vec()
names(use.dist) <- dist[, "name"] %>% as.vec()
use.dist.sel <- use.dist[which(dist[, "use"] %>% as.vec())]



shinyUI(fluidPage(

  sidebarLayout(
    sidebarPanel(

      #ファイル選択
      fileInput("file", "Data file(.csv)",
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
                           #https://code.i-harness.com/ja/q/227d360
                           tags$head(
                             tags$link(rel = "stylesheet", type = "text/css", href = "my.css")
                           ),
                           DT::dataTableOutput("result")
                           ),
                  
                  tabPanel("Result",
                           plotOutput("result.plot"),
                           verbatimTextOutput("summary"),
                           verbatimTextOutput("gofstat")
                           ),
                  
                  tabPanel("Setting",
                           checkboxGroupInput("use", label = "Use distibution",
                                       choices = use.dist, 
                                       selected = use.dist.sel,
                                       inline = TRUE))
      )
      
      
      
      
      
      
    )
  )
))
