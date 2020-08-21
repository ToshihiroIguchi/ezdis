#ライブラリ読み込み
library(shiny)
library(DT)

#分布関数のデータ読み込み
dist <- read_csv("dist.csv")


#使用する手法
use.dist <- dist[, "distr"] %>% as.vec()
names(use.dist) <- dist[, "name"] %>% as.vec()
use.dist.sel <- use.dist[which(dist[, "use"] %>% as.vec())]

#A character string coding for the fitting method #フィッティング法
fitdist.method <- c("mle", "mme", "qme", "mge", "mse")
names(fitdist.method) <- c("maximum likelihood estimation", 
                           "moment matching estimation",
                           "quantile matching estimation",
                           "maximum goodness-of-fit estimation",
                           "maximum spacing estimation")


shinyUI(fluidPage(

  sidebarLayout(
    sidebarPanel(

      #ファイル選択
      fileInput("file", "Data file(.csv, .xls, .xlsx, .xlsm)",
                accept = c("csv", "xls", "xlsx", "xlsm")
                ),
      
      #変数の選択
      htmlOutput("colname"),
      
      
      #分布の選択
      htmlOutput("distr.sel")

    ),
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  
                  tabPanel("Histgram", 
                           plotOutput("gg.hist"),
                           verbatimTextOutput("vec.summary")
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
                           verbatimTextOutput("summary")
                           ),
                  
                  tabPanel("Paper",
                           fluidRow(
                             column(6, htmlOutput("paper.method")),
                             column(6, htmlOutput("paper.rank"))
                           ),
                          plotOutput("plot_paper"),
                          verbatimTextOutput("summary2")
                          ),
                  
                  tabPanel("Calc",
                           fluidRow(
                             column(6, 
                                    htmlOutput("input.q"),
                                    h4(textOutput("output.p"))
                                    ),
                             
                             column(6, 
                                    htmlOutput("input.p"), 
                                    h4(htmlOutput("output.q"))
                                    )
                           )
                           
                           
                           ),
                  
                  
                  tabPanel("Setting",
                           checkboxGroupInput("use", label = "Use distibution",
                                       choices = use.dist, 
                                       selected = use.dist.sel,
                                       inline = TRUE),
                           
                           selectInput("fitdist.method", label = "Fitting method",
                                       choices = fitdist.method),
                           
                           numericInput("timeout", label = "Seconds to timeout", 
                                        value = 30, min = 1, step = 1)
                           
                           ),
                  tabPanel("Debug",
                           verbatimTextOutput("fit.dist.res")
                           )
      )

    )
  )
))
