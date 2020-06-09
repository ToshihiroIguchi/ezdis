#ライブラリ読み込み
library(shiny)
library(fitdistrplus)
library(dplyr)

library(DT)


#分布関数のデータ読み込み
dist <- read.csv("dist.csv")



shinyServer(function(input, output, session) {
  
  observeEvent(input$file, {
    
    #データ読み込み
    raw.data <- reactive({
      read.csv(input$file$datapath) %>%
        select_if(is.numeric) #数値のみ選択
    })
    
    #データ選択
    output$colname <- renderUI({
      selectInput("colname", label = "Use data",
                choices = c(NA, colnames(raw.data())))
      })
    

    #データ選択されたら
    observeEvent(is.null.na.null(input$colname), {

      #ヒストグラム表示
      output$gg.hist <- renderPlot({gg.hist(raw.data()[, input$colname])})
      
      #一通り分布にあてはめる
      result <- reactive({
        
        #解析するデータを準備
        data.vec <-  raw.data()[, input$colname] %>% na.omit() %>% as.vec()
        
        #戻り値をリストにする
        ret <- list()
        
        
        #プログレスバー
        withProgress(message = "Fit of univariate distributions", {
          
          #各分布関数にあてはめる
          for(i in 1:nrow(dist)){
            
            #分布関数名
            distr.name <- dist[i, "distr"] %>% as.vec()
            
            #分布関数の正式名
            distr.original.name <- dist[i, "name"] %>% as.vec()
            
            #インジケータを進ませる
            incProgress(1/nrow(dist), detail = distr.original.name)
            
            #あてはめ
            ret[[distr.name]] <- fit.dist(data = data.vec, distr = distr.name)
            
            #nameをリストに
            ret[[distr.name]]$name <- distr.original.name
            
            
          }
          
        })
        

        
        #クラスを追加
        class(ret) <- c("fit.dist", class(ret))
        
        #戻り値
        return(ret)
        
        
      })
      
      #結果一覧をdata.frameで作成
      dt.result <- reactive({summary(result())})
      
      #結果一覧表示
      output$result <- DT::renderDataTable(dt.result()[, -1], options = list(autoWidth = TRUE))
      
      #変数選択
      output$distr.sel <- renderUI({
        
        selectInput("distr.sel", label = "Distribution",
                    choices = dt.result()[, "distr"] %>% as.vec())
        
      })
      
      #結果表示
      output$result.plot <- renderPlot({plot(result()[[input$distr.sel]])})
 
      #結果のまとめ表示
      output$summary <- renderText({
        summary(result()[[input$distr.sel]]) %>% 
          capture.output() %>% paste(by = "\n") #summaryはlistなので、表示できる形に変換
        
        })
      
      #結果のまとめ表示2
      output$gofstat <- renderText({
        gofstat(result()[[input$distr.sel]]) %>% 
          capture.output() %>% paste(by = "\n")
      })
      
      
      
      
    })
    
    
    
    
    
  })



})
