#ライブラリ読み込み
library(shiny)
library(fitdistrplus)
library(dplyr)

library(DT)

#分布関数のデータ読み込み
dist <- read.csv("dist.csv")

#Maximum upload size exceededを回避
#100MB設定
#https://github.com/rstudio/shiny-examples/blob/master/066-upload-file/server.R
#https://stackoverflow.com/questions/18037737/how-to-change-maximum-upload-size-exceeded-restriction-in-shiny-and-save-user
options(shiny.maxRequestSize = 100*1024^2)

#確率紙
paper.method <- c("norm", "lnorm", "gumbel", "weibull", "exp")
names(paper.method) <- c("Normal probability plot",
                         "Log normal probability plot",
                         "Gumbel probability plot",
                         "Weibull probability plot",
                         "Exponential probability paper")
paper.method.pos <- paper.method[-c(2, 4)]

#ランクの計算方法
paper.rank <- c("median", "mean")
names(paper.rank) <- c("Median rank", "Mean rank")


shinyServer(function(input, output, session) {
  
  observeEvent(input$file, {
    
    #データ読み込み
    raw.data <- reactive({
      try(
        read.data(input$file$datapath) %>%
          select_if(is.numeric) #数値のみ選択
        , silent = TRUE)
      })
    
    #データ選択
    output$colname <- renderUI({
      if(class(raw.data())[1] != "try-error"){
        selectInput("colname", label = "Use data", choices = c(NA, colnames(raw.data())))
      }else{
        NULL
      }
    })
    
    #ファイルの内容がおかしかったらエラーメッセージ
    if(class(raw.data())[1] == "try-error"){
      showModal(modalDialog(
        title = "Error",
        "Failed to read the numerical data from the file. Check the contents of the file.",
        easyClose = TRUE
      ))
    }
    
    #データ選択されたら
    observeEvent(is.null.na.null(input$colname), {
      
      #データを作る
      vec.data <- reactive({
        ret <- try(raw.data()[, input$colname] %>% na.omit() %>% as.vec(), silent = TRUE)
        if(class(ret)[1] == "try-error"){ret <- NULL}
        return(ret)
      })
      
      #どの確率紙を使うか
      output$paper.method <- renderUI({
        
        if(is.null(vec.data() )){return(NULL)}
        
        if(min(vec.data()) > 0){
          selectInput("paper.method", label = "Probability plot",
                      choices = paper.method)
        }else{
          selectInput("paper.method", label = "Probability plot",
                      choices = paper.method.pos)
        }
      })
      
      #確率紙のプロッティング公式
      output$paper.rank <- renderUI({
        selectInput("paper.rank", label = "Rank calculation",
                    choices = paper.rank)
      })
      
      #ヒストグラム表示
      output$gg.hist <- renderPlot({gg.hist(vec.data() )})
      
      #数値のまとめ表示
      output$vec.summary <- renderText({vec.summary(vec.data())})
      
      #一通り分布にあてはめる
      result <- reactive({
        
        #解析するデータを準備
        data.vec <-  vec.data()
        
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

            #計算リストにあるか確認
            if(match(distr.name, input$use) %>% is.na() == FALSE){
              
              #あてはめ
              ret[[distr.name]] <- fit.dist(data = data.vec, distr = distr.name, 
                                            method = input$fitdist.method,
                                            timeout = input$timeout)
              
            }

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
      
      #デバッグ用の表示
      output$fit.dist.res <- renderPrint({
        result()
      })
      
      #結果一覧表示
      output$result <- DT::renderDataTable({
        datatable(
          dt.result()[, -1],
          extensions = c('Buttons'), 
          filter = 'top',
          selection = "none",
          options = list(
            autoWidth = TRUE, pageLength = 100,
            dom = 'Blfrtip',buttons = c("copy", "csv"))
          ) %>% formatRound(columns = c(2:(ncol(dt.result()) - 1)), digits = 3)
        }, server = FALSE)
      
      #変数選択
      output$distr.sel <- renderUI({
        
        choices.distr <- dt.result()[, "distr"] %>% as.vec()
        names(choices.distr) <- dt.result()[, "name"] %>% as.vec()
        
        selectInput("distr.sel", label = "Distribution",
                    choices = choices.distr)
        
      })
      
      #結果表示
      output$result.plot <- renderPlot({
        plot.obj <- try(result()[[input$distr.sel]], silent = TRUE)
        if(class(plot.obj)[1] == "try-error"){return(NULL)}
        plot(plot.obj)
        })
 
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
      
      #確率紙プロット
      output$plot_paper <- renderPlot({
        plot_paper(result()[[input$distr.sel]], 
                   method = input$paper.method,
                   rank = input$paper.rank)
      })
      

    })
    
    
    
    
    
  })



})
