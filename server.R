library(shiny)
library(shinyWidgets)
library(plotly)
library(parallel)
library(DT)
library(colorspace)
library(googleVis)
library(longpower)
library(Hmisc)
library(lme4)
library(nlme)
library(arsenal)
library(dplyr)
library(foreign)

# Parameters ----
## Set colors ----
# https://www.w3.org/TR/css-color-3/#svg-color
confirmed_color <- "purple"
active_color <- "blue" # #1f77b4"
recovered_color <- "green"
death_color <- "red"

# To update data ----
# adnimerge <- ADNIMERGE::adnimerge
# adniDate <- packageDate('ADNIMERGE')
# csf2numeric <- function(x){
#   as.numeric(gsub('>', '', gsub('<', '', x)))
# }
# 
# adnimerge <- adnimerge %>%
#   mutate(across(c(ABETA.bl, TAU.bl), csf2numeric)) %>%
#   filter(Years.bl <= 5.25)
# 
# save(adnimerge, adniDate, file="ADNIMERGE.rdata")

load("./ADNIMERGE.rdata")

methods <- c(
  liuliang = 'Liu and Liang (1997)',
  diggle = 'Diggle et al (2002)',
  edland = 'Edland (2009)'
)

varCov <- c(
  exchangeable = "Compound symmetric, heterogeneous",
  general = "Unstructured, heterogeneous",
  ar1 = "AR1, heterogeneous"
)

# shinyServer ----
shinyServer(
  function(input, output, session){
    shinyalert("Disclaimer!", "This free software service is provided with ABSOLUTELY NO WARRANTY.",
      type = "warning", closeOnClickOutside=T, confirmButtonText = "I accept")
    
    ## Power/Sample size switch ----
    observeEvent(input$analysisType, {
      if(input$analysisType %in% c("Sample size")){
        shinyjs::enable("power")
        shinyjs::disable("sampleSize")
      }else{
        shinyjs::enable("sampleSize")
        shinyjs::disable("power")
      }
    })
    
    ## correlation matrix (R) ----
    observeEvent(input$matrix, {
      size<-length(seq(input$startTime, input$entTime, by=input$timeStep))
      M <- matrix(0,nrow=size, ncol = size, dimnames = list(paste0('t',1:size),paste0('t',1:size)))
      if(input$matrix %in% "correlation"){
        insertUI(
          selector = "#rmoo",
          where="afterEnd",
          ui=tags$div(id="R",tags$b("Pilot estimate of a marginal model working correlation matrix (R)"),
            matrixInput(inputId = "Rmatrix",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
      }else{
        removeUI(
          selector = "#R"
        )
      }
    })
    
    observeEvent(input$estimate, {
      if(input$estimate %in% c("delta")){
        shinyjs::enable("delta")
        shinyjs::disable("beta")
        shinyjs::disable("pct.change")
      }else{
        shinyjs::enable("beta")
        shinyjs::enable("pct.change")
        shinyjs::disable("delta")
      }
    })
    
    observeEvent(input$method, {
      if(input$method %in% c("edland")){
        shinyjs::enable("edlandAllocation")
      }else{
        shinyjs::disable("edlandAllocation")
      }
    })

    # Linear model ----
    ## Diggle ----
    R<-reactive({
      #cov.s.i <- #0.8*sqrt(input$sig2.i)*sqrt(input$sig2.s)
      cov.t <- function(t1, t2, sig2.i, sig2.s, cov.s.i){
        input$sig2.i + t1*t2*input$sig2.s + (t1+t2)*input$cov.s.i 
      }
      
      t = seq(input$startTime,input$entTime,input$timeStep)
      n = length(t)
      R = outer(t, t, function(x,y){cov.t(x,y, sig2.i, sig2.s, cov.s.i)})
      R = R + diag(input$sig2.e, n, n)
    })
    
    output$diggleSS<-renderPlotly({
      if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") &
          input$method %in% c("diggle")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        power<-seq(input$power[1],input$power[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimate %in% c("delta")){
            fit<-diggle.linear.power(delta=input$delta, t=t, R=R(), sig.level=input$alpha, 
              alternative=input$alternative, power=power[i]/100) #, sigma2 = input$sig2.e
          } else{
            delta<-input$beta*input$pct.change
            fit<-diggle.linear.power(delta=delta, t=t, R=R(), sig.level=input$alpha, 
              alternative=input$alternative, power=power[i]/100) #, sigma2 = input$sig2.e
          }
          ssize[i]<-fit$N
        }
        dat<-data.frame(n=ssize, power=power)
        
        plot_ly(data=dat, x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]), 
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else if(input$analysisType %in% c("Power") & #input$matrix %in% c("covariance") &
          input$method %in% c("diggle")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        n<-seq(input$sampleSize[1],input$sampleSize[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimate %in% c("delta")){
            
            fit<-lmmpower(delta=input$delta, t=t, R=R(), sig.level=input$alpha, method="diggle",
              alternative=input$alternative,n=n[i])
          } else{
            delta<-input$beta*input$pct.change
            fit<-lmmpower(delta=delta, t=t, R=R(), sig.level=input$alpha,method="diggle", 
              alternative=input$alternative, n=n[i])
          }
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        
        plot_ly(data=dat, x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]),
            xaxis=list(title="Sample size per group"), yaxis=list(title="Power (%)"))
      
      }else if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") &
          ## Edland Power ----
          input$method %in% c("edland")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        power<-seq(input$power[1],input$power[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimate %in% c("delta")){
            fit<-edland.linear.power(delta=input$delta, t=t, sig2.s = input$sig2.s,sig2.int = input$sig2.i, sig2.e =input$sig2.e,
              sig.level=input$alpha, alternative=input$alternative, power =power[i]/100, 
              lambda = input$edlandAllocation)
          } else {
            delta<-input$beta*input$pct.change
            fit<-edland.linear.power(delta=delta, t=t, sig2.s = input$sig2.s, sig2.e =input$sig2.e,
              sig2.int = input$sig2.i,sig.level=input$alpha, alternative=input$alternative, 
              power =power[i]/100, lambda = input$edlandAllocation)
          }
          ssize[i]<-fit$N
        }
        dat<-data.frame(n=ssize, power=power)
        
        plot_ly(data=dat, x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]), 
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else if(input$analysisType %in% c("Power") & #input$matrix %in% c("covariance") &
          ### Edland SS ----
          input$method %in% c("edland")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        n<-seq(input$sampleSize[1],input$sampleSize[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimate %in% c("delta")){
            
            fit<-edland.linear.power(delta=input$delta, t=t, sig2.s = input$sig2.s,
              sig2.int = input$sig2.i, sig2.e =input$sig2.e,
              sig.level=input$alpha, alternative=input$alternative, n =n[i], 
              lambda = input$edlandAllocation)
          } else{
            delta<-input$beta*input$pct.change
            fit<-edland.linear.power(delta=input$delta, t=t, sig2.s = input$sig2.s,
              sig2.int = input$sig2.i, sig2.e =input$sig2.e,
              sig.level=input$alpha, alternative=input$alternative, n =n[i], 
              lambda = input$edlandAllocation)
          }
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        
        plot_ly(data=dat, x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]),
            xaxis=list(title="Sample size, group 1"), yaxis=list(title="Power (%)"))
        
      }else if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") & 
          ### Liu, Liang Power ----
          input$method %in% c("liuliang")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        n = length(t)
        u = list(u1 = t, u2 = rep(0,n))
        v = list(v1 = cbind(1,1,t),
          v2 = cbind(1,0,t))
        
        power<-seq(input$power[1],input$power[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimate %in% c("delta")){
            fit<-liu.liang.linear.power(delta=input$delta, u=u, v=v, R=R(), sig.level=input$alpha, 
              alternative=input$alternative,power=power[i]/100) 
          }else {
            delta<-input$beta*input$pct.change
            fit<-liu.liang.linear.power(delta=delta, u=u, v=v, R=R(), sig.level=input$alpha, 
              alternative=input$alternative,power=power[i]/100) 
          }
          
          ssize[i]<-fit$N
        }
        dat<-data.frame(n=ssize, power=power)
        
        plot_ly(data=dat, x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else if(input$analysisType %in% c("Power") & #input$matrix %in% c("covariance") & 
          ### Liu, Liang SS ----
          input$method %in% c("liuliang")){
        t = seq(input$startTime,input$entTime,input$timeStep)
        n = length(t)
        u = list(u1 = t, u2 = rep(0,n))
        v = list(v1 = cbind(1,1,t),
          v2 = cbind(1,0,t))
        
        n<-seq(input$sampleSize[1],input$sampleSize[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimate %in% c("delta")){
            fit<-liu.liang.linear.power(delta=input$delta, u=u, v=v, R=R(), sig.level=input$alpha, 
              alternative=input$alternative,N=n[i]) 
          }else {
            delta<-input$beta*input$pct.change
            fit<-liu.liang.linear.power(delta=delta, u=u, v=v, R=R(), sig.level=input$alpha, 
              alternative=input$alternative,N=n[i]) 
          }
          
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        
        plot_ly(data=dat, x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste("Power analysis using", methods[input$method]),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
      }
    })
    
    ## Describe method ----
    output$describeMethod<-renderText({
      if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") &
          input$method %in% c("diggle")){
        "Sample size calculation for difference in slopes between two groups. See
          Diggle et al (2002) for parameter definitions and other details."
      }else if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") &
          input$method %in% c("edland")){
        "Sample size/power calculation for a linear mixed model with only random slope 
        (simply by setting sig2.i = 0). See Edland (2009) for details."
      }else if(input$analysisType %in% c("Sample size") & #input$matrix %in% c("covariance") &
          input$method %in% c("liuliang")){
        "Sample size calculation for a linear mixed model. See Liu and Liang
        (1997) for parameter definitions and other details and the 'longpower' package vignette for more details."
      }
    })
    
    ## Summary selection ----
    output$summarySelection<-renderTable({
      if(input$estimate %in% c("delta")){
        if(input$method %in% c("edland")){
          a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error", "Placebo change",
            "Change in the estimate of the parameter of interest:", "Allocation ratio"),
            Value=as.character(c(input$method,input$alternative, input$alpha, input$beta, input$delta,
              input$edlandAllocation))) 
        }else{
          a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error", "Placebo change",
            "Change in the estimate of the parameter of interest:"),
            Value=as.character(c(input$method,input$alternative, input$alpha, input$beta, input$delta)))
        }
      } else{
        if(input$method %in% c("edland")){
          a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error", "Placebo change",
            "Change in the estimate of the parameter of interest:", "Allocation ratio"),
            Value=as.character(c(input$method,input$alternative, input$alpha, input$beta, 
              input$beta*input$pct.change,input$edlandAllocation))) 
        }else{
          a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error", "Placebo change",
            "Change in the estimate of the parameter of interest:"),
            Value=as.character(c(input$method,input$alternative, input$alpha, input$beta, 
              input$beta*input$pct.change)))
        }
        
      }
    })

    # MMRM ----
    ## Power/Sample size switch ----
    observeEvent(input$analysisTypeMMRM, {
      if(input$analysisTypeMMRM %in% c("Sample size")){
        shinyjs::enable("powerMMRM")
        shinyjs::disable("sampleSizeMMRM")
      }else{
        shinyjs::enable("sampleSizeMMRM")
        shinyjs::disable("powerMMRM")
      }
    })
    
    size<-reactive({
      input$timePoints
    })

    observeEvent(input$matrixMMRM, {
      M <- matrix(input$rhoMMRM,nrow=size(), ncol = size(), 
        dimnames = list(paste0('t',1:size()),paste0('t',1:size())))
      diag(M)<-1
      m <- matrix(1,nrow=size(), ncol = 1, dimnames = list(paste0('t',1:size()),NULL))
      
      if(input$matrixMMRM %in% c("exchangeable")){
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        # removeUI(
        #   selector = "#lambdaExMMRM"
        # )
        insertUI(
          selector = "#RaMMRM",
          where="afterEnd",
          ui=tags$div(id="RaExMMRM",tags$b("Working correlation matrix (Ra)"),
            matrixInput(inputId = "RamatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#RbMMRM",
          where="afterEnd",
          ui=tags$div(id="RbExMMRM",tags$b("Working correlation matrix (Rb)"),
            matrixInput(inputId = "RbmatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
      }else if(input$matrixMMRM %in% c("general")){
        M <- matrix(0,nrow=size(), ncol = size(), dimnames = list(paste0('t',1:size()),paste0('t',1:size())))
        diag(M)<-1
        
        #m <- matrix(1,nrow=size, ncol = 1, dimnames = list(paste0('t',1:size),NULL))
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        # removeUI(
        #   selector = "#lambdaExMMRM"
        # )
        insertUI(
          selector = "#RaMMRM",
          where="afterEnd",
          ui=tags$div(id="RaExMMRM",tags$b("Working correlation matrix (Ra)"),
            matrixInput(inputId = "RamatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#RbMMRM",
          where="afterEnd",
          ui=tags$div(id="RbExMMRM",tags$b("Working correlation matrix (Rb)"),
            matrixInput(inputId = "RbmatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
      }else{
        m2 <- matrix(0,nrow=size()*2, ncol = 1, dimnames = list(paste0('t',1:(size()*2)),NULL))
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        # removeUI(
        #   selector = "#lambdaExMMRM"
        # )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
        
        # removeUI(
        #   selector = "#lambdaMMRM"
        # )
      }
    })

    observeEvent(input$updateMMRM, {
      M <- matrix(input$rhoMMRM,nrow=size(), ncol = size(), 
        dimnames = list(paste0('t',1:size()),paste0('t',1:size())))
      diag(M)<-1
      m <- matrix(1,nrow=size(), ncol = 1, dimnames = list(paste0('t',1:size()),NULL))
      
      if(input$matrixMMRM %in% c("exchangeable")){
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        insertUI(
          selector = "#RaMMRM",
          where="afterEnd",
          ui=tags$div(id="RaExMMRM",tags$b("Working correlation matrix (Ra)"),
            matrixInput(inputId = "RamatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#RbMMRM",
          where="afterEnd",
          ui=tags$div(id="RbExMMRM",tags$b("Working correlation matrix (Rb)"),
            matrixInput(inputId = "RbmatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
      }else if(input$matrixMMRM %in% c("general")){
        M <- matrix(0,nrow=size(), ncol = size(), dimnames = list(paste0('t',1:size()),paste0('t',1:size())))
        diag(M)<-1
        
        #m <- matrix(1,nrow=size, ncol = 1, dimnames = list(paste0('t',1:size),NULL))
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        # removeUI(
        #   selector = "#lambdaExMMRM"
        # )
        insertUI(
          selector = "#RaMMRM",
          where="afterEnd",
          ui=tags$div(id="RaExMMRM",tags$b("Working correlation matrix (Ra)"),
            matrixInput(inputId = "RamatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#RbMMRM",
          where="afterEnd",
          ui=tags$div(id="RbExMMRM",tags$b("Working correlation matrix (Rb)"),
            matrixInput(inputId = "RbmatrixMMRM",value = M, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
      }else{
        m2 <- matrix(0,nrow=size()*2, ncol = 1, dimnames = list(paste0('t',1:(size()*2)),NULL))
        removeUI(
          selector = "#RaExMMRM"
        )
        removeUI(
          selector = "#raExMMRM"
        )
        removeUI(
          selector = "#RbExMMRM"
        )
        removeUI(
          selector = "#rbExMMRM"
        )
        # removeUI(
        #   selector = "#lambdaExMMRM"
        # )
        insertUI(
          selector = "#raMMRM",
          where="afterEnd",
          ui=tags$div(id="raExMMRM",tags$b("Retention in group a (ra)"),
            matrixInput(inputId = "ramatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        insertUI(
          selector = "#rbMMRM",
          where="afterEnd",
          ui=tags$div(id="rbExMMRM",tags$b("Retention in group b (rb)"),
            matrixInput(inputId = "rbmatrixMMRM",value = m, rows = list(names = T),
              cols = list(names = T), class = 'numeric'))
        )
        # insertUI(
        #   selector = "#lambdaMMRM",
        #   where="afterEnd",
        #   ui=tags$div(id="lambdaExMMRM",tags$b("Allocation ratio(lambda)"),
        #               matrixInput(inputId = "lambdaVecMMRM",value = t(m2), rows = list(names = T),
        #                           cols = list(names = T), class = 'numeric'))
        # )
        
        # removeUI(
        #   selector = "#lambdaMMRM"
        # )
      }
    })
    
    observeEvent(input$estimateMMRM, {
      if(input$estimateMMRM %in% c("delta")){
        shinyjs::enable("deltaMMRM")
        shinyjs::disable("betaMMRM")
        shinyjs::disable("pct.changeMMRM")
      }else{
        shinyjs::enable("betaMMRM")
        shinyjs::enable("pct.changeMMRM")
        shinyjs::disable("deltaMMRM")
      }
    })
    
    observeEvent(input$matrixMMRM, {
      if(input$matrixMMRM %in% c("exchangeable")){
        shinyjs::enable("rhoMMRM")
        shinyjs::enable("lambdaMMRM")
        #shinyjs::enable("Rmatrix")
      }else if(input$matrixMMRM %in% c("ar1")){
        shinyjs::enable("rhoMMRM")
        shinyjs::enable("lambdaMMRM")
        #shinyjs::enable("RaExMMRM")
      } else if(input$matrixMMRM %in% c("general")){
        shinyjs::disable("rhoMMRM")
        shinyjs::enable("lambdaMMRM")
      }
    })

    ## Output ----
    output$mmrmplot<-renderPlotly({
      if(input$matrixMMRM %in% c("exchangeable") & input$analysisTypeMMRM %in% c("Sample size")){
        power<-seq(input$powerMMRM[1],input$powerMMRM[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =input$deltaMMRM, power =power[i]/100,
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =delta, power =power[i]/100,
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          }
          ssize[i]<-fit$n1+fit$n2
        }
        dat<-data.frame(n=ssize, power=power)
        plot_ly(data=dat[-c(1:2,nrow(dat)),], x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else  if(input$matrixMMRM %in% c("exchangeable") & input$analysisTypeMMRM %in% c("Power")){
        n<-seq(input$sampleSizeMMRM[1],input$sampleSizeMMRM[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =input$deltaMMRM, N =n[i],
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =delta,N =n[i],
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          }
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        plot_ly(data=dat[-c(1:2,nrow(dat)),], x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else if(input$matrixMMRM %in% c("general") & input$analysisTypeMMRM %in% c("Sample size")){
        power<-seq(input$powerMMRM[1],input$powerMMRM[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =input$deltaMMRM, power =power[i]/100,
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =delta, power =power[i]/100,
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          }
          ssize[i]<-fit$n1+fit$n2
        }
        dat<-data.frame(n=ssize, power=power)
        plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      }else if(input$matrixMMRM %in% c("general") & input$analysisTypeMMRM %in% c("Power")){
        n<-seq(input$sampleSizeMMRM[1],input$sampleSizeMMRM[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =input$deltaMMRM, N =n[i],
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm(Ra = input$RamatrixMMRM, ra = input$ramatrixMMRM, sigmaa = input$sigmaaMMRM, 
              sigmab = input$sigmabMMRM,
              lambda =input$lambdaMMRM,
              delta =delta, N =n[i],
              sig.level=input$alphaMMRM, alternative=input$alternativeMMRM)
          }
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      } else if(input$matrixMMRM %in% c("ar1")& input$analysisTypeMMRM %in% c("Sample size")){
        power<-seq(input$powerMMRM[1],input$powerMMRM[2],1)
        ssize<-rep(0,length(power))
        for(i in 1:length(power)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm.ar1(rho=input$rhoMMRM, ra=input$ramatrixMMRM, sigmaa=input$sigmaaMMRM, 
              rb = input$rbmatrixMMRM,alternative = input$alternativeMMRM,
              lambda =input$lambdaMMRM, power = power[i]/100, delta =input$deltaMMRM,
              sig.level=input$alphaMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm.ar1(rho=input$rhoMMRM, ra=input$ramatrixMMRM, sigmaa=input$sigmaaMMRM, 
              rb = input$rbmatrixMMRM,alternative = input$alternativeMMRM,
              lambda =input$lambdaMMRM, power = power[i]/100, delta =delta,
              sig.level=input$alphaMMRM)
          }
          ssize[i]<-fit$n1+fit$n2
        }
        dat<-data.frame(n=ssize, power=power)
        plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
        
      } else if(input$matrixMMRM %in% c("ar1")& input$analysisTypeMMRM %in% c("Power")){
        n<-seq(input$sampleSizeMMRM[1],input$sampleSizeMMRM[2],1)
        pw<-rep(0,length(n))
        for(i in 1:length(n)){
          if(input$estimateMMRM %in% c("delta")){
            fit<-power.mmrm.ar1(rho=input$rhoMMRM, ra=input$ramatrixMMRM, sigmaa=input$sigmaaMMRM, 
              rb = input$rbmatrixMMRM,alternative = input$alternativeMMRM,
              lambda =input$lambdaMMRM, N = n[i], delta =input$deltaMMRM,
              sig.level=input$alphaMMRM)
          } else {
            delta<-input$betaMMRM*input$pct.changeMMRM
            fit<-power.mmrm.ar1(rho=input$rhoMMRM, ra=input$ramatrixMMRM, sigmaa=input$sigmaaMMRM, 
              rb = input$rbmatrixMMRM,alternative = input$alternativeMMRM,
              lambda =input$lambdaMMRM, N = n[i], delta =delta,
              sig.level=input$alphaMMRM)
          }
          pw[i]<-fit$power
        }
        dat<-data.frame(n=n, power=pw)
        plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power*100, type = 'scatter', mode = 'lines') %>%
          layout(title=paste0("Power analysis for MMRM"),
            xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
      }
    })
    
    output$describeMethodMMRM<-renderText({
      if(input$matrixMMRM %in% c("exchangeable")){
        "Sample size/power calculation for a mixed model of repeated measures with
       exhangeable correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other
           details."
      } else if(input$matrixMMRM %in% c("general")){
        "Sample size/power calculation for a mixed model of repeated measures with
       general correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other
           details."
      }else{
        "Sample size/power calculation for a mixed model of repeated measures with
        AR(1) correlation structure. See Lu, Luo, & Chen (2008) for parameter definitions and other details"
      }
    })
    
    output$summarySelectionMMRM<-renderTable({
      if(input$estimateMMRM %in% c("deltaMMRM")){
        a<-data.frame(Summary=c("Var-cov structure", "Number of time points","Type of test", "Type I error", "Placebo change",
          "Effect size:", "Allocation ratio: "),
          Value=as.character(c(varCov[input$matrixMMRM], input$timePoints, input$alternative, input$alpha, input$beta,
            input$delta, input$lambdaMMRM)))
      } else{
        a<-data.frame(Summary=c("Var-cov structure", "Number of time points","Type of test", "Type I error", "Placebo change",
          "Effect size:", "Allocation ratio"),
          Value=as.character(c(varCov[input$matrixMMRM], input$timePoints, input$alternative, input$alpha, input$beta,
            input$beta*input$pct.change, input$lambdaMMRM)))
      }
    })
    
    # ADNI GENERATOR (slope model) ----
    output$selectAge<-renderUI({
      age<-unique(adnimerge$AGE)
      sliderInput("age", label = "Age", min = 50, max = 100, 
        value = c(53, ceiling(max(age, na.rm = T))), step = 1)
    })
    
    output$selectAV45<-renderUI({
      AV45.bl<-round(unique(adnimerge$AV45.bl),2)
      sliderInput("av45bl", label = "AV45 (SUVR)", min = 0, max = max(AV45.bl, na.rm = T), 
        value = c(0, round(max(AV45.bl, na.rm = T),2)), step = 0.01)
    })
    
    output$selectAbeta<-renderUI({
      ABETA.bl<-unique(as.numeric(adnimerge$ABETA.bl))
      sliderInput("abetabl", label = "Aβ42", min = min(ABETA.bl, na.rm = T), max = ceiling(max(ABETA.bl, na.rm = T)), 
        value = c(min(ABETA.bl, na.rm = T), round(max(ABETA.bl, na.rm = T),1)))
    })
    
    output$selectTau<-renderUI({
      TAU.bl<-unique(as.numeric(adnimerge$TAU.bl))
      sliderInput("taubl", label = "Tau", min = min(TAU.bl, na.rm = T), max = ceiling(max(TAU.bl, na.rm = T)), 
        value = c(min(TAU.bl, na.rm = T), round(max(TAU.bl, na.rm = T),1)))
    })
    
    output$selectRatio<-renderUI({
      ratio.bl<-unique(as.numeric(adnimerge$TAU.bl)/as.numeric(adnimerge$ABETA.bl))
      sliderInput("tauRatio", label = "Tau/Aβ42 ratio", min = 0, 
        max = round(max(ratio.bl, na.rm = T),1), step = 0.01,
        value = c(0, round(max(ratio.bl, na.rm = T),1)))
    })
    
    output$selectDx<-renderUI({
      Dx<-levels(adnimerge$DX.bl)
      checkboxGroupInput("dx", label =HTML(paste0("Diagnosis (",tags$em("Do not deselect all"),")")),
        choices = Dx, selected = Dx, inline = T)
    })
    
    output$selectMMSE<-renderUI({
      MMSE.bl<-adnimerge$MMSE.bl
      MMSE.bl<-unique(as.numeric(MMSE.bl))
      sliderInput("mmsebl", label = "MMSE", min = min(MMSE.bl, na.rm = T), max = ceiling(max(MMSE.bl, na.rm = T)), 
        value = c(min(MMSE.bl, na.rm = T), round(max(MMSE.bl, na.rm = T),1)))
    })
    
    output$selectCDRSB<-renderUI({
      CDRSB.bl<-adnimerge$CDRSB.bl
      CDRSB.bl<-unique(as.numeric(CDRSB.bl))
      sliderInput("cdrsbbl", label = "CDRSB", min = min(CDRSB.bl, na.rm = T), max = ceiling(max(CDRSB.bl, na.rm = T)), 
        value = c(min(CDRSB.bl, na.rm = T), round(max(CDRSB.bl, na.rm = T),1)))
    })
    
    output$selectLogMem<-renderUI({
      LDELTOTAL.bl<-adnimerge$LDELTOTAL.bl
      LDELTOTAL.bl<-unique(as.numeric(LDELTOTAL.bl))
      sliderInput("logmembbl", label = "Logical Memory", min = min(LDELTOTAL.bl, na.rm = T), max = ceiling(max(LDELTOTAL.bl, na.rm = T)), 
        value = c(min(LDELTOTAL.bl, na.rm = T), round(max(LDELTOTAL.bl, na.rm = T),1)))
    })
    
    output$selectDuration<-renderUI({
      sliderInput("studyDuration", label = "Select study duration (years)", min = 1, 
        max = 5, value = 1.5, step=0.5)
    })
    
    ## Grayout criteria selection ----
    observeEvent(input$submitCriteria, {
      if(c("AGE") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("age")
      }else{
        shinyjs::disable("age")
      }
      if(c("ABETA.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("abetabl")
      }else{
        shinyjs::disable("abetabl")
      }
      
      if(c("APOE4") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("APOE4")
      }else{
        shinyjs::disable("APOE4")
      }
      
      if(c("TAU.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("taubl")
      }else{
        shinyjs::disable("taubl")
      }
      
      if(c("DX.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("dx")
      }else{
        shinyjs::disable("dx")
      }
      
      if(c("TAU/ABETA") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("tauRatio")
      }else{
        shinyjs::disable("tauRatio")
      }
      
      if(c("AV45.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("av45bl")
      }else{
        shinyjs::disable("av45bl")
      }
      
      if(c("LDELTOTAL.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("logmembbl")
      }else{
        shinyjs::disable("logmembbl")
      }
      
      if(c("MMSE.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("mmsebl")
      }else{
        shinyjs::disable("mmsebl")
      }
      
      
      if(c("CDRSB.bl") %in% input$criteria | is.null(input$criteria)){
        shinyjs::enable("cdrsbbl")
      }else{
        shinyjs::disable("cdrsbbl")
      }
      
    })
    
    observeEvent(input$methodADNI, {
      if(input$methodADNI %in% c("edland")){
        shinyjs::enable("edlandAllocationADNI")
      }else{
        shinyjs::disable("edlandAllocationADNI")
      }
    })

    ## Default calculation ----
    subdata<-reactive({
      a<-c("","AGE","ABETA.bl","APOE4","TAU.bl","DX.bl","TAU/ABETA","AV45.bl","LDELTOTAL.bl",
        "MMSE.bl","CDRSB.bl")
      
      if(identical(as.character(input$criteria), a)){
        if("Positive" %in% as.character(input$APOE4) & "Negative" %in% as.character(input$APOE4)){
          apoe<-c(0,1,2)
        }else if("Positive" %in% as.character(input$APOE4) & !"Negative" %in% as.character(input$APOE4)){
          apoe<-c(1,2)
        }else if(!"Positive" %in% as.character(input$APOE4) & "Negative" %in% as.character(input$APOE4)){
          apoe<-c(0)
        }
        adni<-adnimerge %>%
          subset(AGE >= input$age[1] & AGE <= input$age[2])%>%
          subset(ABETA.bl >= input$abetabl[1] & ABETA.bl <= input$abetabl[2] | ABETA.bl %in% NA_real_) %>%
          subset(APOE4 %in% c(0,1,2) |APOE4 %in% NA_integer_)%>%
          subset(TAU.bl >= input$taubl[1] & TAU.bl <= input$taubl[2]| TAU.bl %in% NA_real_) %>%
          subset(DX.bl %in% as.character(input$dx)| DX.bl %in% NA_character_) %>%
          subset(round(TAU.bl/ABETA.bl,1) >= input$tauRatio[1] & round(TAU.bl/ABETA.bl,1) <= input$tauRatio[2]| 
              round(TAU.bl/ABETA.bl,1) %in% NA_real_) %>%
          subset(AV45.bl >= input$av45bl[1] & AV45.bl <= input$av45bl[2]| AV45.bl %in% NA_real_) %>%
          subset(LDELTOTAL.bl >= input$logmembbl[1] & LDELTOTAL.bl <= input$logmembbl[2]| LDELTOTAL.bl %in% NA_real_) %>%
          subset(MMSE.bl >= input$mmsebl[1] & MMSE.bl <= input$mmsebl[2]| MMSE.bl %in% NA_real_) %>%
          subset(CDRSB.bl >= input$cdrsbbl[1] & CDRSB.bl <= input$cdrsbbl[2]| CDRSB.bl %in% NA_real_) %>%
          subset(Years.bl <= input$studyDuration + 0.25| Years.bl %in% NA_real_) %>%
          filter(M <=  input$studyDuration*12) %>%
          filter(APOE4 %in% apoe| APOE4 %in% NA_integer_)
        
      }else {
        adni<-adnimerge %>%
          filter(Years.bl <= input$studyDuration + 0.25 | Years.bl %in% NA_real_) %>%
          filter(M <=  input$studyDuration*12)
        if("AGE" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(AGE >= input$age[1] & AGE <= input$age[2])
        }
        if("ABETA.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(ABETA.bl >= input$abetabl[1] & ABETA.bl <= input$abetabl[2]| ABETA.bl %in% NA_real_)
        }
        if("APOE4" %in% as.character(input$criteria)){
          if("Positive" %in% as.character(input$APOE4) & "Negative" %in% as.character(input$APOE4)){
            apoe<-c(0,1,2)
          }else if("Positive" %in% as.character(input$APOE4) & !"Negative" %in% as.character(input$APOE4)){
            apoe<-c(1,2)
          }else if(!"Positive" %in% as.character(input$APOE4) & "Negative" %in% as.character(input$APOE4)){
            apoe<-c(0)
          }
          adni<-adni %>%
            filter(APOE4 %in% apoe| APOE4 %in% NA_integer_)
        }
        if("TAU.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(TAU.bl >= input$taubl[1] & TAU.bl <= input$taubl[2]| TAU.bl %in% NA_real_) 
        }
        if("DX.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(DX.bl %in% as.character(input$dx)| DX.bl %in% NA_character_)
        }
        if("TAU/ABETA" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(round(TAU.bl/ABETA.bl,1) >= input$tauRatio[1] & round(TAU.bl/ABETA.bl,1) <= input$tauRatio[2]| 
                round(TAU.bl/ABETA.bl,1) %in% NA_real_)
        }
        if("AV45.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(AV45.bl >= input$av45bl[1] & AV45.bl <= input$av45bl[2]| AV45.bl %in% NA_real_)
        }
        if("LDELTOTAL.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(LDELTOTAL.bl >= input$logmembbl[1] & LDELTOTAL.bl <= input$logmembbl[2]| LDELTOTAL.bl %in% NA_real_)
        }
        if("MMSE.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(MMSE.bl >= input$mmsebl[1] & MMSE.bl <= input$mmsebl[2]| MMSE.bl %in% NA_real_)
        }
        if("CDRSB.bl" %in% as.character(input$criteria)){
          adni<-adni %>%
            filter(CDRSB.bl >= input$cdrsbbl[1] & CDRSB.bl <= input$cdrsbbl[2]| CDRSB.bl %in% NA_real_)
        }
      }
      adni<-adni %>%
        mutate(RID=factor(RID))
    })
    
    ## Baseline summaries ----
    observeEvent(input$summaryby,{
      output$baselineSummary<-renderTable({
        dat<-subdata() %>%
          arrange(Years.bl)
        a<-as.character(input$criteria)[-1]
        a[a=="TAU/ABETA"]<-"TAU_ABETARatio"
        if("TAU/ABETA" %in% input$criteria){
          dat$TAU_ABETARatio<-with(dat, TAU.bl/ABETA.bl)
        }
        if("APOE4" %in% a){
          dat <-dat %>%
            mutate(APOE4=if_else(APOE4 %in% c(1,2), "Positive",
              if_else(APOE4 %in% c(0), "Negative", NA_character_)))
        }
        
        dat<-dat[!duplicated(dat$RID),]
        label(dat$AV45.bl)<-"AV45.bl"
        
        cov<-paste(a, collapse="+")
        cov<-paste0("~",cov)
        
        regFormula<- reactive({
          as.formula(paste(input$summaryby,cov))
        })
        tbl<-arsenal::tableby(regFormula(),data=dat,
          numeric.stats = c("N","Nmiss","range","meansd", "medianq1q3"),
          cat.stats = c("N","Nmiss", "countpct"),
          digits=1)
        as.data.frame(summary(tbl, text = "html"))
      }, sanitize.text.function = function(x) x)
    })
    
    output$indPlot<-renderPlotly({
      p <- ggplot(data = subdata(), aes(x = Years.bl, 
        y = get(input$longout), group = RID))+
        geom_line(alpha=0.25)+theme_bw()+xlab("Years from baseline")+
        ylab(as.character(input$longout))
      ggplotly(p)
    })

    output$plotProfile<-renderPlotly({
      p<-ggplot(subdata(), aes(x=Years.bl,y=get(input$longout))) +
        geom_point(alpha=0.25) +
        geom_smooth()+theme_bw()+xlab("Years from baseline")+
        ylab(as.character(input$longout))
      ggplotly(p)
    })
    
    ## Model fit ----
    modelFit<-reactive({
      if(is.null(input$selectCovariates)){
        m1 <-lmer(formula(paste(input$longout, "~ Years.bl + (Years.bl |RID)")), subdata())
        m1
      }else{
        cov<-paste(input$selectCovariates,collapse = "+")
        m1 <-lmer(formula(paste(input$longout,"~",cov,"+ Years.bl + (Years.bl |RID)")), subdata())
        m1
      }
    })
    
    slope<-reactive({
      round(summary(modelFit())$coef['Years.bl', 'Estimate'],4)
    })
    
    output$modelFitSummary<-renderPrint({
      print(summary(modelFit()))
    })
    
    output$digglePlot<-renderPlotly({
      if(is.null(input$selectCovariates)){
        m1 <-lmer(formula(paste(input$longout, "~ Years.bl + (Years.bl |RID)")), subdata())
      }else{
        cov<-paste(input$selectCovariates,collapse = "+")
        m1 <-lmer(formula(paste(input$longout,"~",cov,"+ Years.bl + (Years.bl |RID)")), subdata())
      }
      
      t<-seq(0, input$studyDuration, by=input$timeStepADNI)
      
      pct.change<-input$pchangeADNI
      power<- seq(input$powerADNI[1],input$powerADNI[2],1)
      
      ssize<-rep(0,length(power))
      for(i in 1:length(power)){
        if(input$methodADNI %in% c("edland")){
          fit<-lmmpower(m1, beta=summary(m1)$coef['Years.bl', 'Estimate'], t=t, sig.level=input$alphaADNI, 
            alternative=input$alternativeADNI, power=power[i]/100, method=as.character(input$methodADNI),
            pct.change=pct.change/100, lambda=input$edlandAllocationADNI)
        } else{
          fit<-lmmpower(m1, beta=summary(m1)$coef['Years.bl', 'Estimate'], t=t, sig.level=input$alphaADNI, 
            alternative=input$alternativeADNI, power=power[i]/100, method=as.character(input$methodADNI),
            pct.change=pct.change/100) 
        }
        ssize[i]<-fit$n[1]+fit$n[2]
      }
      
      dat<-data.frame(n=ssize, power=power)
      plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
        layout(title=paste("Power analysis using", methods[input$methodADNI]),
          xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
    })
    
    output$describeMethodADNI<-renderText({
      if(input$methodADNI %in% c("diggle")){
        "Sample size calculation for difference in slopes between two groups using ADNI data to estimate pilot parameters. See
          Diggle et al (2002) for parameter definitions and other details."
      }else if(input$methodADNI %in% c("edland")){
        "Sample size/power calculation for a linear mixed model with only random slope 
        (simply by setting sig2.i = 0). See Edland (2009) for details."
      }else if(input$methodADNI %in% c("liuliang")){
        "Sample size calculation for a linear mixed model. See Liu and Liang
        (1997) for parameter definitions and other details and the 'longpower' package vignette for more details."
      }
    })
    
    ## summary table ----
    output$summarySelectionADNI<-renderTable({
      
      if(input$methodADNI %in% c("edland")){
        a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error",
          "Estimated rate of change in placebo group", 
          "Effect size (% of placebo group change)", 
          "Effect size (raw scale)",
          "Allocation ratio",
          "Observation times"),
          Value=as.character(c(input$methodADNI,input$alternativeADNI, input$alphaADNI, 
            slope(), 
            input$pchangeADNI,
            slope()*input$pchangeADNI/100, 
            input$edlandAllocationADNI,  
            paste(seq(0, input$studyDuration, by=input$timeStepADNI), collapse=', ')))) 
      }else{
        a<-data.frame(Summary=c("Sample size method", "Type of test", "Type I error",
          "Estimated rate of change in placebo group", 
          "Effect size (% of placebo group change)", 
          "Effect size (raw scale)",
          "Observation times"),
          Value=as.character(c(input$methodADNI,
            input$alternativeADNI, 
            input$alphaADNI,
            slope(), 
            input$pchangeADNI,
            slope()*input$pchangeADNI/100,  
            paste(seq(0, input$studyDuration, by=input$timeStepADNI), collapse=', '))))
      }
    })
    
    # ADNI GENERATOR (MMRM) ----
    output$selectAgeMMRM<-renderUI({
      age<-round(unique(adnimerge$AGE),1)
      sliderInput("ageMMRM", label = "Age", min =50, max = 100, 
        value = c(53, ceiling(max(age, na.rm = T))), step = 1)
    })
    
    output$selectAV45MMRM<-renderUI({
      AV45.bl<-round(unique(adnimerge$AV45.bl),2)
      sliderInput("av45blMMRM", label = "AV45 (SUVR)", min = 0, max = max(AV45.bl, na.rm = T), 
        value = c(0, round(max(AV45.bl, na.rm = T),2)), step = 0.01)
    })
    
    output$selectAbetaMMRM<-renderUI({
      ABETA.bl<-unique(as.numeric(adnimerge$ABETA.bl))
      sliderInput("abetablMMRM", label = "Aβ42 (pg/ml)", min = min(ABETA.bl, na.rm = T), max = ceiling(max(ABETA.bl, na.rm = T)), 
        value = c(min(ABETA.bl, na.rm = T), round(max(ABETA.bl, na.rm = T),1)))
    })
    
    output$selectTauMMRM<-renderUI({
      TAU.bl<-unique(as.numeric(adnimerge$TAU.bl))
      sliderInput("taublMMRM", label = "Tau (pg/ml)", min = min(TAU.bl, na.rm = T), max = ceiling(max(TAU.bl, na.rm = T)), 
        value = c(min(TAU.bl, na.rm = T), round(max(TAU.bl, na.rm = T),1)))
    })
    
    output$selectRatioMMRM<-renderUI({
      ratio.bl<-unique(as.numeric(adnimerge$TAU.bl)/as.numeric(adnimerge$ABETA.bl))
      sliderInput("tauRatioMMRM", label = "Tau/Aβ42 ratio", min = min(0, na.rm = T), 
        max = round(max(ratio.bl, na.rm = T),1), step = 0.01,
        value = c(0, round(max(ratio.bl, na.rm = T),1)))
    })
    
    output$selectDxMMRM<-renderUI({
      Dx<-levels(adnimerge$DX.bl)
      checkboxGroupInput("dxMMRM", label = HTML(paste0("Diagnosis (",tags$em("Do not deselect all"),")")), 
        choices = Dx, selected = Dx, inline = T)
    })
    
    output$selectMMSEMMRM<-renderUI({
      MMSE.bl<-adnimerge$MMSE.bl
      MMSE.bl<-unique(as.numeric(MMSE.bl))
      sliderInput("mmseblMMRM", label = "MMSE", min = min(MMSE.bl, na.rm = T), max = ceiling(max(MMSE.bl, na.rm = T)), 
        value = c(min(MMSE.bl, na.rm = T), round(max(MMSE.bl, na.rm = T),1)))
    })
    
    output$selectCDRSBMMRM<-renderUI({
      CDRSB.bl<-adnimerge$CDRSB.bl
      CDRSB.bl<-unique(as.numeric(CDRSB.bl))
      sliderInput("cdrsbblMMRM", label = "CDRSB", min = min(CDRSB.bl, na.rm = T), max = ceiling(max(CDRSB.bl, na.rm = T)), 
        value = c(min(CDRSB.bl, na.rm = T), round(max(CDRSB.bl, na.rm = T),1)))
    })
    
    output$selectLogMemMMRM<-renderUI({
      LDELTOTAL.bl<-adnimerge$LDELTOTAL.bl
      LDELTOTAL.bl<-unique(as.numeric(LDELTOTAL.bl))
      sliderInput("logmembblMMRM", label = "Logical Memory", min = min(LDELTOTAL.bl, na.rm = T), max = ceiling(max(LDELTOTAL.bl, na.rm = T)), 
        value = c(min(LDELTOTAL.bl, na.rm = T), round(max(LDELTOTAL.bl, na.rm = T),1)))
    })
    
    output$selectDurationMMRM<-renderUI({
      sliderInput("studyDurationMMRM", label = "Select study duration (years)", min = 1, 
        max = 5, value = 1.5, step=0.5)
    })
    
    ## Grayout criteria selection ----
    observeEvent(input$submitCriteriaMMRM, {
      if(c("AGE") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("ageMMRM")
      }else{
        shinyjs::disable("ageMMRM")
      }
      if(c("ABETA.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("abetablMMRM")
      }else{
        shinyjs::disable("abetablMMRM")
      }
      
      if(c("APOE4") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("APOE4MMRM")
      }else{
        shinyjs::disable("APOE4MMRM")
      }
      
      if(c("TAU.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("taublMMRM")
      }else{
        shinyjs::disable("taublMMRM")
      }
      
      if(c("DX.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("dxMMRM")
      }else{
        shinyjs::disable("dxMMRM")
      }
      
      if(c("TAU/ABETA") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("tauRatioMMRM")
      }else{
        shinyjs::disable("tauRatioMMRM")
      }
      
      if(c("AV45.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("av45blMMRM")
      }else{
        shinyjs::disable("av45blMMRM")
      }
      
      if(c("LDELTOTAL.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("logmembblMMRM")
      }else{
        shinyjs::disable("logmembblMMRM")
      }
      
      if(c("MMSE.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("mmseblMMRM")
      }else{
        shinyjs::disable("mmseblMMRM")
      }
      
      
      if(c("CDRSB.bl") %in% input$criteriaMMRM | is.null(input$criteriaMMRM)){
        shinyjs::enable("cdrsbblMMRM")
      }else{
        shinyjs::disable("cdrsbblMMRM")
      }
      
    })
    
    observeEvent(input$methodADNIMMRM, {
      if(input$methodADNIMMRM %in% c("edland")){
        shinyjs::enable("edlandAllocationADNIMMRM")
      }else{
        shinyjs::disable("edlandAllocationADNIMMRM")
      }
    })
    
    ## Output (MMRM) ----
    ### MMRM data -----
    subdataMMRM<-reactive({
      a<-c("","AGE","ABETA.bl","APOE4","TAU.bl","DX.bl","TAU/ABETA","AV45.bl","LDELTOTAL.bl",
        "MMSE.bl","CDRSB.bl")
      if(identical(as.character(input$criteriaMMRM), a)){
        if("Positive" %in% as.character(input$APOE4MMRM) & "Negative" %in% as.character(input$APOE4MMRM)){
          apoe<-c(0,1,2)
        }else if("Positive" %in% as.character(input$APOE4MMRM) & !"Negative" %in% as.character(input$APOE4MMRM)){
          apoe<-c(1,2)
        }else if(!"Positive" %in% as.character(input$APOE4MMRM) & "Negative" %in% as.character(input$APOE4MMRM)){
          apoe<-c(0)
        }
        adni<-adnimerge %>%
          subset(AGE >= input$ageMMRM[1] & AGE <= input$ageMMRM[2] )%>%
          subset(ABETA.bl >= input$abetablMMRM[1] & ABETA.bl <= input$abetablMMRM[2]| ABETA.bl %in% NA_real_) %>%
          subset(APOE4 %in% c(0,1,2) |APOE4 %in% NA_integer_)%>%
          subset(TAU.bl >= input$taublMMRM[1] & TAU.bl <= input$taublMMRM[2]| TAU.bl %in% NA_real_) %>%
          subset(DX.bl %in% as.character(input$dxMMRM) | DX.bl %in% NA_character_) %>%
          subset(round(TAU.bl/ABETA.bl,1) >= input$tauRatioMMRM[1] & round(TAU.bl/ABETA.bl,1) <= input$tauRatioMMRM[2]| 
              round(TAU.bl/ABETA.bl,1) %in% NA_real_) %>%
          subset(AV45.bl >= input$av45blMMRM[1] & AV45.bl <= input$av45blMMRM[2]| AV45.bl %in% NA_real_) %>%
          subset(LDELTOTAL.bl >= input$logmembblMMRM[1] & LDELTOTAL.bl <= input$logmembblMMRM[2]| LDELTOTAL.bl %in% NA_real_) %>%
          subset(MMSE.bl >= input$mmseblMMRM[1] & MMSE.bl <= input$mmseblMMRM[2]| MMSE.bl %in% NA_real_) %>%
          subset(CDRSB.bl >= input$cdrsbblMMRM[1] & CDRSB.bl <= input$cdrsbblMMRM[2]| CDRSB.bl %in% NA_real_) %>%
          subset(Years.bl <= input$studyDurationMMRM + 0.25 | Years.bl %in% NA_real_) %>%
          filter(M <=  input$studyDurationMMRM*12) %>%
          filter(APOE4 %in% apoe| APOE4 %in% NA_integer_)
      }else {
        adni<-adnimerge %>%
          filter(Years.bl <= input$studyDurationMMRM + 0.25 | Years.bl %in% NA_real_) %>%
          filter(M <=  input$studyDurationMMRM*12)
        if("AGE" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(AGE >= input$ageMMRM[1] & AGE <= input$ageMMRM[2]| AGE %in% NA_real_)
        }
        if("ABETA.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(ABETA.bl >= input$abetablMMRM[1] & ABETA.bl <= input$abetablMMRM[2]| ABETA.bl %in% NA_real_)
        }
        if("APOE4" %in% as.character(input$criteriaMMRM)){
          if("Positive" %in% as.character(input$APOE4MMRM) & "Negative" %in% as.character(input$APOE4MMRM)){
            apoe<-c(0,1,2)
          }else if("Positive" %in% as.character(input$APOE4MMRM) & !"Negative" %in% as.character(input$APOE4MMRM)){
            apoe<-c(1,2)
          }else if(!"Positive" %in% as.character(input$APOE4MMRM) & "Negative" %in% as.character(input$APOE4MMRM)){
            apoe<-c(0)
          }
          adni<-adni %>%
            filter(APOE4 %in% apoe| APOE4 %in% NA_integer_)
        }
        if("TAU.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(TAU.bl >= input$taublMMRM[1] & TAU.bl <= input$taublMMRM[2]| TAU.bl %in% NA_real_) 
        }
        if("DX.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(DX.bl %in% as.character(input$dxMMRM)| DX.bl %in% NA_character_)
        }
        if("TAU/ABETA" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(round(TAU.bl/ABETA.bl,1) >= input$tauRatioMMRM[1] & round(TAU.bl/ABETA.bl,1) <= input$tauRatioMMRM[2]| 
                round(TAU.bl/ABETA.bl,1) %in% NA_real_)
        }
        if("AV45.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(AV45.bl >= input$av45blMMRM[1] & AV45.bl <= input$av45blMMRM[2]| AV45.bl %in% NA_real_)
        }
        if("LDELTOTAL.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(LDELTOTAL.bl >= input$logmembblMMRM[1] & LDELTOTAL.bl <= input$logmembblMMRM[2]| 
                LDELTOTAL.bl %in% NA_real_)
        }
        if("MMSE.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(MMSE.bl >= input$mmseblMMRM[1] & MMSE.bl <= input$mmseblMMRM[2]| MMSE.bl %in% NA_real_)
        }
        if("CDRSB.bl" %in% as.character(input$criteriaMMRM)){
          adni<-adni %>%
            filter(CDRSB.bl >= input$cdrsbblMMRM[1] & CDRSB.bl <= input$cdrsbblMMRM[2]| CDRSB.bl %in% NA_real_)
        }
      }

      adni$VISCODE[adni$VISCODE=="bl"]<-"bl0"
      adni$VISCODE2<-as.numeric(gsub("[^0-9.]", "",  adni$VISCODE))
      
      adni$Y <- adni[, input$longoutMMRM]
      
      adni <- left_join(
        adni %>% filter(M==0) %>% select(RID, Y),
        adni %>% filter(M>0),
        suffix = c('.bl', ''), 
        by='RID'
      )
      
      adni<-adni %>%
        mutate(Y.ch = Y-Y.bl) %>%
        filter(!is.na(Y.ch)) %>%
        mutate(
          Yr = as.factor(M/12),
          t.index=as.numeric(Yr),
          RID=factor(RID))
    })
    
    ### Baseline summaries ----
    observeEvent(input$summarybyMMRM,{
      output$baselineSummaryMMRM<-renderTable({
        dat<-subdataMMRM() %>%
          arrange(Years.bl)
        a<-as.character(input$criteriaMMRM)[-1]
        a[a=="TAU/ABETA"]<-"TAU_ABETARatio"

        if("TAU/ABETA" %in% input$criteriaMMRM){
          dat$TAU_ABETARatio<-with(dat, TAU.bl/ABETA.bl)
        }
        if("APOE4" %in% a){
          dat <-dat %>%
            mutate(APOE4=if_else(APOE4 %in% c(1,2), "Positive",
              if_else(APOE4 %in% c(0), "Negative", NA_character_)))
        }
        
        dat<-dat[!duplicated(dat$RID),]
        label(dat$AV45.bl)<-"AV45.bl"
        
        cov<-reactive({
          cov<-paste(a, collapse="+")
          cov<-paste0("~",cov)
        })

        regFormulaMMRM<- reactive({
          as.formula(paste(input$summarybyMMRM,cov()))
        })
        
        tbl<-arsenal::tableby(regFormulaMMRM(),data=dat,
          numeric.stats = c("N","Nmiss","range","meansd", "medianq1q3"),
          cat.stats = c("N","Nmiss", "countpct"),
          digits=1)
        as.data.frame(summary(tbl, text = "html"))
      }, sanitize.text.function = function(x) x)
    })
    ### model fitting ----
    modelFitMMRM<-reactive({
      if(input$matrixADNIMMRM %in% c("exchangeable")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),
            na.action = na.omit,
            correlation = corCompSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+")
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),
            na.action = na.omit,
            correlation = corCompSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
        
      }else if(input$matrixADNIMMRM %in% c("general")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),
            na.action = na.omit,
            correlation = corSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+") 
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),
            na.action = na.omit,
            correlation = corSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
      }else if(input$matrixADNIMMRM %in% c("ar1")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),na.action = na.omit,
            correlation = corAR1(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+")
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),na.action = na.omit,
            correlation = corAR1(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
      }
    })
    
    changeMMRM<-reactive({
      summary(modelFitMMRM())$coef[paste0('Yr', as.character(input$studyDurationMMRM))]
    })

    sigmaaMMRM<-reactive({
      WEIGHTS <- coef(modelFitMMRM()$modelStruct$varStruct, unconstrained = FALSE)
      modelFitMMRM()$sigma  * last(WEIGHTS)
    })
    
    output$modelFitSummaryMMRM<-renderPrint({
      print(summary(modelFitMMRM()))
    })
    
    output$indPlotMMRM<-renderPlotly({
      p <- ggplot(data = subdataMMRM(), aes(x = Years.bl, 
        y = Y.ch, group = RID))+
        geom_line(alpha=0.25)+theme_bw()+xlab("Years from baseline")+
        ylab(paste(as.character(input$longoutMMRM), 'change from baseline'))
      ggplotly(p)
    })
    
    output$plotProfileMMRM<-renderPlotly({
      p<-ggplot(subdataMMRM(), aes(x=Years.bl,y=Y.ch)) +
        geom_point(alpha=0.25) +
        geom_smooth()+theme_bw()+xlab("Years from baseline")+
        ylab(paste(as.character(input$longoutMMRM), 'change from baseline'))
      ggplotly(p)
    })
    
    output$digglePlotMMRM<-renderPlotly({
      if(input$matrixADNIMMRM %in% c("exchangeable")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),na.action = na.omit,
            correlation = corCompSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+")
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),na.action = na.omit,
            correlation = corCompSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
        
      }else if(input$matrixADNIMMRM %in% c("general")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),na.action = na.omit,
            correlation = corSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+") 
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),na.action = na.omit,
            correlation = corSymm(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
      }else if(input$matrixADNIMMRM %in% c("ar1")){
        if(is.null(input$selectCovariatesMMRM)){
          m1 <- gls(Y.ch ~ Y.bl + Yr, 
            subdataMMRM(),
            na.action = na.omit,
            correlation = corAR1(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }else{
          cov<-paste(input$selectCovariatesMMRM,collapse = "+")
          m1 <- gls(formula(paste("Y.ch","~",cov,"+ Y.bl + Yr")), 
            subdataMMRM(),na.action = na.omit,
            correlation = corAR1(form = ~ t.index | RID),
            weights = varIdent(form = ~ 1 | Yr))
        }
      }
      
      C <- corMatrix(m1$modelStruct$corStruct)[[1]]
      ra <- seq(1,input$percRetentionA/100,length=nrow(C))
      rb <- seq(1,input$percRetentionB/100,length=nrow(C))
      power<- seq(input$powerADNIMMRM[1],input$powerADNIMMRM[2],1)
      ssize<-rep(0,length(power))
      delta<-changeMMRM()*input$pchangeADNIMMRM/100
      for(i in 1:length(power)){
        fit<-power.mmrm(Ra = C, ra = ra, sigmaa = sigmaaMMRM(), rb=rb,
          sig.level=input$alphaADNIMMRM,alternative=input$alternativeADNIMMRM, 
          power=power[i]/100, delta=delta, lambda=input$edlandAllocationADNIMMRM)
        ssize[i]<-fit$n1+fit$n2
      }
      
      dat<-data.frame(n=ssize, power=power)
      plot_ly(data=dat[-c(1:5,nrow(dat)),], x=~n, y=~power, type = 'scatter', mode = 'lines') %>%
        layout(title="Power analysis for MMRM",
          xaxis=list(title="Total sample size"), yaxis=list(title="Power (%)"))
      #} 
    })
    
    output$describeMethodADNIMMRM<-renderText({
      if(input$matrixADNIMMRM %in% c("general")){
        "Sample size/power calculation for a mixed model of repeated measures with general correlation structure.
        See Lu, Luo, & Chen (2008) for parameter definitions and other details"
      }else if(input$matrixADNIMMRM %in% c("exchangeable")){
        "Sample size/power calculation for a mixed model of repeated measures with compound symmetric correlation structure.
        See Lu, Luo, & Chen (2008) for parameter definitions and other details."
      }else if(input$matrixADNIMMRM %in% c("ar1")){
        "Sample size/power calculation for a mixed model of repeated measures with ar1 correlation structure.
        See Lu, Luo, & Chen (2008) for parameter definitions and other details"
      }
    })
    
    output$summarySelectionADNIMMRM<-renderTable({
      a<-data.frame(Summary=c("Var-cov structure", "Type of test", "Type I error", 
        "Residual standard deviation at last visit",
        "Pilot estimate of placebo group change",
        "Effect size (raw scale)", "Effect size (% of placebo change)", "Allocation ratio", "Percent retention in group a",
        "Percent retention in group b", "Observation times (years)"),
        Value=as.character(c(varCov[input$matrixADNIMMRM],
          input$alternativeADNIMMRM, 
          input$alphaADNIMMRM,
          round(sigmaaMMRM(), 4),
          round(changeMMRM(), 4),
          round(changeMMRM()*input$pchangeADNIMMRM/100, 2),
          paste0(input$pchangeADNIMMRM,'%'), 
          input$edlandAllocationADNIMMRM,
          paste0(input$percRetentionA,'%'), 
          paste0(input$percRetentionB,'%'),
          paste0(levels(subdataMMRM()$Yr), collapse = ", ")
          )))
    })
    
    # About ----
    output$aboutApp = renderText({
      HTML(paste(tags$div(
        tags$b(h3("Sample Size and Power Analysis Dashboard")),
        tags$p("Sample size calculation or power analysis is an important 
               requirement when designing a new trial or study. For many 
               non-statisticians, the processes involved in performing these 
               types of analyses can be daunting. A major hurdle to overcome is 
               the availability of pilot parameters which are required inputs 
               for generating sample size and power outputs. "
        ), 
        tags$p("This ShinyApp dashboard is developed to easily generate sample
               size and conduct power analysis for a longitudinal study design
               with two-group comparisons for a continuous outcome. The App
               implements the sample size formulae of Liu and Liang (1997) and
               Diggle et al (1994) using functions developed in the R",tags$code("longpower"),
          "package. The",tags$code('longpower'),"package handles cases where time is
               treated either as continuous or categorical. The former approach
               uses the linear mixed model with random intercept and slope while
               the later leads to the well-known Mixed Model of Repeated Measures
               (MMRM) used in many clinical trial applications for conducting
               statistical analysis. "),
        tags$p("The dashboard is in two parts. The first part accepts user inputs 
               to generate sample size when time is treated as both categorical 
               and continuous. Thus, this part assumes that the user already has 
               pilot parameter estimates including effect size and known variance. 
               Users can generate similar sample sizes and perform power analysis 
               using different sample size methods (diggle, liuliang, and edland). 
               The second part of the dashboard uses the methodology that is similar 
               to the first part except that Alzheimer's Disease Neuroimaging Initiative 
               (ADNI)-based pilot parameters are generated based on user-selected
               inclusion and exclusion criteria, primary outcome, duration of the study, 
               and covariate options to be included in the linear mixed and MMRM models. 
               Some descriptive summaries and plots of the selected sample used for performing
               the analysis are also dynamically generated. "),
        tags$b(h3("Data")),
        tags$p("Data used as inputs in the ADNI-based generator of this dashboard were obtained from the Alzheimer's Disease Neuroimaging Initiative
(ADNI)", tags$a(href="adni.loni.usc.edu","database"), paste0("on ", adniDate, ". As such, the investigators within the ADNI contributed to the design
and implementation of ADNI and/or provided data but did not participate in the development of this dashboard.
A complete listing of ADNI investigators can be found at:"), 
          tags$a(href="http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf",
            "Acknowledgement List"),". The ADNI was launched in 2003 as a public-private
partnership, led by Principal Investigator Michael W. Weiner, MD. The primary goal of ADNI has been to
test whether serial magnetic resonance imaging (MRI), positron emission tomography (PET), other
biological markers, and clinical and neuropsychological assessment can be combined to measure the
progression of mild cognitive impairment (MCI) and early Alzheimer's disease (AD)"),
        tags$b(h3("Packages")),
        tags$p(
          tags$ul(
            tags$li("Sample size/power calculations - the", 
              tags$a(href="https://cran.r-project.org/web/packages/longpower/longpower.pdf", 
                "longpower"),"package"), 
            tags$li("The ADNI data package - the", 
              tags$a(href="https://adni.bitbucket.io/", 
                "ADNIMERGE"),"package"),
            tags$li("Dashboard interface - the ", 
              tags$a(href="https://rstudio.github.io/shinydashboard/", "shinydashboard"),"package"), 
            tags$li("Visualization - the ", 
              tags$a(href="https://plotly.com/r/", "plotly"),"  and ", 
              tags$a(href="https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf", "ggplot")," packages for the plots"), 
            tags$li("Data manipulation - ",tags$a(href="https://cran.r-project.org/web/packages/dplyr/dplyr.pdf", "dplyr"),"and ", 
              tags$a(href="https://www.tidyverse.org/packages/", "tidyverse"))
          )
        ),
        tags$b(h3("Authors")),
        tags$p("This shinydashboard is developed by",
          tags$ul(
            tags$li(tags$strong("Samuel Iddi"),"(PhD), Department of Statistics and Actuarial Science, University of Ghana"),
            tags$li(tags$strong("Michael C. Donohue"),"(PhD),  Alzheimer's Therapeutic Research Institute (ATRI), Keck School of Medicine, University of Southern California")
          ),
          tags$b(h3("Contacts")), 
          tags$p(
            "If you have any questions, comments, and/or suggestions, please email Samuel (siddi@aphrc.org) or Michael (mdonohue@usc.edu.)" 
          ),
          tags$b(h3("Last Update")), 
          tags$p(
            adniDate
          ),
          br()
        )
      )))
    })
    
    
    
    
    
    
  }
)