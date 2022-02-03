rm(list=ls())
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyalert)
library(plotly)
library(shinyjs)
library(shinyMatrix)
library(shinyWidgets)

loadingHelper <- function() {
  conditionalPanel("$('html').hasClass('shiny-busy')",
    id = "waitIndicator", "Running...please wait",
    HTML("<br><img src = 'ajax_loader.gif' alt='loading' >"))
  
}

shinyUI(dashboardPage(skin = 'blue',
  dashboardHeader(title="Power and Sample Size for Longitudinal Study Designs",
    #"Power and Sample Size for Longitudinal Study Designs", 
    titleWidth = 450),
  dashboardSidebar(
    shinyjs::useShinyjs(),
    # select app section ----
    sidebarMenu(id="design",
      menuItem(HTML("Power Analysis Based on  <br/> Linear Mixed Model (LMM)"),
        menuSubItem(text = HTML("Sample size calculations for <br/> linear mixed models of<br/>  rate of change"),
          tabName = "lmmpower"),
        menuSubItem(text = HTML("Mixed Model Repeated Measures <br/> (MMRM) sample size calculations"),
          tabName = "powermmrm")),
      menuItem(text=HTML('Power Analysis with ADNI-based <br/> pilot estimate generator'),
        menuSubItem(text = HTML("Sample size calculations for <br/> linear mixed models of<br/>  rate of change"), 
          tabName = "lmmpowerADNI"),
        menuSubItem(text = HTML("Mixed Model Repeated Measures <br/> (MMRM) sample size calculations"), 
          tabName = "powermmrmADNI")),
      menuItem("About", tabName = 'about')
      
    )
  ),
  dashboardBody(
    tags$style(type="text/css",
      ".shiny-output-error { visibility: hidden; }",
      ".shiny-output-error:before { visibility: hidden; }"
    ),
    # lmm ----
    tabItems(
      tabItem(tabName = "lmmpower",
        fluidRow(
          box(status = 'danger', title='Inputs', solidHeader = T,width = 12,
            fluidRow(column(3,
              selectInput(inputId = "analysisType",label = "Analysis type", choices = list("Sample size", "Power"),
                selected = "sample size")),
              column(3, 
                #numericInput(inputId="sampleSize", label = "Sample size (n)", value = 100, min=30)),
                sliderInput(inputId="sampleSize", label = "Sample size (n)", min=0, max=5000, value = c(100,1000))),
              column(3,
                sliderInput("power",label="Power (%)", min=0, max=100, step = 1, value = c(50,95))),
              column(3,
                numericInput(inputId = "startTime",label="Start time (years)", value = 0, min = 0))
            ),
            fluidRow(
              column(3,
                numericInput(inputId = "entTime",label='End time (years)', value = 1.5, min=0)),
              column(3, 
                numericInput(inputId = "timeStep",label="Time step (years)",value = 0.5, min=0)),
              column(3,
                sliderInput("alpha",label="Type I error (sig.level)", min=0, max=1, step = 0.01, value = 0.05)),
              column(3,
                radioButtons("alternative", label="Type of test", choices = list("two.sided","one.sided"), selected = "two.sided"))
            ),
            fluidRow(
              column(3,
                radioButtons(inputId="estimate", 
                  label = "Effect size scale", 
                  choiceNames = c('Raw scale of the outcome', 'Percent of placebo change'),
                  choiceValues = c("delta","percent change"),
                  selected = "delta")),
              column(3,
                numericInput("beta",
                  label="Pilot estimate of the placebo change (beta)", 
                  value=0.5, min = 0, step = 0.5)),
              column(3,
                sliderInput(inputId="pct.change",label="Percent change in the pilot estimate of the parameter of interest (pct.change)",
                  min=0, max = 1, step = 0.01, value = 0.30)),
              column(3,
                numericInput("delta","Target treatment effect size (delta)",min = 0, step = 0.01, value = 1.5)),
              # column(3,
              #    sliderInput("beta.CI",label="95% confidence limits of the pilot estimate of beta (beta.CI)", value=95, min=0, max=100, step=5)),
              # column(3,
              #   sliderInput("delta.CI",label="95% confidence limits of the effect size (delta.CI)", value=95, min=0, max=100, step=5)),
              
            ),
            fluidRow(
              column(3,
                numericInput("sig2.i","Pilot estimate of variance of random intercept (sig2.i)", value=55,min = 0, step = 0.01)),
              column(3,
                numericInput("sig2.s","Pilot estimate of variance of random slope (sig2.s)", value=22,min = 0, step = 0.01)),
              column(3,
                numericInput("cov.s.i","Pilot estimate of covariance of random slope and intercept(cov.s.i)", value=29,min = 0, step = 0.01)),
              column(3,
                numericInput("sig2.e","Pilot estimate of residual variance (sig2.e)", value=10,min = 0, step = 0.01))
              
            ),
            fluidRow(
              column(4,
                radioButtons(inputId="method", label = "Method", 
                  choiceNames = c("Liu and Liang (1997)", "Diggle et al (2002)", "Edland (2009)"),
                  choiceValues = c("diggle","liuliang","edland"), selected = "diggle")),
              # column(3,
              #        radioButtons(inputId="matrix", label = "Association structure", choices = list("covariance"), selected = "covariance")),
              
              # column(2,tags$div(id="rmoo"
              #   #tags$b("Pilot estimate of a marginal model working correlation matrix (R)")
              # )
              #        #tags$b("Pilot estimate of a marginal model working correlation matrix (R)"),
              #        # matrixInput("Rmatrix",value =uiOutput("R") , rows = list(names = T),
              #        #             cols = list(names = T),copy = TRUE,paste = TRUE, class = 'numeric')
              #        # uiOutput("R")
              #        ),
              column(4, sliderInput("edlandAllocation", "Allocation ratio (lambda)", min = 0, max = 5, value = 1))
            )
          )),
        fluidRow(
          box(status = 'success', title='Outputs', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("diggleSS")),column(6, textOutput("describeMethod"),br(),
                tableOutput("summarySelection")),
              
            )    
          ))
      ),
      # mmrm ----
      tabItem(tabName = "powermmrm",
        fluidRow(
          box(status = 'danger', title='Inputs', solidHeader = T,width = 12,
            fluidRow(column(3,
              selectInput(inputId = "analysisTypeMMRM",label = "Analysis type", choices = list("Sample size", "Power"),
                selected = "sample size")),
              column(3, 
                #numericInput(inputId="sampleSize", label = "Sample size (n)", value = 100, min=30)),
                sliderInput(inputId="sampleSizeMMRM", label = "Total sample size (N)", min=0, max=2000, value = c(50,100))),
              column(3,
                sliderInput(inputId = "powerMMRM",label="Power (%)", min=0, max=100, step = 1, value = c(50,95))),
              column(3,
                sliderInput(inputId = "timePoints",label="Number of time pionts", min = 1, max = 50, value = 4))
            ),
            fluidRow(
              column(3,
                numericInput(inputId ="sigmaaMMRM",label="Standard deviation of observation of interest in group a (sigmaa)", value=1,min = 0, step = 0.01)),
              column(3,
                numericInput(inputId ="sigmabMMRM","Standard deviation of observation of interest in group b. (sigmab)", value=1,min = 0, step = 0.01)),
              column(3,
                sliderInput(inputId ="alphaMMRM",label="Type I error (sig.level)", min=0, max=1, step = 0.01, value = 0.05)),
              column(3,
                radioButtons(inputId ="alternativeMMRM", label="Type of test", choices = list("two.sided","one.sided"), selected = "two.sided"))
            ),
            fluidRow(
              column(3,
                radioButtons(inputId="matrixMMRM", label = "Select association structure", 
                  choices = list("exchangeable","general","ar1"), selected = "exchangeable")),
              #column(3,actionButton("updateMMRM","Update/Enter correlation structure"), icon = icon("refresh")),
              column(3,
                numericInput("rhoMMRM",label="Exchangeable correlation (rho)", value=0.25,min = 0, step = 0.05)),
              column(3,
                radioButtons(inputId="estimateMMRM", 
                  label = "Effect size scale", 
                  choiceNames = c('Raw scale of the outcome', 'Percent of placebo change'),
                  choiceValues = c("delta","percent change"), 
                  selected = "delta")),
              column(3,
                numericInput("betaMMRM", label="Estimate of the placebo change (beta)", value=0.5,min = 0, step = 0.5))
            ),
            fluidRow(
              column(4,
                numericInput(inputId ="deltaMMRM","Effect size (delta)",min = 0, step = 0.01, value = 0.5)),
              column(4,
                numericInput(inputId ="lambdaMMRM","Allocation ratio (lambda)",min = 1, step = 0.5, value = 2)),
              column(4,
                sliderInput(inputId="pct.changeMMRM",label="Percent change in the pilot estimate of the parameter of interest (pct.change)",
                  min=0, max = 1, step = 0.01, value = 0.30))
              #,column(6,tags$div(id="lambdaMMRM"))
              
            ),
            fluidRow(
              column(2,tags$div(id="raMMRM")),
              column(2,tags$div(id="rbMMRM")),
              column(3,tags$div(id="RaMMRM")),
              column(3,tags$div(id="RbMMRM"))
              ,column(2,actionButton("updateMMRM","Update/Enter"),
                p("Use this action button to update size of  correlation matrix"))
              #submitButton("Update!")
              
              
            )
          )),
        fluidRow(
          
          box(status = 'success', title='Outputs', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("mmrmplot")),column(6, textOutput("describeMethodMMRM"), br(),
                tableOutput("summarySelectionMMRM"))
            )    
          ))
      ), #tabItem2 end
      # lmm ADNI ----
      tabItem(tabName = "lmmpowerADNI",
        ## baseline selections ----
        fluidRow(
          box(status = 'danger', title='Baseline selections', solidHeader = T,width = 12,
            fluidRow(
              column(12, multiInput("criteria", "Inc/Exc criteria", 
                choiceNames = c("","Age","CSF Aβ42",
                  "APOE4","CSF Tau","Diagnosis",
                  "CSF Tau/Aβ42", "AV45 (SUVR)", "Logical Memory", "MMSE", "CDRSB"),
                choiceValues = c("", "AGE", "ABETA.bl",
                "APOE4", "TAU.bl", "DX.bl",
                "TAU/ABETA", "AV45.bl", "LDELTOTAL.bl", "MMSE.bl", "CDRSB.bl"),
                selected = c("", "AGE", "ABETA.bl",
                  "APOE4", "TAU.bl", "DX.bl",
                  "TAU/ABETA", "AV45.bl", "LDELTOTAL.bl", "MMSE.bl", "CDRSB.bl"),
                options = list(
                  enable_search = T,
                  non_selected_header = "Choose to select:",
                  selected_header = "You have selected:"
                )),
                actionButton("submitCriteria","Submit selected criteria"))
            )  
          )
        ),
        ## Inclusion/Exclusion criteria ----
        fluidRow(
          box(status = 'danger', title='Inclusion/Exclusion criteria', solidHeader = T,width = 12,
            fluidRow(
              column(3, uiOutput("selectAge")), column(3, checkboxGroupInput('APOE4', HTML(paste0('APOE4 status (', 
                tags$em("Do not deselect both"),")")), 
                choices = c("Positive", "Negative"), selected =c("Positive", "Negative"))),
              column(3, uiOutput("selectAV45")),
              column(3, uiOutput("selectDuration"))
            ),
            fluidRow(
              column(4, uiOutput("selectAbeta")),column(4, uiOutput("selectTau")), column(4, uiOutput("selectRatio"))
            ),
            fluidRow(
              column(3, uiOutput("selectDx")), #helpText("You should not deselect all")), 
              column(3, uiOutput("selectMMSE")), column(3, uiOutput("selectCDRSB")),
              column(3, uiOutput("selectLogMem"))
            )
          ),
          ## baseline summary ----
          box(status = 'success', title='Baseline summary', solidHeader = T,width = 12,
            fluidRow(
              column(6, radioButtons(inputId = "summaryby", "Sumarize by", 
                choices = c("Gender"="PTGENDER","Education"="PTEDUCAT",
                  "Ethnicity"="PTETHCAT","Race"="PTRACCAT"),selected =c("Gender"="PTGENDER"),
                inline = T))
            ),
            fluidRow(
              column(12, tableOutput("baselineSummary"))
            )
          ),
          ## covariate options ----
          box(status = 'danger', title='Covariate Options', solidHeader = T,width = 5,
            fluidRow(
              column(12,checkboxGroupInput('selectCovariates', 'Select covariate', 
                choices = c("Age"="AGE","Gender"="PTGENDER","Education"="PTEDUCAT",
                  "Ethnicity"="PTETHCAT","Race"="PTRACCAT"), 
                selected = c(""), inline = T))
            )
          ),
          ## outcome options ----
          box(status = 'danger', title='Outcome Options', solidHeader = T,width = 7,
            fluidRow(
              column(12, radioButtons("longout",label = "Longitudinal Outcomes", 
                choices = c("FDG PET"="FDG","AV45 PET"="AV45","ADAS11"="ADAS11","ADAS13"="ADAS13","MMSE"="MMSE",
                  "CDRSB"="CDRSB","mPACCtrailsB"="mPACCtrailsB"), selected = c("CDRSB"="CDRSB"), inline=T))
              
            )
          ),
          ## Plots of pilot data ----
          box(status = 'success', title='Plots of pilot data', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("indPlot")),
              column(6, plotlyOutput("plotProfile"))
            )
          ),
          ## Other design options ----
          box(status = 'danger', title='Other design options', solidHeader = T,width = 12,
            fluidRow(
              column(4, radioButtons(inputId="methodADNI", label = "Method", 
                choiceNames = c("Liu and Liang (1997)", "Diggle et al (2002)", "Edland (2009)"),
                choiceValues = c("diggle","liuliang","edland"), 
                selected = "diggle", inline =T)),
              column(3, radioButtons(inputId ="alternativeADNI", label="Type of test", 
                choices = list("two.sided","one.sided"), selected = "two.sided", inline = T)),
              column(5,
                sliderInput("powerADNI",label="Power (%)", min=0, max=100, step = 1, value = c(50,95)))
            ),
            fluidRow(
              column(4,  sliderInput(inputId ="alphaADNI",label="Type I error (sig.level)",
                min=0, max=1, step = 0.01, value = 0.05)),
              # column(3, radioButtons("directionChange",label = "Direction of change in slope",choices = c("Decrease","Increase"), 
              #                        selected = "Increase", inline = T)),
              column(4, sliderInput(inputId ="pchangeADNI",label="Effect size (% of estimated placebo slope)",
                min=0, max=100, value = 75)),
              column(4, sliderInput("edlandAllocationADNI", "Allocation ratio (lambda)", min = 0, max = 5, value = 1))
            ),
            fluidRow(
              column(4,
                sliderInput("timeStepADNI", label="Time step (years)", min=0.5, max=2, step = 0.5, value = 0.5))
            )
          ),
          ## Power analysis results ----
          box(status = 'success', title='Power analysis results', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("digglePlot")), 
              column(6, uiOutput("describeMethodADNI"), 
                tableOutput("summarySelectionADNI"))
            )
          ),
          ## Model fit to pilot data ----
          box(status = 'success', title='Model fit to pilot data', solidHeader = T,width = 12,
            fluidRow(
              column(12, verbatimTextOutput("modelFitSummary"))
            )
          )
        )
      ), #tabitem3 ends
      # mmrm ADNI ----
      tabItem(tabName = "powermmrmADNI",
        ## Baseline selections ----
        fluidRow(
          box(status = 'danger', title='Baseline selections', solidHeader = T,width = 12,
            fluidRow(
              column(12, multiInput("criteriaMMRM", "Inc/Exc criteria", 
                choiceNames = c("","Age","CSF Aβ42",
                  "APOE4","CSF Tau","Diagnosis",
                  "CSF Tau/Aβ42", "AV45 (SUVR)", "Logical Memory", "MMSE", "CDRSB"),
                choiceValues = c("", "AGE", "ABETA.bl",
                  "APOE4", "TAU.bl", "DX.bl",
                  "TAU/ABETA", "AV45.bl", "LDELTOTAL.bl", "MMSE.bl", "CDRSB.bl"),
                selected = c("", "AGE", "ABETA.bl",
                  "APOE4", "TAU.bl", "DX.bl",
                  "TAU/ABETA", "AV45.bl", "LDELTOTAL.bl", "MMSE.bl", "CDRSB.bl"),
                options = list(
                  enable_search = T,
                  non_selected_header = "Choose to select:",
                  selected_header = "You have selected:"
                )),
                actionButton("submitCriteriaMMRM","Submit selected criteria"))
            )  
          )
        ),
        fluidRow(
          ## Inclusion/Exclusion criteria ----
          box(status = 'danger', title='Inclusion/Exclusion criteria', solidHeader = T,width = 12,
            fluidRow(
              column(3, uiOutput("selectAgeMMRM")), column(3, checkboxGroupInput('APOE4MMRM', label = HTML(paste0('APOE4 status (', 
                tags$em("Do not deselect both"),")")), 
                choices = c("Positive", "Negative"), selected =c("Positive", "Negative"))),
              #             helpText("You should not deselect both")),
              column(3, uiOutput("selectAV45MMRM")),
              column(3, uiOutput("selectDurationMMRM")),
              
            ),
            fluidRow(
              column(4, uiOutput("selectAbetaMMRM")),column(4, uiOutput("selectTauMMRM")), column(4, uiOutput("selectRatioMMRM"))
            ),
            fluidRow(
              column(3, uiOutput("selectDxMMRM")), #helpText("You should not deselect all")), 
              column(3, uiOutput("selectMMSEMMRM")), 
              column(3, uiOutput("selectCDRSBMMRM")),
              column(3, uiOutput("selectLogMemMMRM"))
            )
          ),
          ## Baseline summary ----
          box(status = 'success', title='Baseline summary', solidHeader = T,width = 12,
            fluidRow(
              column(6, radioButtons(inputId = "summarybyMMRM", "Sumarize by", 
                choices = c("Gender"="PTGENDER","Education"="PTEDUCAT",
                  "Ethnicity"="PTETHCAT","Race"="PTRACCAT"),
                selected =c("Gender"="PTGENDER"),
                inline = T))
            ),
            fluidRow(
              column(12, tableOutput("baselineSummaryMMRM"))
            )
          ),
          ## Covariate Options ----
          box(status = 'danger', title='Covariate Options', solidHeader = T,width = 5,
            fluidRow(
              # column(4, uiOutput("selectGender"),uiOutput("selectEduc")), #column(3, uiOutput("selectEduc")), 
              # column(4, uiOutput("selectEth"),checkboxGroupInput('APOE42', 'APOE4 status', 
              #       choices = c("APOE+"=0, "APOE-"=1), selected = c("APOE+"=0, "APOE-"=1))),column(4, uiOutput("selectRace"))
              # 
              column(12,checkboxGroupInput('selectCovariatesMMRM', 'Select covariate', 
                choices = c("Age"="AGE","Gender"="PTGENDER","Education"="PTEDUCAT",
                  "Ethnicity"="PTETHCAT","Race"="PTRACCAT"), 
                selected = c(""), inline = T))
            )
          ),
          ## Outcome Options ----
          box(status = 'danger', title='Outcome Options', solidHeader = T,width = 7,
            fluidRow(
              column(12, radioButtons("longoutMMRM",label = "Longitudinal Outcomes", 
                choices = c("FDG PET"="FDG","AV45 PET"="AV45","ADAS11"="ADAS11","ADAS13"="ADAS13","MMSE"="MMSE",
                  "CDRSB"="CDRSB","mPACCtrailsB"="mPACCtrailsB"), selected = c("CDRSB"="CDRSB"), inline=T))
              
            )
          ),
          ## Plots of pilot data ----
          box(status = 'success', title='Plots of pilot data', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("indPlotMMRM")),
              column(6, plotlyOutput("plotProfileMMRM"))
            )
          ),
          ## Other design options ----
          box(status = 'danger', title='Other design options', solidHeader = T,width = 12,
            fluidRow(
              column(3,
                radioButtons(inputId="matrixADNIMMRM", label = "Select var-cov structure", 
                  choiceNames = c("Compound symmetric, heterogeneous", 
                    "Unstructured, heterogeneous",
                    "AR1, heterogeneous"),
                  choiceValues = c("exchangeable", "general", "ar1"), 
                  selected = "exchangeable", 
                  inline = T)),
              column(3, radioButtons(inputId ="alternativeADNIMMRM", label="Type of test", 
                choices = list("two.sided","one.sided"), selected = "two.sided", inline = T)),
              column(3,
                sliderInput("powerADNIMMRM",label="Power (%)", min=0, max=100, step = 1, value = c(50,95))),
              column(3,  sliderInput(inputId ="alphaADNIMMRM",label="Type I error (sig.level)",
                min=0, max=1, step = 0.01, value = 0.05))
            ),
            fluidRow(
              # column(3, radioButtons("directionChangeMMRM",label = "Direction of change in slope",choices = c("Decrease","Increase"), 
              #                        selected = "Increase", inline = T)),
              column(3, sliderInput(inputId ="pchangeADNIMMRM",label="Effect size (% of estimated placebo change)",
                min=0, max=100, value = 75)),
              column(3, sliderInput("edlandAllocationADNIMMRM", "Allocation ratio (lambda)", min = 0, max = 5, value = 1)),
              column(3, sliderInput(inputId ="percRetentionA",label="Percent retention in pilot group a (ra)",
                min=0, max=100, value = 80)),
              column(3, sliderInput(inputId ="percRetentionB",label="Percent retention in group b (rb)",
                min=0, max=100, value = 80))
            )
          ),
          ## Power analysis results ----
          box(status = 'success', title='Power analysis results', solidHeader = T,width = 12,
            fluidRow(
              column(6, plotlyOutput("digglePlotMMRM")), 
              column(6, uiOutput("describeMethodADNIMMRM"), 
              tableOutput("summarySelectionADNIMMRM"))
            )
          ),
          ## Model fit to pilot data ----
          box(status = 'success', title='Model fit to pilot data', solidHeader = T,width = 12,
            fluidRow(
              column(12, verbatimTextOutput("modelFitSummaryMMRM"))
            )
          )
        )
      ), #tabitem 4 ends
      tabItem(tabName = "about",
        fluidRow(
          column(12,
            htmlOutput("aboutApp")
          )
        )
      )#tabItem About Ends
    )#tabItems end
  ) #dashboard body end
)
)
