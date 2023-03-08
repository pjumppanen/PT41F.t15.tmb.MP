# Define UI ----
ui <- shiny::fluidPage(

  # App title 
  shiny::titlePanel("IOTC Big Eye Tuna Management Procedure"),
  
  # Sidebar layout
  shiny::sidebarLayout(
    # Sidebar panel for inputs
    shiny::sidebarPanel(
      shiny::fileInput("file", 
                shiny::h3("Catch, CPUE and TAC file"),
                accept = c("text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
      shiny::textOutput("dummy"),
      shiny::tags$head(shiny::tags$script('
                            var PageWidth  = 0;
                            var PageHeight = 0;
                            $(document).on("shiny:connected", function(e) {
                                PageWidth  = window.innerWidth;
                                PageHeight = window.innerHeight;
                                Shiny.setInputValue("PageWidth", PageWidth);
                                Shiny.setInputValue("PageHeight", PageHeight);
                            });
                            $(window).resize(function(e) {
                                PageWidth  = window.innerWidth;
                                PageHeight = window.innerHeight;
                                Shiny.setInputValue("PageWidth", PageWidth);
                                Shiny.setInputValue("PageHeight", PageHeight);
                            });
                        ')),

      width=3                           
    ),

    # Main panel for outputs
    shiny::mainPanel(
      shiny::tabsetPanel(
        shiny::tabPanel("Recommendation",
          shiny::htmlOutput("TAC"),
          shiny::plotOutput(outputId = "cobsPlot", inline=TRUE)
        ),
        shiny::tabPanel("Model Diagnostics",
          shiny::plotOutput(outputId = "cpuePlot", inline=TRUE)
        )
      ),
      width=9
    )
  )
)

# Define server logic ----
server <- function(input, output)
{
  check <- function(input)
  {
    if (is.null(input$file)          || 
        is.null(input$file$datapath) ||
        (nchar(input$file$datapath) <= 0))
    {
      return("please select a Catch and CPUE history .csv file")
    }

    CatchAndCPUE <- read.csv(input$file$datapath, header = TRUE)
    Names        <- names(CatchAndCPUE)

    # Check for Year column
    if (!any(Names == "Year"))
    {
      return("Year Column missing (year)")
    }

    # Check for Catch column
    if (!any(Names == "Catch"))
    {
      return("Catch Column missing (observed catch)")
    }

    # Check for CPUE column
    if (!any(Names == "CPUE"))
    {
      return("CPUE Column missing (observed CPUE)")
    }

    # Check for TAC column
    if (!any(Names == "TAC"))
    {
      return("TAC Column missing (historic TAC)")
    }

    # Check for ascending contiguous years
    MinYear <- min(CatchAndCPUE$Year)
    MaxYear <- max(CatchAndCPUE$Year)
    y       <- MinYear

    for (idx in 1:length(CatchAndCPUE$Year))
    {
      if (CatchAndCPUE$Year[idx] != y)
      {
        return("Year must be in contiguous ascending years")
      }

      y <- y + 1
    }

    return (NULL)
  }
  
  data <- shiny::reactive({
    shiny::validate(
      check(input)
    )

    CatchAndCPUE <- read.csv(input$file$datapath, header = TRUE)
    results      <- assessMP(input$file$datapath)

    return (list(CE=CatchAndCPUE, TAC=results$TAC, B=results$B, Depletion=results$Depletion, q=results$q, plots=results$plots))
  })

  graphWidthCatch <- shiny::reactive({
    input$PageWidth * 0.43 + 58
  })

  graphWidthCPUE <- shiny::reactive({
    input$PageWidth * 0.43
  })

  graphHeight <- shiny::reactive({
    input$PageHeight * 0.66
  })

  graphWidthDiag <- shiny::reactive({
    input$PageWidth * 0.5
  })

  graphHeightDiag <- shiny::reactive({
    input$PageHeight * 0.70
  })

  output$dummy <- shiny::reactive({
    if (!is.null(input$file))
    {
      Data <- data()

      output$cobsPlot <- shiny::renderPlot({
          data <- Data$CE
          data <- rbind(data, list(Year=max(data$Year) + 1, CPUE=NA, Catch=NA, TAC=NA))

          RecommendedTAC                          <- rep(NA, nrow(data))
          RecommendedTAC[length(RecommendedTAC)]  <- Data$TAC
          data                                    <- cbind(data, "Recommended TAC"=RecommendedTAC)

          data_melt <- reshape2::melt(data[, c("Year", "Catch", "TAC", "Recommended TAC")], id.vars='Year', value.name='Catch')
          colors    <- c("Catch"="#00345D", "TAC"="#00A9CE", "Recommended TAC"="#000080")
          shapes    <- c("Catch"=NA,        "TAC"=NA,        "Recommended TAC"=1)
          types     <- c("Catch"=1,         "TAC"=1,         "Recommended TAC"=0)

          g1  <- ggplot2::ggplot(data_melt, ggplot2::aes(x=Year, y=Catch/1000, linetype=variable, color=variable, shape=variable)) +
            ggplot2::geom_line(size=2) +
            ggplot2::geom_point(size=6) + 
            ggplot2::scale_linetype_manual(values=types) + 
            ggplot2::scale_shape_manual(values=shapes) + 
            ggplot2::scale_color_manual(values=colors) + 
            ggplot2::scale_y_continuous(limits=c(0, 1.2*max(data_melt$Catch/1000)), expand=c(0, 0)) +
            ggplot2::xlab('') + 
            ggplot2::ylab('Catch (x 1000)\n ') + 
            ggplot2::theme_bw() + 
            ggplot2::theme(legend.title=ggplot2::element_blank())

          data_melt <- reshape2::melt(Data$CE[, c("Year", "CPUE")], id.vars='Year', value.name='CPUE')
          colors    <- c("CPUE"="#00345D")

          g2  <-  ggplot2::ggplot(data_melt, ggplot2::aes(x=Year, y=CPUE, color=variable)) +
            ggplot2::geom_line(size=2) +
            ggplot2::scale_color_manual(values=colors) + 
            ggplot2::scale_y_continuous(limits=c(0, NA)) + 
            ggplot2::theme_bw() +
            ggplot2::xlab(' \n Year') + 
            ggplot2::theme(legend.title=ggplot2::element_blank())          

          g1 + g2 + patchwork::plot_layout(ncol=1)
        },
        width=graphWidthCatch,
        height=graphHeight)

      output$iobsPlot <- shiny::renderPlot({
          data_melt <- reshape2::melt(Data$CE[, c("Year", "CPUE")], id.vars='Year', value.name='CPUE')
          colors    <- c("CPUE"="#00345D")

          ggplot2::ggplot(data_melt, ggplot2::aes(x=data_melt$Year, y=data_melt$CPUE, color=data_melt$variable)) +
            ggplot2::geom_line(size=2) +
            ggplot2::scale_color_manual(values=colors) + 
            ggplot2::scale_y_continuous(limits=c(0, NA)) +
            ggplot2::xlab(' \n Year') + 
            ggplot2::ylab('CPUE \n ') + 
            ggplot2::theme_bw() + 
            ggplot2::theme(legend.title=ggplot2::element_blank())
        },
        width=graphWidthCPUE,
        height=graphHeight)

      output$cpuePlot <- shiny::renderPlot({
          Data$plots[["cpue_plot"]]
        },
        width=graphWidthDiag,
        height=graphHeightDiag)

      output$TAC <- shiny::renderUI({
        shiny::HTML(sprintf("<H3 class='well' style='background-color:#001a76;border-color:#000f43;color:#fff;'>Recommended TAC: %3g</H3>", Data$TAC))
      })
    }

    return ("")
  })
}

runShinyMP <- function()
{
  # Run the app ----
  shiny::shinyApp(ui = ui, server = server)
}

 