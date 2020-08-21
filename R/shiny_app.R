#' Visualisation Tool for Kernel Density Estimation
#'
#' Open the visualisation app for the kernel density estimation package.
#'
#' @importFrom shiny fluidPage fluidRow column wellPanel selectInput
#'   conditionalPanel numericInput sliderInput checkboxInput actionButton tags
#'   tableOutput plotOutput reactiveValues reactive observeEvent renderPlot
#'   renderTable shinyApp
#' @importFrom stats runif dunif dnorm dbeta dexp rbeta rexp rnorm
#' @importFrom graphics par legend lines points
#'
#' @export
shiny_kde <- function(){
  ui <- fluidPage(
    fluidRow(
      style="margin-top: 15px;",
      # Settings Column
      column(align="center",
        width = 4,
    #Density selection
    wellPanel(selectInput(
      "density",
      "Choose a density:",
      c(
        "Normal distribution" = "dnorm",
        "Beta distribution" = "dbeta",
        "Exponential distribution" = "dexp",
        "Uniform distribution" = "dunif",
        "Custom Density 1" = "custom_density_1",
        "Custom Density 2" = "custom_density_2"
      )
    ),
    fluidRow(
      column(
        width = 6,
        conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_min", "Min", "0"))
      ),
      column(
        width = 6,
        offset = 0,
        conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_max", "Max", "1", min = 1e-3))
      )
    ),
    fluidRow(
      column(
        width = 6,
        conditionalPanel(condition = "input.density == 'dnorm'", numericInput("dnorm_mean", "mean", "0"))
      ),
      column(
        width = 6,
        offset = 0,
        conditionalPanel(condition = "input.density == 'dnorm'", sliderInput("dnorm_sd", "sd", 0.5, 2, 0.5, 0.1))
      )
    ),
    fluidRow(
      column(
        width = 6,
        conditionalPanel(condition = "input.density == 'dbeta'", sliderInput("dbeta_alpha", "alpha", 0.5, 5, 0.5, 0.1))
      ),
      column(
        width = 6,
        offset = 0,
        conditionalPanel(condition = "input.density == 'dbeta'", sliderInput("dbeta_beta", "beta:", 0.5, 5, 0.5, 0.1))
      )
    ),
    conditionalPanel(condition = "input.density == 'dexp'", sliderInput("dexp_rate", "rate:", 0.01, 10, 1, 0.1)),
    ),

    # Kernel selection
    wellPanel(selectInput(
      "kernel",
      "Choose a kernel:",
      c(
        "Gaussian" = "gaussian",
        "Biweight" = "biweight",
        "Cosine" = "cosine",
        "Epanechnikov" = "epanechnikov",
        "Logistic" = "logistic",
        "Rectangular" = "rectangular",
        "Sigmoid" = "sigmoid" ,
        "Silverman" = "silverman",
        "Triangular" = "triangular",
        "Tricube" = "tricube",
        "Triweight" = "triweight"
      )
    ),
    checkboxInput("show_kernel", "Plot kernel", 0),
    ),
    wellPanel(fluidRow(
      numericInput("subdivisions", "Number subdivisions for integration:", value = 1000L)),

    # parameter tweaking
    fluidRow(
      numericInput("num_samples", "Number of samples:", value = 100L)),
    fluidRow(class="row_refresh_button",
    actionButton("refresh_samples", "refresh data"),
    tags$head(tags$style("
      .row_refresh_button{height:45px;}"
    ))),
    fluidRow(
      sliderInput("bandwidth", "Choose a bandwidth:", 0.001, 1, 0.5, 0.01)),

    # bandwidth estimation
    fluidRow(
      checkboxInput("suggestions", "get bandwidth suggestions")),

    fluidRow(tableOutput("bandwidth_table"))),
      ),
    # Plotting Column
    column(width = 8,align = "center",
           plotOutput("plot"),

           wellPanel(fluidRow(
             column(width = 3,
                    numericInput("xlim_1", "x lower limit:", value = -5L)),
             column(
               width = 3,
               numericInput("xlim_2", "x upper limit:", value = 5L)
             ),

             column(width = 3,
                    numericInput("ylim_1", "y lower limit:", value = 0L)),
             column(
               width = 3,
               numericInput("ylim_2", "y upper limit:", value = 1L)
             )
           ))))
  )

  server <- function(input, output, session) {
    reactive_density <- reactiveValues()
    reactive_density$support <- function() {
      f <- function(x, default) ifelse(is.na(x), default, x)

      switch(
        input$density,
        "dunif" = c(f(input$dunif_min,0), max(f(input$dunif_max, 0), f(input$dunif_min, 0) + 0.01)),
        "dnorm" = c(f(input$dnorm_mean, 0) - 15, f(input$dnorm_mean, 0) + 15),
        "dbeta" = c(0,1),
        "dexp" = NULL,
        "custom_density_1" = c(0,1),
        "custom_density_2" = c(0,1)
      )
    }
    reactive_density$custom_density_1 <- function(x){
      ret <- 1 + sin(2*pi*x)
      ret[x < 0 | 1 < x] <- 0
      ret
    }

    reactive_density$custom_density_2 <- function(x){
      ret <- 0.75 + x^3
      ret[x < 0 | x > 1] <- 0
      ret
    }
    reactive_density$fun <- function(x) {
      f <- function(x, default) ifelse(is.na(x), default, x)
      switch(
        input$density,
        "dunif" = dunif(x, f(input$dunif_min,0), max(f(input$dunif_max, 0), f(input$dunif_min, 0) + 0.01)),
        "dnorm" = dnorm(x, f(input$dnorm_mean, 0), max(f(input$dnorm_sd, 0.1), 0.5)),
        "dbeta" = dbeta(x, max(c(0.3, input$dbeta_alpha), na.rm=TRUE), max(c(0.3, input$dbeta_beta), na.rm=TRUE)),
        "dexp" = dexp(x, max(c(0.1,input$dexp_rate), na.rm=TRUE)),
        "custom_density_1" = reactive_density$custom_density_1(x),
        "custom_density_2" = reactive_density$custom_density_2(x)

      )
    }
    reactive_density$object <-
      function() {
        switch(
          input$density,
          "dunif" = Density(reactive_density$fun, reactive_density$support(), subdivisions=1000),
          "dnorm" = Density(reactive_density$fun, reactive_density$support(), subdivisions=1000),
          "dbeta" = Density(reactive_density$fun, reactive_density$support(), subdivisions=50000),
          "dexp" = Density(reactive_density$fun, reactive_density$support(), subdivisions=10000),
          "custom_density_1" = Density(reactive_density$fun, reactive_density$support(), subdivisions=1000),
          "custom_density_2" = Density(reactive_density$fun, reactive_density$support(), subdivisions=1000)
        )
      }

    reactive_samples <- reactiveValues()
    reactive_samples$values <- function() {
      num_of_samples <- input$num_samples
      if (is.na(input$num_samples)) {
        num_of_samples <- 1L
      }
      sampler_fun <- function(x) {
        f <- function(x, default) ifelse(is.na(x), default, x)

        switch(
          input$density,
          "dunif" = runif(x, f(input$dunif_min,0), max(f(input$dunif_max, 0),f(input$dunif_min, 0) + 0.01)),
          "dnorm" = rnorm(x, f(input$dnorm_mean, 0),  max(f(input$dnorm_sd, 0.1), 0.5)),
          "dbeta" = rbeta(x, max(c(0.3,input$dbeta_alpha), na.rm=TRUE), max(c(0.3,input$dbeta_beta), na.rm=TRUE)),
          "dexp" = rexp(x, max(c(0.1,input$dexp_rate), na.rm=TRUE)),
          "custom_density_1" = runif(x,min=0, max=1),
          "custom_density_2" = runif(x,min=0, max=1)
        )
      }

      if (input$density == 'custom_density_1' || input$density == 'custom_density_2') {
        M <- 2
        helper_density <- Density(dunif, support=c(0,1))
      }
      else{
        M <- 1
        helper_density <- reactive_density$object()
      }

      sampler <- rejection_sampling(reactive_density$object(),
                                    helper_density,
                                    sampler_fun,
                                    M)
      sampler(num_of_samples)
    }

    reactive_kernel <- reactiveValues()
    reactive_kernel$object <- function() {
      switch(
        input$kernel,
        "rectangular" = rectangular,
        "triangular" = triangular,
        "epanechnikov" = epanechnikov,
        "biweight" = biweight,
        "triweight" = triweight,
        "tricube" = tricube,
        "gaussian" = gaussian,
        "cosine" = cosine,
        "logistic" = logistic,
        "sigmoid" = sigmoid,
        "silverman" = silverman
      )
    }
    reactive_kernel$fun <-
      function(x) {
        reactive_kernel$object()$fun(x)
      }

    reactive_kde <- reactiveValues()
    reactive_kde$kde <-
      function(samples, bandwidth) {
        kernel_density_estimator(reactive_kernel$object(), samples, bandwidth, input$subdivisions)
      }
    reactive_kde$fun <-
      function(x, samples, bandwidth) {
        kernel_density_estimator(reactive_kernel$object(), samples, bandwidth, input$subdivisions)$fun(x)
      }

    x_grid <-
      reactive({
        seq(
          from = input$xlim_1,
          to = input$xlim_2,
          length.out = 1000
        )
      })
    input_density_parameters <- reactiveValues()
    input_density_parameters$params <- function(){c(input$dunif_min, input$dunif_max,
                                                    input$dnorm_mean, input$dnorm_sd,
                                                    input$dbeta_alpha, input$dbeta_beta,
                                                    input$dexp_rate)}


      observeEvent(input$density, {
        observeEvent(c(input_density_parameters$params(),input$num_samples, input$refresh_samples), {
          #get new samples

          samples <- reactive_samples$values()
              output$plot <- renderPlot({
                par(mar=c(3, 2, 0.1, 0.1))
                plot(
                  x_grid(),
                  reactive_density$fun(x_grid()),
                  xlim = c(input$xlim_1, input$xlim_2),
                  ylim = c(input$ylim_1, input$ylim_2),
                  xlab = "",
                  ylab = "",
                  col = "dark red",
                  type = "l",
                  lwd = 2

                )
                legend("topright", legend = c("density", "estimator", "samples"), col = c("dark red","black", "royal blue"), lty = c(1,1,NA),pch=c(NA,NA,"."), lwd = c(2,1,1), cex = 1.2)
                if (input$show_kernel) {
                  lines(x_grid(), reactive_kernel$fun(x_grid()))
                }
                if(input$bandwidth){
                lines(x_grid(),
                      reactive_kde$fun(x_grid(), samples, input$bandwidth), col = "black")}
                if (input$suggestions) {
                  # bandwidth estimations
                  cv_suggestion <- cross_validation(reactive_kernel$object(), samples, subdivisions=as.integer(input$subdivisions))
                  pco_suggestion <- pco_method(reactive_kernel$object(), samples, subdivisions=as.integer(input$subdivisions))
                  gl_suggestion <- goldenshluger_lepski(reactive_kernel$object(), samples, subdivisions=as.integer(input$subdivisions))


                  lines(x_grid(),
                        reactive_kde$fun(x_grid(), samples, pco_suggestion), col = "dark green")
                  lines(x_grid(),
                        reactive_kde$fun(x_grid(), samples, cv_suggestion), col = "violet")
                  lines(x_grid(),
                        reactive_kde$fun(x_grid(), samples, gl_suggestion), col = "steelblue2")
                  points(samples,
                         integer(length(samples)),
                         pch = ".",
                         col = "blue")

                  legend("topleft", legend = c("PCO method", "Crossvalidation", "Goldenshluger-Lepski"), col = c("dark green","violet", "steelblue2"), lty = c(1,1,1), lwd = c(1,1,1), cex = 1.2)
                }
                else {
                  points(samples,
                         integer(length(samples)),
                         pch = ".",
                         col = "blue")
                }

                output$bandwidth_table <- renderTable({
                  if (input$suggestions)
                      data.frame(
                        "method" = c("pco_method", "cross_validation", "goldenshluger_lepski"),
                        "value" = c(pco_suggestion, cv_suggestion,gl_suggestion)
                      )

                })
              })
        })
      })
    }

  shinyApp(ui, server)
}
