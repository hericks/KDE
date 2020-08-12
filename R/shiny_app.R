library(shiny)

withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(
      paste0("#", containerId),
      where = "beforeEnd",
      ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}
# withConsoleRedirect("console", {print(sampler(reactive_density))})

ui <- fluidPage(
  pre(id = "console"),
  #Density selection


  selectInput(
    "density",
    "Choose a density:",
    c(
      "Beta distribution" = "dbeta",
      #"Chi-squared distribution" = "dchisq"
      "Exponential distribution" = "dexp",
      #"F distribution" = "df",
      #"Gamma distribution" = "dgamma",
      #"Geometric distribution" = "dgeom"
      #"Hypergeometric distribution" = "dhyper",
      #"Log-normal distribution" = "dlnorm",
      "Normal distribution" = "dnorm",
      #"Poisson distribution" = "dpois"
      #"Student's t distribution" = "dt",
      "Uniform distribution" = "dunif"
      #"Weibull distribution" = "dweibull"
      #"Custom Density" = "custom_dens"
    )
  ),
  fluidRow(
    column(
      width = 2,
      conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_min", "Min", "0"))
    ),
    column(
      width = 2,
      offset = 0,
      conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_max", "Max", "1"))
    )
  ),
  fluidRow(
    column(
      width = 2,
      conditionalPanel(condition = "input.density == 'dnorm'", numericInput("dnorm_mean", "mean", "0"))
    ),
    column(
      width = 2,
      offset = 0,
      conditionalPanel(condition = "input.density == 'dnorm'", numericInput("dnorm_sd", "sd", "1"))
    )
  ),
  fluidRow(
    column(
      width = 2,
      conditionalPanel(condition = "input.density == 'dbeta'", numericInput("dbeta_alpha", "alpha", "0.5", min=1e-5))
    ),
    column(
      width = 2,
      offset = 0,
      conditionalPanel(condition = "input.density == 'dbeta'", numericInput("dbeta_beta", "beta", "0.5", min=1e-5))
    )
  ),
  conditionalPanel(condition = "input.density == 'dexp'", numericInput("dexp_rate", "rate", "1")),


  # Kernel selection
  selectInput(
    "kernel",
    "Choose a kernel:",
    c(
      "Biweight" = "biweight",
      "Cosine" = "cosine",
      "Epanechnikov" = "epanechnikov",
      "Gaussian" = "gaussian",
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

  numericInput("subdivisions", "Number subdivisions for integration:", value = 100L),

  # parameter tweaking
  numericInput("num_samples", "Number of samples:", value = 25L),
  sliderInput("bandwidth", "Choose a bandwidth:", 0.001, 1, 0.5, 0.01),

  # bandwidth estimation
  actionButton("suggestions", "get bandwidth suggestions"),
  tableOutput("bandwidth_table"),

  plotOutput("plot"),

  fluidRow(
    column(width = 2,
           numericInput("xlim_1", "x lower limit:", value = -5L)),
    column(
      width = 2,
      offset = 0,
      numericInput("xlim_2", "x upper limit:", value = 5L)
    )
  ),
  fluidRow(
    column(width = 2,
           numericInput("ylim_1", "y lower limit:", value = 0L)),
    column(
      width = 2,
      offset = 0,
      numericInput("ylim_2", "y upper limit:", value = 1L)
    )
  )
)

server <- function(input, output, session) {
  reactive_density <- reactiveValues()
  reactive_density$support <- function() {
    switch(
      input$density,
      "dunif" = c(input$dunif_min, input$dunif_max),
      "dnorm" = NULL,
      "dbeta" = c(0,1),
      "dexp" = NULL
    )
  }
  reactive_density$fun <- function(x) {
    switch(
      input$density,
      "dunif" = dunif(x, input$dunif_min, input$dunif_max),
      "dnorm" = dnorm(x, input$dnorm_mean, input$dnorm_sd),
      "dbeta" = dbeta(x, max(c(0.3,input$dbeta_alpha), na.rm=TRUE), max(c(0.3,input$dbeta_beta), na.rm=TRUE)),
      "dexp" = dexp(x, max(c(0.1,input$dexp_rate), na.rm=TRUE))
    )
  }
  reactive_density$object <-
    function() {
      Density(reactive_density$fun, reactive_density$support())
    }

  reactive_samples <- reactiveValues()
  reactive_samples$values <- function() {
    num_of_samples <- input$num_samples
    if (is.na(input$num_samples)) {
      num_of_samples <- 1L
    }
    sampler_fun <- function(x) {
      switch(
        input$density,
        "dunif" = runif(x, input$dunif_min, input$dunif_max),
        "dnorm" = rnorm(x, input$dnorm_mean, input$dnorm_sd),
        "dbeta" = rbeta(x, max(c(0.3,input$dbeta_alpha), na.rm=TRUE), max(c(0.3,input$dbeta_beta), na.rm=TRUE)),
        "dexp" = rexp(x, max(c(0.1,input$dexp_rate), na.rm=TRUE))
      )
    }

    if (input$density == 'custom_dens') {
      M <- input$M
    }
    else{
      M <- 1
    }

    sampler <- rejection_sampling(reactive_density$object(),
                                  reactive_density$object(),
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
                                                  input$dcauchy_location, input$dcauchy_scale,
                                                  input$dexp_rate)}

  observeEvent(input_density_parameters$params(), {
    observeEvent(input$density, {
      observeEvent(input$num_samples, {
        samples <- isolate(reactive_samples$values())
        observeEvent(input$show_kernel, {
          observeEvent(input$show_kernel, {
            output$plot <- renderPlot({
              plot(
                x_grid(),
                reactive_density$fun(x_grid()),
                xlim = c(input$xlim_1, input$xlim_2),
                ylim = c(input$ylim_1, input$ylim_2),
                main = "graphic representation",
                xlab = "",
                ylab = "",
                col = "dark red",
                type = "l",
                lwd = 2

              )
              legend("topleft", legend = c("density", "KDE", "samples"), col = c("dark red","black", "royal blue"), lty = c(1,1,1), lwd = c(2,1,1), cex = 1.2)

              if (input$show_kernel) {
                lines(x_grid(), reactive_kernel$fun(x_grid()))
              }
              if (input$suggestions) {
                # bandwidth estimations
                cv_suggestion <- cross_validation(reactive_kernel$object(), samples, subdivisions=as.integer(input$subdivisions))
                pco_suggestion <- pco_method(reactive_kernel$object(), samples, subdivisions=as.integer(input$subdivisions))
                gl_suggestion <- goldenschluger_lepski(reactive_kernel$object(), samples, c(1,0.5,0.25), subdivisions=as.integer(input$subdivisions))

                lines(x_grid(),
                      reactive_kde$fun(x_grid(), samples, input$bandwidth), col = "black")
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

                bandwidth_tb <-
                  data.frame(
                    "pco_method" = pco_suggestion,
                    "crossvalidation_method" = cv_suggestion,
                    "goldenshluger_lepski_method" = gl_suggestion
                  )
                output$bandwidth_table <- renderTable(bandwidth_tb)
                legend("topright", legend = c("PCO method", "Crossvalidation", "Goldenshluger-Lepski"), col = c("dark green","violet", "steelblue2"), lty = c(1,1,1), lwd = c(2,1,1), cex = 1.2)
              }
              else {
                lines(x_grid(),
                      reactive_kde$fun(x_grid(), samples, input$bandwidth))
                points(samples,
                       integer(length(samples)),
                       pch = ".",
                       col = "blue")

              }
            })
          })
        })
      })
    })
  })
}

shinyApp(ui, server)
