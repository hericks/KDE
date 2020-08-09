library(shiny)
library(KDE)
library(ggplot2)

withConsoleRedirect <- function(containerId, expr) {
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "output")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}
# withConsoleRedirect("console", {print(sampler(reactive_density))})


ui <- fluidPage(
  pre(id = "console"),
  #Density selection


  selectInput("density", "Choose a density:", c("Uniform distribution" = "dunif", "Normal distribution" = "dnorm", "Custom Density" = "custom_dens")),
  fluidRow(
    column(width=2,
      conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_min", "Min", "0"))),
    column(width=2, offset=0,
           conditionalPanel(condition = "input.density == 'dunif'", numericInput("dunif_max","Max", "1")))
           ),
  fluidRow(
    column(width=2,
           conditionalPanel(condition = "input.density == 'dnorm'", numericInput("dnorm_mean","mean", "0"))),
    column(width=2, offset=0,
           conditionalPanel(condition = "input.density == 'dnorm'", numericInput("dnorm_sd", "sd", "1")))),

  conditionalPanel(condition = "input.density == 'custom_dens'", textInput("Custom Density Term", "custom_density_term")),


  # Kernel selection
  selectInput("kernel", "Choose a kernel:", c("Rectangular" = "rectangular", "Triangular" = "triangular",
                                              "Epanechnikov" = "epanechnikov", "Biweight" = "biweight", "Triweight" = "triweight",
                                              "Tricube" = "tricube","Gaussian" = "gaussian", "Cosine" = "cosine",
                                              "Logistic" = "logistic", "Sigmoid" = "sigmoid" , "Silverman" = "silverman")),

  # parameter tweaking
  numericInput("num_samples", "Number of samples:", value = 25L),
  sliderInput("bandwidth", "Choose a bandwidth:", 0.001, 1, 0.5, 0.01),

  plotOutput("plot"),

  # bandwidth estimation
  actionButton("suggestions", "get bandwidth suggestions")
)

server <- function(input, output, session){
  dens <- reactive({
    dens_fun <- function(x){switch(input$density,
                             "dunif" = dunif(x, input$dunif_min, input$dunif_max),
                              "dnorm" = dnorm(x, input$dnorm_mean, input$dnorm_sd))
                            }
  })

  dens_fun <- function(x){switch(input$density,
                             "dunif" = dunif(x, input$dunif_min, input$dunif_max),
                             "dnorm" = dnorm(x, input$dnorm_mean, input$dnorm_sd))}

  reactive_density <- reactiveValues()
  reactive_density$support <- function(){switch(input$density,
                                                 "dunif" = c(input$dunif_min, input$dunif_max),
                                                 "dnorm" = NULL)}
  reactive_density$fun <- function(x){switch(input$density,
                                             "dunif" = dunif(x, input$dunif_min, input$dunif_max),
                                             "dnorm" = dnorm(x, input$dnorm_mean, input$dnorm_sd))}
  reactive_density$object <- function(){Density(reactive_density$fun, reactive_density$support())}


  reactive_samples <- reactiveValues()
  reactive_samples$values <- function(){
    sampler_fun <- function(x){switch(input$density,
                          "dunif" = runif(x, input$dunif_min, input$dunif_max),
                          "dnorm" = rnorm(x, input$dnorm_mean, input$dnorm_sd))}
    sampler <- rejection_sampling(reactive_density$object(), reactive_density$object(),sampler_fun, 1)
    sampler(input$num_samples)
  }


  reactive_kernel <- reactiveValues()
  reactive_kernel$object <- function(){switch(input$kernel,
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
                                              "silverman" = silverman)}
  reactive_kernel$fun <- function(x){reactive_kernel$object()$fun(x)}


  reactive_kde <- reactiveValues()
  reactive_kde$kde <- function(){kernel_density_estimator(reactive_kernel$object(), reactive_samples$values(), input$bandwidth)}
  reactive_kde$fun <- function(x){kernel_density_estimator(reactive_kernel$object(), reactive_samples$values(), input$bandwidth)$fun(x)}

  x_lims <- c(-10,10)
  y_lims <- c(-1,3)
  x_grid <- seq(from=x_lims[1], to=x_lims[2], length.out=1000)




  output$plot <- renderPlot({
    plot(x_grid, reactive_density$fun(x_grid))
    lines(x_grid, reactive_kernel$fun(x_grid))
    lines(x_grid, reactive_kde$fun(x_grid))
  }
  )

}

shinyApp(ui, server)




