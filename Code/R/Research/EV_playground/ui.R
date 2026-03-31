library(shiny)
library(bslib)

fluidPage(
  title = ("T-SSA Eigenvalues Playground"),
  navset_tab(
    nav_panel(
      "Signal:Noise EV ratio, MSE",
      sidebarLayout(
        sidebarPanel(
          layout_columns(
            card(
              card_header("Common"),
              numericInput(
                "N",
                "N (series length)",
                value = 99,
                min = 3,
              ),
              numericInput(
                "epsilon",
                "Floating point error",
                value = sqrt(.Machine$double.eps),
                min = 0,
                step = sqrt(.Machine$double.eps) / 10
              ),
              textInput(
                "series",
                "Time Series",
                "cos(2 * pi * 1:N / 5)",
                placeholder = "Enter time series formula..."
              ),
              textInput(
                "noise",
                "Noise",
                "rnorm(N, 0, 1)",
                placeholder = "Enter noise formula..."
              ),
            ),
            card(
              card_header("SSA"),
              uiOutput("L_slider_ssa"),
            ),
            card(
              card_header("T-SSA"),
              uiOutput("L_slider_tssa"),
              uiOutput("K_slider_tssa"),
            )
          ),
        ),
        mainPanel(
          layout_column_wrap(column(
            width = 12,
            card(
              card_header("SSA Results"),
              card_header("Series:"),
              tableOutput("ssa_result"),
              card_header("Noise:"),
              tableOutput("ssa_noise_result"),
              card_header("Min TS EV / Max Noise EV:"),
              textOutput("ssa_ratio")
            )
          ), column(
            width = 12,
            card(
              card_header("T-SSA Results"),
              card_header("Series:"),
              tableOutput("tssa_result"),
              card_header("Noise:"),
              tableOutput("tssa_noise_result"),
              card_header("Min TS EV / Max Noise EV:"),
              textOutput("tssa_ratio")
            )
          ))
        )
      ),
      absolutePanel(
        top = "40%",
        left = "50px",
        fluidRow(h3("MSE")),
        fluidRow(layout_column_wrap(
          width = 1 / 12,
          card(
            numericInput(
              "repeats",
              "Repeats",
              value = 100,
              min = 1,
              step = 1
            ),
            numericInput(
              "rank",
              "Rank",
              value = 2,
              min = 1,
              step = 1
            ),
            actionButton("calc_mses", "Calculate MSEs"),
          ),
          card(
            h5("SSA:"),
            textOutput("ssa_mse"),
            h5("T-SSA:"),
            textOutput("tssa_mse")
          )
        ))
      ),
      fluidRow(column(width = 8, plotOutput("sq_err", height = "600px")))
    ),
    nav_panel(
      "Noise EV Distribution",
      sidebarLayout(
        sidebarPanel(
          layout_columns(
            card(
              card_header("Common"),
              numericInput(
                "distr_N",
                "N (series length)",
                value = 99,
                min = 3,
              ),
              textInput(
                "distr_noise",
                "Noise",
                "rnorm(N, 0, 1)",
                placeholder = "Enter noise formula..."
              ),
              numericInput(
                "distr_repeats",
                "Repeats (EVs sample size)",
                value = 500,
                min = 3,
              ),
              numericInput(
                "distr_bins",
                "Number of Hist bins",
                value = 20,
                min = 2,
              ),
              numericInput(
                "distr_seed",
                "RNG Seed",
                value = 5
              ),
            ),
            card(
              card_header("SSA"),
              uiOutput("distr_L_slider_ssa"),
            ),
            card(
              card_header("T-SSA"),
              uiOutput("distr_L_slider_tssa"),
              uiOutput("distr_K_slider_tssa"),
            )
          ),
          actionButton("calc_distr", "Calculate Histograms")
        ),
        mainPanel(
          fluidRow(plotOutput("distr_hist_ssa", height = "600px")),
          fluidRow(plotOutput("distr_hist_tssa", height = "600px"))
        )
      )
    )
  )
)
