library(shiny)
library(Rssa)

if (file.exists("TSSA.R")) {
    source("TSSA.R")
} else if (file.exists("../../Source/TSSA.R")) {
    source("../../Source/TSSA.R")
}

function(input, output, session) {
    output$L_slider_ssa <- renderUI({
        sliderInput("L_ssa", "L:", 2, input$N - 1, (input$N + 1) %/% 2, step = 1)
    })
    
    create_sliders <- reactive({
        if (is.null(input$K_tssa)) {
            output$K_slider_tssa <- renderUI({
                sliderInput("K_tssa",
                            "K:",
                            (input$N + 2) %/% 3,
                            input$N %/% 2,
                            (input$N + 2) %/% 3,
                            step = 1)
            })
            output$L_slider_tssa <- renderUI({
                sliderInput("L_tssa",
                            "L:",
                            2,
                            (input$N + 2) %/% 3,
                            (input$N + 2) %/% 3,
                            step = 1)
            })
        }
        else {
            output$K_slider_tssa <- renderUI({
                sliderInput("K_tssa",
                            "K:",
                            input$L_tssa,
                            min(input$N - input$L_tssa - input$K_tssa + 2, input$N %/% 2),
                            input$K_tssa,
                            step = 1)
            })
            output$L_slider_tssa <- renderUI({
                sliderInput("L_tssa",
                            "L:",
                            2,
                            min(input$K_tssa, (input$N + 2) %/% 3),
                            input$L_tssa,
                            step = 1)
            })
        }
    })
    
    get_res_ssa <- function(tsid) {
        N <- input$N
        eps <- input$epsilon
        L <- input$L_ssa
        if (is.null(L))
            L <- (N + 1) %/% 2
        K <- N - L + 1
        ts <- eval(str2expression(input[[tsid]]))
        evs <- ssa(ts, L = L, neig = min(L, K))$sigma ^ 2
        min_ev <- min(evs[evs >= eps]) 
        max_ev <- max(evs)
        
        data.frame(
            name = c("L", "K", "Min EV", "Max EV"),
            value = c(L, K, min_ev, max_ev)
        )
    }
    
    get_res_tssa <- function(tsid) {
        create_sliders()
        N <- input$N
        eps <- input$epsilon
        L <- input$L_tssa
        K <- input$K_tssa
        if (is.null(L)) {
            L <- (N + 2) %/% 3
            K <- (N + 2) %/% 3
        }
        M <- N - L - K + 2
        ts <- eval(str2expression(input[[tsid]]))
        h <- tens_ssa_decompose(ts, L = c(L, K), status = FALSE, neig = c(L, K, M))
        evs_list <- lapply(get_tensor_eigenvalues(h), \(x) x^2)
        min_ev <- lapply(evs_list, function(evs) min(evs[evs >= eps]))
        max_ev <- lapply(evs_list, function(evs) max(evs))
        evs_bound_L <- ssa(ts, L = L)$sigma ^ 2 * K
        evs_bound_K <- ssa(ts, L = K)$sigma ^ 2 * L
        evs_bound_M <- ssa(ts, L = M)$sigma ^ 2 * L
        bounds1 <- round(c(min(evs_bound_L[evs_bound_L >= eps]), max(evs_bound_L)), 2)
        bounds2 <- round(c(min(evs_bound_K[evs_bound_K >= eps]), max(evs_bound_K)), 2)
        bounds3 <- round(c(min(evs_bound_M[evs_bound_M >= eps]), max(evs_bound_M)), 2)
        
        data.frame(
            name = c(
                "L",
                "K",
                "M",
                "Min EV (dim 1)",
                "Max EV (dim 1)",
                "Min EV (dim 2)",
                "Max EV (dim 2)",
                "Min EV (dim 3)",
                "Max EV (dim 3)"
            ),
            value = c(
                L,
                K,
                M,
                min_ev[[1]],
                max_ev[[1]],
                min_ev[[2]],
                max_ev[[2]],
                min_ev[[3]],
                max_ev[[3]]
            ),
            bound = c(
                "",
                "",
                "Upper bounds:",
                bounds1[1],
                bounds1[2],
                bounds2[1],
                bounds2[2],
                bounds3[1],
                bounds3[2]
            )
        )
    }
    
    get_ratio_ssa <- function() {
        ts_table <- get_res_ssa("series")
        noise_table <- get_res_ssa("noise")
        
        ts_table$value[3] / noise_table$value[4]
    }
    
    get_ratio_tssa <- function() {
        ts_table <- get_res_tssa("series")
        noise_table <- get_res_tssa("noise")
        
        paste0(round(ts_table$value[c(4, 6, 8)] / noise_table$value[c(5, 7, 9)], 4), sep = ", ")
    }

    result_ssa <- reactive({get_res_ssa("series")})
    result_noise_ssa <- reactive({get_res_ssa("noise")})
    result_tssa <- reactive({get_res_tssa("series")})
    result_noise_tssa <- reactive({get_res_tssa("noise")})  
    ratio_ssa <- reactive({get_ratio_ssa()})
    ratio_tssa <- reactive({get_ratio_tssa()})
    
    output$ssa_result <- renderTable(result_ssa(), colnames = FALSE)
    output$tssa_result <- renderTable(result_tssa(), colnames = FALSE)
    output$ssa_noise_result <- renderTable(result_noise_ssa(), colnames = FALSE)
    output$tssa_noise_result <- renderTable(result_noise_tssa(), colnames = FALSE)
    
    output$ssa_ratio <- renderText(ratio_ssa())
    output$tssa_ratio <- renderText(ratio_tssa())
}
