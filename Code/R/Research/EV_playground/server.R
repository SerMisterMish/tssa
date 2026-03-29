library(ggplot2)
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
    
    mse <- function(true, pred) mean(abs(true - pred)^2)
    squared_errors <- function(true, pred) abs(true - pred)^2
    
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
    
    results <- list()
    
    calc_res_ssa <- function(tsid) {
        N <- input$N
        eps <- input$epsilon
        L <- input$L_ssa
        if (is.null(L))
            L <- (N + 1) %/% 2
        K <- N - L + 1
        set.seed(5)
        ts <- eval(str2expression(input[[tsid]]))
        evs <- ssa(ts, L = L, neig = min(L, K))$sigma ^ 2
        min_ev <- min(evs[evs >= eps]) 
        max_ev <- max(evs)
        
        results$ssa[[tsid]] <<- data.frame(
            name = c("L", "K", "Min EV", "Max EV"),
            value = c(L, K, min_ev, max_ev)
        )
    }
    
    calc_res_tssa <- function(tsid) {
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
        set.seed(5)
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
        
        results$tssa[[tsid]] <<- data.frame(
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
    
    calc_ratio_ssa <- function() {
        # ts_table <- get_res_ssa("series")
        # noise_table <- get_res_ssa("noise")
        N <- input$N; L_ssa <- input$L_ssa
        ts_table <- results$ssa[["series"]]
        noise_table <- results$ssa[["noise"]]
        
        results$ssa_ratio <<- ts_table$value[3] / noise_table$value[4]
    }
    
    calc_ratio_tssa <- function() {
        # ts_table <- get_res_tssa("series")
        # noise_table <- get_res_tssa("noise")
        N <- input$N; L_tssa <- input$L_tssa; K_tssa <- input$K_tssa
        ts_table <- results$tssa[["series"]]
        noise_table <- results$tssa[["noise"]]
        
        results$tssa_ratio <<- paste0(
            round(ts_table$value[c(4, 6, 8)] / 
                  noise_table$value[c(5, 7, 9)], 
                  4),
            collapse = ", ")
    }
    
    get_res <- function(id1, id2=NULL) {
        switch(id1,
               ssa = calc_res_ssa(id2),
               tssa = calc_res_tssa(id2),
               ssa_ratio = calc_ratio_ssa(),
               tssa_ratio = calc_ratio_tssa())
        return(results[[id1]][[ifelse(is.null(id2), TRUE, id2)]])
    }
    
    result_ssa <- reactive({get_res("ssa", "series")}) 
    result_noise_ssa <- reactive({get_res("ssa", "noise")})
    result_tssa <- reactive({get_res("tssa", "series")})
    result_noise_tssa <- reactive({get_res("tssa", "noise")})
    ratio_ssa <- reactive({get_res("ssa_ratio")})
    ratio_tssa <- reactive({get_res("tssa_ratio")})
    
    output$ssa_result <- renderTable(result_ssa(), colnames = FALSE)
    output$tssa_result <- renderTable(result_tssa(), colnames = FALSE)
    output$ssa_noise_result <- renderTable(result_noise_ssa(), colnames = FALSE)
    output$tssa_noise_result <- renderTable(result_noise_tssa(), colnames = FALSE)
    
    output$ssa_ratio <- renderText(ratio_ssa())
    output$tssa_ratio <- renderText(ratio_tssa())
    
    observe({
        withProgress( 
            message = 'MSE calculation in progress',
            value = 0, 
            {
                mses_ssa <- numeric(input$repeats)
                mses_tssa <- numeric(input$repeats)
                
                N <- input$N
                sq_err_ssa <- matrix(numeric(input$repeats * N), nrow = input$repeats)
                sq_err_tssa <- matrix(numeric(input$repeats * N), nrow = input$repeats)
                
                L_ssa <- input$L_ssa
                if (is.null(L_ssa))
                    L_ssa <- (N + 1) %/% 2
                
                L_tssa <- input$L_tssa
                K_tssa <- input$K_tssa
                if (is.null(L_tssa)) {
                    L_tssa <- (N + 2) %/% 3
                    K_tssa <- (N + 2) %/% 3
                }
                
                ts <- eval(str2expression(input$series))

                set.seed(5)
                for (i in 1:input$repeats) {
                    noise <- eval(str2expression(input$noise))
                    ts_noised <- ts + noise

                    ssa_rec <- reconstruct(
                                  ssa(ts_noised, L = L_ssa),
                                  groups = list(1:input$rank)
                               )$F1
                    mses_ssa[i] <- mse(ts, ssa_rec)
                    sq_err_ssa[i, ] <- squared_errors(ts, ssa_rec)
                    
                    tssa_rec <- tens_ssa_reconstruct(ts_noised, L = c(L_tssa, K_tssa),
                        groups = list(1:input$rank), status = FALSE)$F1
                    sq_err_tssa[i, ] <- squared_errors(ts, tssa_rec)
                    mses_tssa[i] <- mse(ts, tssa_rec)
                    
                    incProgress(1/input$repeats)
                }
                output$ssa_mse <- renderText(mean(mses_ssa))
                output$tssa_mse <- renderText(mean(mses_tssa))
                
                output$sq_err <- renderPlot({
                    ggplot(data = rbind(
                        data.frame(
                            mse = colMeans(sq_err_ssa),
                            x = 1:N,
                            method = "SSA"
                        ),
                        data.frame(
                            mse = colMeans(sq_err_tssa),
                            x = 1:N,
                            method = "T-SSA"
                        )
                    ), aes(x = x, y = mse, color = method)) +
                        geom_line()
                }) 
            } 
        ) 
    }) |> bindEvent(input$calc_mses)
}
