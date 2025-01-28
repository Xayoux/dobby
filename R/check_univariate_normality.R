#' @title Check the normality of a univariate serie
#'
#' @description Perform statistical tests and draw graphs to test whether a serie
#' follow a normal distribution or not.
#'
#' @details 2 type of graphs are drawn :
#' - A density graph wich compare the empirical density of the serie to the
#' normal density where the normal distribution is characterized by the same mean and
#' standard deviation than the real serie.
#' - A QQplot graph wich compare The empirical quantiles distribution to
#' theoretical quantiles if the serie follow a normal distribution.
#'
#' 3 tests statistical tests are performed :
#' - Agostino's test : H0 : the serie has a skewness of 0 ; H1 : the serie
#' has a skewness different of 0.
#' - Anscombe's test : H0 : the serie has a kurtosis of 3 ; H1 : the serie has
#' a kurtosis different of 3
#' - Jarque-Berra's tests : H0 : the serie follow a normal distribution ; H1 :
#' the serie does not follow the normal distribution.
#'
#'
#' @importFrom stats dnorm
#' 
#' @param serie A numerical vector containing the serie to be tested out.
#' @param name Character indicating the name of the serie.
#' @param print_tests Logical indicating whether statistical tests should be
#' printed or not. TRUE by default.
#' @param return_output Logical indicating whether the graphs output should be
#' returned or not. TRUE by default.
#' @param path_output NULL or character indicating the path to save the graphs.
#' If NULL graphs will not be saved. NULL by default.
#'
#' @return Return (if wanted) the normality and qqplot graph in one unique ggplot
#' object.
#'
#' @examples
#' normal_serie <- stats::rnorm(1000, mean = 0, sd = 1)
#' check_univariate_normality(
#'   normal_serie,
#'   "Normal serie"
#' )
#' 
#' student_serie <- stats::rt(1000, df = 5)
#' check_univariate_normality(
#'   student_serie,
#'   "Student serie"
#' )
#'
#' @seealso
#' - \code{\link[moments]{agostino.test}}() For the Agostino's skewness test.
#' 
#' - \code{\link[moments]{anscombe.test}}() For the Anscombe's kurtosis test.
#' 
#' - \code{\link[moments]{jarque.test}}() For the Jarque-Berra's normality test.
#' 
#' @export
check_univariate_normality <- function(serie, name, print_tests = TRUE,
                                       return_output = TRUE, path_output = NULL){
  # Remove Na's
  serie <- serie[which(!is.na(serie))]

  # Compute statistics to calibrate the normal density
  mean_serie <- mean(serie, na.rm = TRUE)
  median_serie <- stats::median(serie, na.rm = TRUE)
  sd_serie <- stats::sd(serie, na.rm = TRUE)

  ## Create density graph
  graph_densite <-
    serie |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(x = serie)) +
    ## Add empirical density
    ggplot2::geom_density(ggplot2::aes(color = "Empirical density"), fill = "lightgray", na.rm = TRUE) +
    ## Vertical line at the median
    ggplot2::geom_vline(xintercept = median_serie, color = "red", linetype = "solid") +
    ## Draw the normal density
    ggplot2::stat_function(
      fun = dnorm, 
      args = list(mean = mean_serie, sd = sd_serie), 
      ggplot2::aes(color = "Normal density"), 
      linetype = "dashed"
    ) +
    ggplot2::scale_color_manual(
      values = c("Empirical density" = "black", "Normal density" = "blue")
    ) + 
    ggplot2::labs(
      color = "Density",
      title = glue::glue("Normal and {name} empirical densities"),
      y = "",
      x = glue::glue("{name}")
    ) +  
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = 
        ggplot2::element_text(
          color = "black",
          size = 18,
          hjust = 0.5
        ),
      axis.title.x = 
        ggplot2::element_text(
          color = "black",
          size = 22,
          vjust = -0.5
        ),
      axis.text.y = 
        ggplot2::element_text(
          color = "black",
          size = 18
        ),
      axis.title.y =
        ggplot2::element_text(
          color = "black",
          size = 22
        ),
      plot.title = ggplot2::element_text(color = "black", hjust = 0.5, size = 24),
      legend.text =
        ggplot2::element_text(
          color = "black",
          size = 18
        ),
      legend.title =
        ggplot2::element_text(
          color = "black",
          size = 22
        ),
      legend.position = "right",
      legend.key.size = ggplot2::unit(1, "cm")
    )


  ## QQ-plot
  graph_qqplot <-
    serie |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(sample = serie)) +
    ggplot2::stat_qq(color = "black", na.rm = TRUE) +  # QQ-plot points
    ggplot2::stat_qq_line(color = "red", na.rm = TRUE) +  # Theoretical reference line
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = glue::glue("{name} QQplot"), 
      x = "Theoretical quantiles", 
      y = "Data quantiles"
    ) +
    ggplot2::theme(
      axis.text.x = 
        ggplot2::element_text(
          color = "black",
          size = 18,
          hjust = 0.5
        ),
      axis.title.x = 
        ggplot2::element_text(
          color = "black",
          size = 22,
          vjust = -0.5
        ),
      axis.text.y = 
        ggplot2::element_text(
          color = "black",
          size = 18
        ),
      axis.title.y =
        ggplot2::element_text(
          color = "black",
          size = 22
        ),
      plot.title = ggplot2::element_text(color = "black", hjust = 0.5, size = 24)
    )

  # Combine the 2 plots together
  graph <-
    patchwork::wrap_plots(graph_densite, graph_qqplot)
   

  # Print statistical tests of normality
  if (print_tests == TRUE){
    message_test <- glue::glue("# Normality tests for {name} serie #")
    diez_line <- strrep("#", nchar(message_test))

    message(glue::glue("\n\n\n{diez_line}\n{message_test}\n{diez_line}\n\n"))
    
    print(moments::agostino.test(serie))
    print(moments::anscombe.test(serie))
    print(moments::jarque.test(serie))
  }
  

  # Save graph if needed
  if (!is.null(path_output)){
    ggplot2::ggsave(
      path_output,
      width = 15,
      height = 8
    )
  }

  # Return the graph
  if (return_output == TRUE){
    return(graph)
  }
}
