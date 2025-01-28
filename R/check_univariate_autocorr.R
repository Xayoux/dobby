#' @title Check for autocorrelation in a univariate time series
#'
#' @description Create ACF and PACF graphs for a univariate serie with
#' a confidence intervall of 5%. Also
#' perform a Ljung-Box or Box-Pierce test for 1 and `n_lag`. The p.values
#' of these tests can be returned in a tibble.
#'
#' For more informations on statistical tests see :
#' \code{\link[stats]{Box.test}}(). The null hypothesis is the abscence of
#' autocorrelation for a given number of lags. 
#'
#' 
#' @param serie Numerical vector containing the time series.
#' @param serie_name Character indicating the name of the serie.
#' @param test_type Character indicating the test to perform. Either "Ljung-Box"
#' (the default) or "Box-Pierce".
#' @param n_lag Numeric (vector) indicating the number of lags to be used for the second
#' test. If it is a vector, then the test will be done for each lag.
#' @param print_tests Logical indicating if the results of the tests should be
#' printed out. True by default.
#' @param return_output Logical indicating if the graph and the result's tibble
#' should be returned. TRUE by default.
#' @param path_output NULL or character indicating the path where to save
#' the graph. If NULL (the default), graph will not be saved.
#' 
#' @return A list with 2 elements :
#' - "graph" : the ACF and PACF graph
#' - "df_res" : the tibble containing the tests p.values
#'
#' @examples
#' set.seed(123)  # Pour assurer la reproductibilité
#' 
#' # Normal independant distribution
#' normal_serie = stats::rnorm(1000, mean = 3, sd = 7)
#' 
#' # Create an AR1 process
#' n <- 1000       
#' phi <- 0.7      
#' epsilon <- stats::rnorm(n)  
#' 
#' serie_autocorr <- numeric(n)
#' serie_autocorr[1] <- epsilon[1] 
#' 
#' for (t in 2:n) {
#'   serie_autocorr[t] <- phi * serie_autocorr[t-1] + epsilon[t]
#' }
#' 
#' # Test the autocorrelation for one serie with a Ljung-Box test
#' check_univariate_autocorr(serie_autocorr, "AR(1)")
#' 
#' # Test the autocorrelation for the 2 series with the Box-Pierce test
#' # Only print test results
#' purrr::walk2(
#'   list(normal_serie, serie_autocorr),
#'   c("Normal", "AR(1)"),
#'   \(serie, name) check_univariate_autocorr(serie, name, test_type = "Box-Pierce",
#'                                            return_output = FALSE)
#' )
#' 
#' # Return only the res on tibble format and with lag 1 to 20 tested
#' df_res <-
#'   purrr::map2(
#'     list(normal_serie, serie_autocorr),
#'     c("Normal", "AR(1)"),
#'     \(serie, name) check_univariate_autocorr(serie, name, test_type = "Box-Pierce",
#'                                              return_output = FALSE,
#'                                              n_lag = 1:20)[["df_res"]]
#'   ) |>
#'   purrr::list_rbind() |>
#'   print()
#'
#' @seealso
#' \code{\link[stats]{Box.test}}() For more information on statistical tests.
#'
#' @export
check_univariate_autocorr <- function(serie, serie_name,
                                      test_type = c("Ljung-Box", "Box-Pierce"),
                                      n_lag = c(1, 20), print_tests = TRUE,
                                      return_output = TRUE, path_output = NULL){
  # Check if serie is a numeric vector
  dobby::.check_numeric(serie, "serie")

  # Check if serie_name is unique
  dobby::.check_length(serie_name, "serie_name", 1)
  
  # Check if serie_name is a character
  dobby::.check_character(serie_name, "serie_name")

  # Check test_type matching args 
  test_type <- match.arg(test_type)

  # Check if n_lag is numeric
  dobby::.check_numeric(n_lag, "n_lag")

  # Check if n_lag is positive
  if (TRUE %in% (n_lag < 1)){
    stop(glue::glue("All elements in n_lag must be > 1."))
  }

  # Check if print_tests is logical
  dobby::.check_logical(print_tests, "print_tests")

  # Check if print_tests is unique
  dobby::.check_length(print_tests, "print_tests", 1)

  # Check if return_output
  dobby::.check_logical( return_output, "return_output")

  # Check if return_output is unique
  dobby::.check_length(return_output, "return_output", 1)

  # Check if path_output is null or character
  dobby::.check_null_character(path_output, "path_output")

  # Check if path_output is unique if character
  if (is.character(path_output)){
    dobby::.check_length(path_output, "path_output")
  }
  

  # Remove NA's from the serie
  if (any(is.na(serie))){
    warnings("NA's had been removed")
    serie <- serie[-which(is.na(serie))]
  }
  
  # Compute ACF and PACF of the serie
  acf_values <- stats::acf(serie, plot = FALSE)
  pacf_values <- stats::pacf(serie, plot = FALSE)

  # Df containing ACF? PACF and confint of the serie
  acf_data <-
    # ACF and ACf confint (remove lag 0 because always = 1)
    dplyr::tibble(
      lag = acf_values$lag[-1],  
      coef = acf_values$acf[-1],  
      conf_high = 1.96 / sqrt(length(serie)),  
      conf_low = -1.96 / sqrt(length(serie))
    ) |>
    # Binary indicating if ACF coef is significant or not
    dplyr::mutate(
      signif = as.factor(dplyr::if_else(coef >= conf_high | coef <= conf_low, 1, 0)),
      type = "ACF"
    ) |>
    # PACF and PACF confint of the serie
    rbind(
      data.frame(
        lag = pacf_values$lag,   
        coef = pacf_values$acf,  
        conf_high = 1.96 / sqrt(length(serie)),  
        conf_low = -1.96 / sqrt(length(serie))   
      ) |>
      # Binary indicating if the PACF coef is significant or not
      dplyr::mutate(
        signif = as.factor(dplyr::if_else(coef >= conf_high | coef <= conf_low, 1, 0)),
        type = "PACF"
      )
    )

  # ACF and PACF graph
  graph_autocorr <-
    ggplot2::ggplot(acf_data, ggplot2::aes(x = lag, y = coef, fill = signif)) +
    ggplot2::geom_bar(stat = "identity", width = 0.1) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = conf_high), linetype = "dashed", color = "blue") + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = conf_low), linetype = "dashed", color = "blue") +  
    ggplot2::scale_fill_manual(values = c("0" = "black", "1" = "cornflowerblue")) +
    ggplot2::labs(
      title = glue::glue("ACF and PACF of {serie_name}"),
      y = "Coefficients",
      x = "Lags"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = 
        ggplot2::element_text(
          color = "black",
          size = 18,
          hjust = 1
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
    plot.title =
      ggplot2::element_text(
        color = "black",
        size = 25,
        hjust = 0.5
      ),
    strip.text =
      ggplot2::element_text(
        color = "black",
        size = 22
      ),
    strip.background =
      ggplot2::element_rect(
        fill = "lightgray"
      )
    ) +
    ggplot2::facet_wrap(~type)
  
  # Save graph if needed
  if (!is.null(path_output)){
    ggplot2::ggsave(
      path_output,
      width = 15,
      height = 8
    )
  }

  # Compute statistical test for 1 and n_lag
  test_autocorr <-
    purrr::map(
      n_lag,
      \(lag) stats::Box.test(serie, lag = lag, type = test_type)
    )

  # Store p_values
  p_values_autocorr <-
    purrr::map_dbl(
      test_autocorr,
      \(test) test$p.value
    )

  # Print tests results if needed
  if (print_tests == TRUE){
    message_test <- glue::glue("# {test_type} test for {serie_name} #")
    ligne_diez <- strrep("#", nchar(message_test))
    message(glue::glue("\n\n\n{ligne_diez}\n{message_test}\n{ligne_diez}\n\n"))

    message(glue::glue("H0 : There is no autocorrelation in '{serie_name}' for a given lag number\n\n"))

    print(test_autocorr)
  }

  name_test <- ifelse(test_type == "Ljung-Box", "LB", "BP")

  # Create a tibble containing pvalues of the tests
  df_res_autocorr <-
    dplyr::tibble(
      serie = serie_name,
      p_value = p_values_autocorr,
      lags = n_lag
    ) |>
    tidyr::pivot_wider(
      names_from = lags,
      values_from = p_value,
      names_prefix = glue::glue("{name_test}_")
    )

  # return graph and df if needed
  if (return_output == TRUE){
    return(
      list(
        graph = graph_autocorr,
        df_res = df_res_autocorr
      )
    )
  }
}

utils::globalVariables(c("coef", "conf_high", "conf_low", "lag", "lags",
                         "p_value"))
