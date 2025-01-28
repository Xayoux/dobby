#' @title Provide a statistical summary for a univariate time series
#'
#' @description Provide statistical information about a univariate time serie.
#' Results contained start and end date, ADF test
#' (see \code{\link{adf_test_auto}}) and the skweness and kurtosis statistics.
#' 
#' @param df Dataframe containing the variable of interest and the date_variable.
#' @param var Character indicating the name of the variable to be studied.
#' @param serie_name Character indicating the how to name the serie.
#' @param date_variable Character indicating the variable containing the dates
#' of the serie. By default it is "date".
#' 
#' @return A tibble containing the name of the variable ; the beggining and
#' ending dates the serie ; ADF test results ; Skewness and Kurtosis statistics.
#'
#' @examples
#' # Create data
#' data <-
#'   dplyr::tibble (
#'     date = seq(from = as.Date("1980-01-01"), by = "day", length.out = 1000),
#'     normal_serie = stats::rnorm(1000, mean = 0, sd = 1),
#'     student_serie = stats::rt(1000, df = 5)
#'   )
#' 
#' # Descriptives statistics for one serie
#' summary_stat(data, "normal_serie", "Normal")
#' 
#' # Descriptives statistics for 2 series
#' purrr::map2(
#'   c("normal_serie", "student_serie"),
#'   c("Normal", "Student"),
#'   \(variable, name) summary_stat(data, variable, name)
#' ) |>
#'   purrr::list_rbind()
#'
#' @seealso
#' - [adf_test_auto()] For more informations on the ADF test.
#'
#' - \code{\link[moments]{skewness}}() For more information on the skewness
#' statistics
#'
#' - \code{\link[moments]{kurtosis}}() For more information on the kurtosis
#' statistics
#' 
#' @export
summary_stat <- function(df, var, serie_name, date_variable = "date"){
  # Check if df is a dataframe
  dobby::.check_data_frame(df, "df")

  # Check if var is a character
  dobby::.check_character(var, "var")

  # Check if var is unique
  dobby::.check_length(var, 'var', 1)

  # Check if serie_name is character
  dobby::.check_character(serie_name, "serie_name")

  # Check if serie_name is unique
  dobby::.check_length(serie_name, "serie_name", 1)

  # Check if var is in df
  .check_var_exist(df, deparse(substitute(df)), var)

  # Check if 'date' variable is in df
  .check_var_exist(df, deparse(substitute(df)), date_variable)
  
  # Create the time serie with date associated to each values
  serie_ts <-
    xts::xts(
      df |> dplyr::select(dplyr::all_of(var))  |> tidyr::drop_na(),
      order.by = df |> tidyr::drop_na(dplyr::all_of(var)) |> dplyr::pull(date)
    )

  # Df containing begginign and ending dates
  df_dates <-
    dplyr::tibble(
      start_date = stats::start(serie_ts),
      end_date = stats::end(serie_ts)
    )
  
  # ADF unit root test
  test_adf <-
    dobby::adf_test_auto(
      y = serie_ts,
      y_name = serie_name,
      message = FALSE
    )

  # Df containing skewness and kurtosis of the serie
  df_moments <-
    dplyr::tibble(
      skewness = moments::skewness(serie_ts),
      kurtosis = moments::kurtosis(serie_ts)
    )

  # Df containg all statistical informations
  df_res <-
    test_adf |>
    dplyr::bind_cols(df_dates) |>
    dplyr::bind_cols(df_moments) |>
    dplyr::mutate(
      skewness = dplyr::if_else(has_ur == TRUE, NA, skewness),
      kurtosis = dplyr::if_else(has_ur == TRUE, NA, kurtosis)
    )  |>
    dplyr::relocate(serie_name, start_date, end_date, model, nb_lags, phi_stat, has_ur)
    
  return(df_res)
}

utils::globalVariables(c("has_ur", "skewness", "kurtosis", "start_date",
                         "end_date", "model", "nb_lags", "phi_stat"))
