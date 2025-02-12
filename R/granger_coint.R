#' @title Perform a cointegration granger test
#'
#' @description The test is conducted on the full sample, but can also be conducted on
#' a rolling window or reccursive way. A long term regression is run with the
#' formula provided. An ADF test (see \code{\link[dobby]{adf_test_auto}}()) is
#' then conducted on the residuals. Is there is no unit root, then the series
#' are said to be cointegrated. 
#' 
#' @param df Dataframe cointaing data and a "date" variable.
#' @param formula Character or formula indicating the formula of the
#' long-term regression.
#' @inheritParams granger_causality
#' 
#' @return A tibble containing results of rolling and full sample test
#' (statistical values and critical values).
#' @examples
#' # No examples yet
#' @export
granger_coint <- function(df, formula, type = c("Moving", "None"),
                          window = integer(0), nb_min_obs = 520, set_seed = TRUE,
                          seed = 1234567) {
  type <- match.arg(type)

  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }

  # Seet seed if needed
  if (set_seed == TRUE){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }

  df <-
    df |>
    dplyr::select(date, dplyr::all_of(all.vars(formula))) |>
    tidyr::drop_na()
  
  reg_LT <- stats::lm(formula, data = df)

  informations <- paste(all.vars(formula), collapse = " / ")

  df_res_all <-
    dobby::adf_test_auto(
      reg_LT$residuals,
      glue::glue("Residuals {informations}"),
      message = FALSE
    )

  if (type != "None") {
    df_res_moving <-
      runner::runner(
        x = df,
        k = window,
        f = \(df_moving) {
          nb_obs <-
            df_moving |>
            nrow()
          
          if (nb_obs >= nb_min_obs){  
            date <-
              df_moving |>
              dplyr::filter(dplyr::row_number() == max(dplyr::row_number())) |>
              dplyr::pull(date)
            
            reg_LT_moving <- stats::lm(formula, data = df_moving)

            informations <- paste(all.vars(formula), collapse = " / ")

            df_res_moving <-
              dobby::adf_test_auto(
                reg_LT_moving$residuals,
                glue::glue("Residuals {informations}"),
                message = FALSE
              ) |>
              dplyr::mutate(date = date) |>
              dplyr::relocate(date)
          }
        }
      ) |>
  purrr::list_rbind()

    df_res_all <-
      df_res_all |>
      dplyr::select(serie_name, phi_stat, phi_crit) |>
      dplyr::rename(
        phi_stat_full = phi_stat,
        phi_crit_full = phi_crit
      )
    
    df_res <-
      df_res_moving |>
      dplyr::full_join(
        df_res_all,
        dplyr::join_by(serie_name == serie_name)
      )
  } else {
    df_res <- df_res_all
  }
  return(df_res)
}

utils::globalVariables(c("serie_name", "phi_crit"))
