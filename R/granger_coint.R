#' @title Perform a cointegration granger test with moving sample
#'
#' @description The Granger cointegration test is a statistical test that can
#' be conducted to tests whether two time series share the same stochastic long
#' term trend not. If this is the case these two series are said to be
#' cointegrated.
#'
#' The Granger cointegration test has two steps :
#' - Estimate a linear regression between a serie \eqn{Y_t} and a serie \eqn{X_t}
#' (both series are non-stationnary)
#' so that : \eqn{Y_t = \alpha + \beta X_t + \epsilon_t}.
#'
#' - Take the serie of estimated residuals : \eqn{\hat{\epsilon}_t} and
#' perform a Unit root test on this serie (see
#' \code{\link[dobby]{adf_test_auto}}()) for more informations on the ADF test).
#'
#' If there is no Unit Root is the residuals serie, then \eqn{Y_t} and \eqn{X_t}
#' are said to be cointegrated. (Reminder that the null hypothesis for the ADF
#' test is the presence of at least one unit root in the serie).
#'
#' This function allow to perform the Granger cointegration test on the full
#' sample, but also on a moving sample (in a recursive or moving window way).
#'
#' 
#' @param df Dataframe containing the data and a "date" variable.
#' @param formula Character or formula indicating the formula of the
#' long-term regression.
#' @inheritParams adf_test_auto
#' @inheritParams granger_causality
#' 
#' @return A tibble containing results of rolling and full sample test
#' (statistical values and critical values).
#' 
#' @examples
#' set.seed(123456)
#' # Create 2 series
#' n <- 300
#' 
#' # Random walk
#' X <- cumsum(rnorm(n))
#' 
#' # Add a white noise to the serie X
#' Y <- X + rnorm(n, sd = 0.5)
#' 
#' # Another random walk
#' Z <- cumsum(rnorm(n))
#' 
#' data <-
#'   dplyr::tibble(
#'     date = 1:n,
#'     X = X,
#'     Y = Y,
#'     Z = Z
#'   )
#' 
#' # Test the cointegration between X and Y on the full sample and with 1%
#' # of significance
#' granger_coint(
#'   df = data,
#'   formula = "X ~ Y",
#'   type = "None",
#'   signif = 0.01
#' )
#' 
#' # Test the cointegration between X and Z with in a recursive way
#' granger_coint(
#'   df = data,
#'   formula = as.formula("X ~ Z"),
#'   type = "Moving",
#'   nb_min_obs = 100
#' )
#' 
#' # Test cointegration between X and Y then between X and Z in a moving window
#' # of 100 obs
#' purrr::map(
#'   list("X ~ Y", "X ~ Z"),
#'   \(LT_formula) granger_coint(
#'     df = data,
#'     formula = LT_formula,
#'     type = "Moving",
#'     window = 100,
#'     nb_min_obs = 100
#'   )
#' ) |>
#'   purrr::list_rbind()
#' @export
granger_coint <- function(df, formula, signif = 0.05, type = c("Moving", "None"),
                          window = integer(0), nb_min_obs = 520, set_seed = TRUE,
                          seed = 1234567) {

  # Check if df is a dataframe
  dobby::.check_data_frame(df, "df")

  # Check if formula is a character or a formula
  if (!is.character(formula) & !inherits(formula, "formula")) {
    class_formula <- class(formula)
    stop(glue::glue("\'formula\' must be a character or a formula, not a {class_formula}"))
  }

  # If formula is a character, then transform it into a formula
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }

  # Check if signif is numeric
  dobby::.check_numeric(signif, "signif")

  # Check if signif is length 1
  dobby::.check_length(signif, "signif", 1)

  # Check if signif is 0.01, 0.05 or 0.1
  if(!signif %in% c(0.01, 0.05, 0.1)){
    stop(glue::glue("signif must be one of : 0.01 ; 0.05 ; 0.1 ; not {signif}"))
  }  

  # Check if type is Moving or None
  type <- match.arg(type)

  # Check if window is > 0
  if (length(window) > 0) {
    if (window < 0){
      stop(glue::glue("window must be strictly positive. Not {window}."))
    }
  }

  # Check if nb_min_obs is numeric
  dobby::.check_numeric(nb_min_obs, "nb_min_obs")

  # Check if nb_min_obs is unique
  dobby::.check_length(nb_min_obs, "nb_min_obs", 1)

  # Check if nb_min_obs is > 0
  if (nb_min_obs < 0){
    stop(glue::glue("nb_min_obs must be strictly positive. Not {nb_min_obs}."))
  } 

  # Check if window >= nb_min_obs
  if (length(window) > 0){
    if (window < nb_min_obs){
      stop("window must be >= nb_min_obs")
    }
  }

  # Check if set_seed is logical
  dobby::.check_logical(set_seed, "set_seed")

  # Check if set_seed is unique
  dobby::.check_length(set_seed, "set_seed", 1)

  # Check if seed is numeric
  dobby::.check_numeric(seed, "seed")

  # Check if seed is unique
  dobby::.check_length(seed, "seed", 1)

  # Seet seed if needed
  if (set_seed == TRUE) {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }

  # Keep only date and needed variables
  df <-
    df |>
    dplyr::select(date, dplyr::all_of(all.vars(formula))) |>
    tidyr::drop_na()

  # Perform the long-term regression
  reg_LT <- stats::lm(formula, data = df)

  # Store variable used in the formula (allow for a for loop)
  informations <- paste(all.vars(formula), collapse = " / ")

  # Df containing results on the full sample
  df_res_all <-
    # Perform the ADF test on residuals serie
    dobby::adf_test_auto(
      reg_LT$residuals,
      glue::glue("Residuals {informations}"),
      message = FALSE,
      signif = signif
    )

  # If a moving sample is wanted
  if (type != "None") {
    # Df containing results on each sample
    df_res_moving <-
      runner::runner(
        x = df,
        k = window,
        f = \(df_moving) {
          # Store the number of observations available in the sample
          nb_obs <-
            df_moving |>
            nrow()
            
            # Do the test only if the sample size if large enough
            if (nb_obs >= nb_min_obs) {  
              # Keep the latest date on the sample (this is the date of the test)
              date <-
                df_moving |>
                dplyr::filter(dplyr::row_number() == max(dplyr::row_number())) |>
                dplyr::pull(date)

              # Perform the long-term regression
              reg_LT_moving <- stats::lm(formula, data = df_moving)

              # Store the variables used in formula (allow for a for loop)
              informations <- paste(all.vars(formula), collapse = " / ")

              # Df containing results for a particular sample
              df_res_moving <-
                dobby::adf_test_auto(
                  reg_LT_moving$residuals,
                  glue::glue("Residuals {informations}"),
                  message = FALSE,
                  signif = signif
                ) |>
                dplyr::mutate(date = date) |>
                dplyr::relocate(date)
            }
        }
      ) |>
  # Bind all individuals dataframe together
  purrr::list_rbind()

    # Keep only some variables on the full sample dataframe
    df_res_all <-
      df_res_all |>
      dplyr::select(serie_name, phi_stat, phi_crit) |>
      dplyr::rename(
        phi_stat_full = phi_stat,
        phi_crit_full = phi_crit
      )

    # Full join moving and full dataframe to have the final dataframe
    df_res <-
      df_res_moving |>
      dplyr::full_join(
        df_res_all,
        dplyr::join_by(serie_name == serie_name)
      )
  } else {
    # If there is no moving sample, final dataframe is the full sample df
    df_res <- df_res_all
  }

  # Return the final dataframe
  return(df_res)
}

utils::globalVariables(c("serie_name", "phi_crit"))
