rolling_granger <- function(df, endog_var, exog_var, selectlag, ic, lags, boot, boot.runs, lag_max,
                            set_seed, seed){
  # Seet seed if needed
  if (set_seed == TRUE){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  # Minimal number of observations for computing the VAR
  nb_obs_min <-
    df |>
    dplyr::pull(obs) |>
    unique()

  # Compute only if there is a sufficient number of obs
  if (nrow(df) >= nb_obs_min){
    date <-
      df |>
      dplyr::filter(dplyr::row_number() == max(dplyr::row_number())) |>
      dplyr::select(date)

    df_endog <-
      df |>
      dplyr::select(dplyr::all_of(endog_var))

    if (!is.null(exog_var)){
      df_exog <-
        df |>
        dplyr::select(dplyr::all_of(exog_var))
    } else {
      df_exog <- NULL
    }
    

    # Restimate the number of lags for each period if wanted
    if (selectlag == "Moving"){
      n_lag <- vars::VARselect(df_endog, lag.max = lag_max, exogen = df_exog)$selection[ic]
    } else {
      n_lag <- lags # Fixed lags
    }
    
    VAR <- vars::VAR(y = df_endog, p = n_lag, exogen = df_exog)

    # Granger test for each cause
    p_value_granger_reccursive <-
      purrr::map_dbl(
        endog_var,
        \(cause) vars::causality(VAR, cause = cause,
                                 boot = boot, boot.runs = boot.runs)$Granger$p.value
      ) |>
      dplyr::as_tibble()  |>
      dplyr::mutate(
        cause = endog_var
      ) |>
      dplyr::rename(p_value = value) |>
      dplyr::bind_cols(date) |>
      dplyr::relocate(date, cause)
  }
}


#' Perform a causality granger causality test 
#'
#' The test is conducted on the full sample, but can also be conducted on
#' a rolling window or reccursive way. Reminder that the null hypothesis of
#' a Granger test is  : the variable X do not Granger cause the variable Y.
#' 
#' @param df Dataframe containing all variables and a "date" variable.
#' @param endog_var Character vector containing the names of endogenous
#' variables.
#' @param exog_var Character vector containing the names of exogeneous
#' variables.
#' @param selectlag Character indicating if the number of lags for the VAR
#' should be reestimated each period ("Moving") or not ("Fixed"). By default
#' "Moving". Usefull only if `type == "Moving"`.
#' @param lag_max Numeric indicating the maximum number of lags to be tested
#' for the VAR.
#' @param ic Name of the information criterion used to determine the optimal
#' number of lags. Either "AIC" (the default), "HQ", "SC" or "FPE". See
#' \code{\link[vars]{VARselect}}() for more informations. 
#' @param boot Logical indicating if critical values should be computed by
#' bootstrap or not. TRUE by default. See \code{\link[vars]{causality}}() for
#' more informations.
#' @param boot.runs Numeric indicating the number of runs for the bootstrap
#' if needed. 
#' @param type Character indicating if the granger tests should be condcuted
#' in a "Moving way" (either reccursive if `window = 0` or rolling if
#' `window > 0`) or juste conducted on the full sample. "Moving" by default.
#' @param window Numeric indicating the size of the window. 0 for reccursive test.
#' It must be equal or greater than nb_min_obs.
#' @param nb_min_obs Numeric indicating the minimal number of observation the
#' window must have before the VAR is computed.
#' @param set_seed Logical indicating whether a seed should be set to always
#' obtain same random numbers. TRUE by default. The seed generation is from the
#' "L'Ecuyer-CMRG" algorithm. See \code{\link[base]{set.seed}}() and
#' \code{\link[furrr]{furrr_options}}().
#' @param seed Numeric indicating the seed number to be used.
#' 
#' @return A tibble with 5 columns :
#' - date : Latest date of the window (only if "Moving")
#' - cause : the name of the variable causing the others
#' - p_value : the p_value of the test for the causing variable at this date.
#' (Only if "Moving")
#' - p_value_full : the p_value of the test for the causing variable on the full
#' sample.
#' - information : the name of each endogeneous variable present in the VAR
#' separated by " / ". Usefull when looping.
#'
#' @examples
# Create 3 series
#' n <- 100
#' 
#' X <- numeric(n)
#' X[1] <- rnorm(1)  # Valeur initiale
#' for (t in 2:n) {
#'   X[t] <- 0.7 * X[t-1] + rnorm(1)
#' }
#' 
#' # Generation of Y_t influenced par X_t-1
#' Y <- numeric(n)
#' Y[1] <- rnorm(1)
#' for (t in 2:n) {
#'   Y[t] <- 0.5 * Y[t-1] + 0.3 * X[t-1] + rnorm(1)
#' }
#' 
#' # Generation of Z_t independant of other variables
#' Z <- numeric(n)
#' Z[1] <- rnorm(1)  # Valeur initiale
#' for (t in 2:n) {
#'   Z[t] <- 0.7 * Z[t-1] + rnorm(1)
#' }
#' 
#' data <-
#'   dplyr::tibble(
#'     date = 1:n,
#'     X = X,
#'     Y = Y,
#'     Z = Z
#'   )
#' 
#' # Granger causality on full sample with X and Y. No bootstrap.
#' granger_causality(
#'   df = data,
#'   endog_var = c("X", "Y"),
#'   ic = "AIC",
#'   boot = FALSE,
#'   type = "None"
#' )
#' 
#' # Granger causality reccursive with X, Y, Z. Bootstrap.
#' granger_causality(
#'   df = data,
#'   endog_var = c("X", "Z", "Y"),
#'   ic = "AIC",
#'   boot = TRUE,
#'   type = "Moving",
#'   nb_min_obs = 50,
#'   lag_max = 10
#' )
#' 
#' # Granger causality with XY and XZ. moving window of 50
#' purrr::map(
#'   list(c("X", "Y"), c("X", "Z")),
#'   \(var_names) granger_causality(
#'     df = data,
#'     endog_var = var_names,
#'     ic = "HQ",
#'     boot = FALSE,
#'     type = "Moving",
#'     window = 50,
#'     nb_min_obs = 50,
#'     lag_max = 10
#'   )
#' )
#' 
#' # Granger causality between X and Y with Z as an exogen variable
#' granger_causality(
#'   df = data,
#'   endog_var = c("X", "Y"),
#'   exog_var = c("Z"),
#'   ic = "AIC",
#'   boot = FALSE,
#'   type = "Moving",
#'   window = 50,
#'   nb_min_obs = 50,
#'   lag_max = 10
#' )
#'
#' @seealso
#' - \code{\link[vars]{VARselect}}() For more informations on lags selection.
#' 
#' - \code{\link[vars]{causality}}() For more informations on Granger test.
#' 
#' - \code{\link[base]{set.seed}}() For more informations on random numbers
#' 
#' - \code{\link[furrr]{furrr_options}}()
#'
#' @export
granger_causality <- function(df, endog_var, exog_var = NULL,
                              selectlag = c("Moving", "Fixed"),
                              ic = c("AIC", "HQ", "SC", "FPE"),
                              lag_max = 20,
                              boot = TRUE, boot.runs = 100,
                              type = c("Moving", "None"),
                              window = integer(0),
                              nb_min_obs = 520, set_seed = TRUE,
                              seed = 1234567
                              ){
  # Check that arguments are matching
  ic <- match.arg(ic)
  type <- match.arg(type)
  selectlag <- match.arg(selectlag)

  # Check if df is a dataframe
  dobby::.check_data_frame(df, "df")

  # Check if lag_max is numeric
  dobby::.check_numeric(lag_max, "lag_max")

  # Check if lag_max is unique
  dobby::.check_length(lag_max, "lag_max", 1)

  # Check if boot is logical
  dobby::.check_logical(boot, "boot")

  # Check if boot is unique
  dobby::.check_length(boot, "boot", 1)

  # Check if boot.runs is numeric
  dobby::.check_numeric(boot.runs, "boot.runs")

  # Check if boot.runs is unique
  dobby::.check_length(boot.runs, "boot.runs", 1)
  
  # Check if boot.run is > 0
  if (boot.runs < 0){
    stop(glue::glue("boot.runs must be strictly positive. Not {boot.runs}."))
  }

  # Check if window is numeric
  dobby::.check_numeric(window, "window")

  # Check if window is > 0
  if (length(window) > 0){
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
  if (set_seed == TRUE){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  
  # Select only wanted variables and dates
  variables <- c(endog_var, exog_var)
  df <-
    df |>
    dplyr::select(dplyr::all_of(variables), date) |>
    tidyr::drop_na() |>
    # Integrate the minimal number observations : usefull when rolling
    # Can't pass this arg otherwise
    dplyr::mutate(obs = nb_min_obs) 

  # Df used for the computation of the VAR (Remove dates)
  df_endog <-
    df |>
    dplyr::select(dplyr::all_of(endog_var))

  if (!is.null(exog_var)){
    df_exog <-
      df |>
      dplyr::select(dplyr::all_of(exog_var))
  } else {
    df_exog <- NULL
  }
  

  # Select the number of lags for the full VAR
  ic <- glue::glue("{ic}(n)")
  lags <- vars::VARselect(df_endog, lag.max = lag_max, exogen = df_exog)$selection[ic]

  # Perform the FULL VAR
  var <- vars::VAR(df_endog, p = lags, exogen = df_exog)

  # Granger test with each variable as a cause
  # Return a dataframe with one row per variable
  p_value_granger_full <-
    purrr::map_dbl(
      endog_var,
      \(cause) vars::causality(var, cause = cause,
                               boot = boot, boot.runs = boot.runs)$Granger$p.value
    ) |>
    dplyr::as_tibble()  |>
    dplyr::mutate(
      cause = endog_var
    ) |>
    dplyr::rename(p_value_full = value)

  # If reccursive or rolling window
  if (type != "None"){
    # Do the granger test for each period and return a df
    p_val_reccursive <-
      runner::runner(
        x = df,
        k = window,
        f = \(df) rolling_granger(df, endog_var = endog_var, exog_var = exog_var,
                                  selectlag = selectlag, ic = ic, lags = lags,
                                  boot = boot, boot.runs = boot.runs, lag_max = lag_max,
                                  set_seed = set_seed, seed = seed)
      ) |>
  purrr::list_rbind()


    # Create the results dataframe  
    df_res <-
      p_val_reccursive |>
      dplyr::full_join(
        p_value_granger_full,
        dplyr::join_by(cause == cause)
      ) 
  } else {
    df_res <-
      p_value_granger_full
  }

  # Add the name of all variable (allow for loop)
  df_res <-
    df_res |>
    dplyr::mutate(
      information = paste(endog_var, collapse = " / ")
    )
  
  return(df_res)
}

utils::globalVariables(c("obs", "value", "cause"))
