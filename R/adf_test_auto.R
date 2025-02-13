#' @title Test the Unit root for on ADF model
#'
#' @description Perform the ADF test for one of the 3 possible ADF models. 
#' 
#' @param y_name String indicating the name of the serie.
#' @param y_diff Numeric vector containing the serie in first difference.
#' @param tt Numeric vector containing the linear trend.
#' @param y_lag_1 Numeric vector containing the serie lagged of one period.
#' @param x Numeric matrix containing the lagged serie in difference for all
#' lags.
#' @param lags Numeric indicating the maximum number of lags to be tested if
#' `selectlags == "AIC"` or if `selectlags == "BIC"`. if `selectlags == "Fixed"`,
#' the it is the number of lag to be include. 
#' @param selectlags String indicating if the number of lags is choose by the
#' user (`selectlags = "Fixed"`) or is choosent by minimizing informations
#' criterion (`selectlags == "AIC"` or `selectlags == "BIC"`).
#' @param num_model String indicating the model to be tested : "1" for the model
#' without trend and drift. "2" for the model with drift and without trend.
#' "3" for the model with drift and trend.
#' @param signif Numeric indicating the significance threshold to be used. Must
#' be one of : 0.01, 0.05 or 0.1.
#' @param n Numeric indicating the number of observations
#' 
#' @return A list containing the message of the result ; a logical indicating
#' if the tested parameter (trend or drift) is significative (always TRUE for
#' model 1) ; A string indicating the number of the model ; The t-stat of the
#' phi parameter (NA if the trend or drift parameter is not significative) ;
#' the number of lags choosen.
#'
#' @examples
#' # Not needed
#'
model_test <- function(y_name, y_diff, tt, y_lag_1, x, lags, selectlags, num_model, signif, n){
  # Define critical values for 1%, 5% and 10%
  # Critical values for n > 500
  df_crit_values_inf <-
    data.frame(
      tt = c(3.46, 2.78, 2.38),
      const = c(3.18, 2.52, 2.16),
      phi1 = c(-2.58, -1.95, -1.62), # modèle sans constante ni tendace
      phi2 = c(-3.43, -2.86, -2.57), # modèle avec constante sans tendance
      phi3 = c(-3.96, -3.41, -3.12) # Modèle avec constante et tendance
    )

  # Critical values for 250 < n <= 500
  df_crit_values_500 <-
    data.frame(
      tt = c(3.48, 2.78, 2.38),
      const = c(3.18, 2.52, 2.16),
      phi1 = c(-2.58, -1.95, -1.62), # modèle sans constante ni tendace
      phi2 = c(-3.44, -2.87, -2.57), # modèle avec constante sans tendance
      phi3 = c(-3.98, -3.42, -3.13) # Modèle avec constante et tendance
    )

  # Critical values for 100 < n <= 250
  df_crit_values_250 <-
    data.frame(
      tt = c(3.49, 2.79, 2.38),
      const = c(3.19, 2.53, 2.16),
      phi1 = c(-2.58, -1.95, -1.62), # modèle sans constante ni tendace
      phi2 = c(-3.46, -2.88, -2.57), # modèle avec constante sans tendance
      phi3 = c(-3.99, -3.43, -3.13) # Modèle avec constante et tendance
    )

  # Critical values for <= 100
  df_crit_values_100 <-
    data.frame(
      tt = c(3.53, 2.79, 2.38),
      const = c(3.22, 2.54, 2.17),
      phi1 = c(-2.60, -1.95, -1.61), # modèle sans constante ni tendace
      phi2 = c(-3.51, -2.89, -2.58), # modèle avec constante sans tendance
      phi3 = c(-4.04, -3.45, -3.15) # Modèle avec constante et tendance
    )

  # Choose the right dataframe to use based on the number of observations
  number_obs <-
    dplyr::case_when(
      n <= 100 ~ "100",
      n <= 250 ~ "250",
      n <= 500 ~ "500",
      .default = "inf"
    )

  df_crit_values <-
    switch(
      number_obs,
      "100" = df_crit_values_100,
      "250" = df_crit_values_250,
      "500" = df_crit_values_500,
      "inf" = df_crit_values_inf,
      )

  crit_value_phi <- NULL

  # Define the number of the column to be used based on the significance level
  index_signif <- ifelse(signif == 0.01, 1, ifelse(signif == 0.05, 2, 3))

  # Define the parameter to be tested based on the model (drift or trend)
  tested_param <-
    switch(
      num_model,
      "2" = "(Intercept)",
      "3" = "tt"
    )

  # Store the critical value for the drift or trend
  crit_value_tested_param <-
    switch(
      num_model,
      "2" = df_crit_values[index_signif,"const"],
      "3" = df_crit_values[index_signif,"tt"]
    )

  # Define the formula based on the model selected
  model_formula <-
    switch(
      num_model,
      "1" = stats::formula("y_diff ~ 0 + y_lag_1 + y_diff_lag"),
      "2" = stats::formula("y_diff ~ y_lag_1 + y_diff_lag"),
      "3" = stats::formula("y_diff ~ tt + y_lag_1 + y_diff_lag"),
      )

  # Initiate the final number of lags to lags (for the Fixed case)
  nb_lags_model <- lags

  # If selection not Fixed, then choose the right number of lags by testing them
  if (selectlags != "Fixed"){
    v_criterion <- rep(NA, lags)
    for (lag in 2:lags){
      y_diff_lag <- x[,2:lag] #From 2nd columns, bc the 1st is the diff not lagged
      result <- stats::lm(model_formula)
      v_criterion[lag] <- stats::AIC(result, k = switch(selectlags, AIC = 2, BIC = log(length(y_diff))))
    }
    nb_lags_model <- which.min(v_criterion) # Define the number of lags to be used
  }

  # Estimate the real model with the right number of lags
  y_diff_lag <- x[,2:nb_lags_model]
  model_adf <- stats::lm(model_formula)
  summary_model_adf <- summary(model_adf)

  # Intiate the result of the significance test on the tested param on TRUE
  # Used if model 1 because no param had to be tested
  signif_tested_param <- TRUE

  # If the model has a drift or a trend, test the significance of this param
  if (num_model != "1"){
    t_stat_tested_param <- summary_model_adf$coefficients[tested_param, "t value"]

    signif_tested_param <- (abs(t_stat_tested_param) >= crit_value_tested_param)
  }

  # Intiate results variables in case of the param is not significant
  logical_ur <- NA
  t_stat_phi <- NA
  message_result <- NA

  # If the param is significant then test the unit root
  if (signif_tested_param == TRUE){
    t_stat_phi <- summary_model_adf$coefficients["y_lag_1", "t value"]

    # Select the critical value for the unit root
    crit_value_phi <-
      switch(
        num_model,
        "1" = df_crit_values[index_signif, "phi1"],
        "2" = df_crit_values[index_signif, "phi2"],
        "3" = df_crit_values[index_signif, "phi3"]
      )

    confidence <- signif * 100
    # H0 is rejected
    if (t_stat_phi <= crit_value_phi){
      message_result <- glue::glue("H0 is rejected with for {confidence}% confidence. {y_name} serie has no unit root.")
      
      logical_ur <- FALSE
    }
    # HO can't be rejected
    else {
      message_result <- glue::glue("H0 can't be rejected for a {confidence}% confidence. {y_name} serie has at least one unit root.")

      logical_ur <- TRUE
    }
  }

  if (is.null(crit_value_phi)) {
    crit_value_phi <- NA
  }

  # Return a list will usefull informations
  return(
    list(
      message_result = message_result,
      signif_tested_param = signif_tested_param,
      logical_ur = logical_ur,
      num_model = num_model,
      t_stat_phi = t_stat_phi,
      nb_lags_model = (nb_lags_model - 1),
      crit_value_phi = crit_value_phi
    )
  )
}




#' @title Perform the Augmented Dickey-Fuller procedure
#'
#' @description The Augmented Dickey-Fuller unit root test is perform
#' automatically.
#'
#' @details
#' # Augmented Dickey-Fuller test
#' The Augmented Dickey-Fuller test is a statistical test used for testing the
#' presence of a unit root inside a univariate time serie. It is a sequential
#' procedure. The procedure is as follow :
#'
#' 1) Estimate the model (3) with drift and trend :
#'
#' \deqn{\Delta Y_t = \alpha + \delta t + \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}
#'
#' With :
#' - \eqn{\Delta} beeing the first difference operator : \eqn{\Delta Y_t = Y_t - Y_{t-1}}
#' - \eqn{Y_t} is the time serie at time t
#' - \eqn{\alpha} is the drift parameter
#' - \eqn{\delta} is the linear trend parameter
#' - \eqn{t} is the linear trend of the serie
#' - \eqn{\phi} is the parameter to be tested for the unit root presence
#' - \eqn{p} is the troncation parameter (i,e., the number of lags to be included
#' to obtain white noise residus)
#' - \eqn{\phi_j} are the parameter for the augmented part of the test
#' - \eqn{\epsilon_t} is the residus at time t. The residus serie must be white
#' noise for the test to be valid.
#'
#' The choice of p is crucial to have a valid test. See next section for more
#' information on the selection of p.
#'
#' The first step is to test the significance of the linear trend. The t-stat
#' of \eqn{\delta} is computed and compared to the values tabulated by
#' Dickey-Fuller. The null hypothesis H0 : \eqn{\delta = 0}.
#'
#' If the null is rejected, then we test for the significance of the
#' \eqn{\phi parameter}. For that we computer the t-stat and compare it to the
#' tabulated values. The Null hypothesis H0 is \eqn{\phi = 0} indicating that
#' there is at least of unit root in the serie. If the null is rejected
#' (i,e., \eqn{\text{t-stat}_\phi \leq \text{DF}_\phi}) then we can conclude that
#' there is no unit root in the serie. If the null can't be rejected
#' (i,e., \eqn{\text{t-stat}_\phi > \text{DF}_\phi}) then we can conclude that
#' the serie has at least one unit root.
#'
#' If the linear trend is not significative, then :
#'
#' 2) Estimate the model (2) with drift and without trend
#'
#' \deqn{\Delta Y_t = \alpha + \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}
#'
#' Repeat the same exercice but the parameter to be tested first is now
#' \eqn{\alpha}.
#'
#' if \eqn{\alpha} is not significative then :
#'
#' 3) Estimate the model without drift and trend :
#'
#' \deqn{\Delta Y_t = \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}
#' 
#' Now, just test the significance of \eqn{\phi} and conclude.
#'
#' # Select the troncation parameter
#'
#' This function implement 2 ways of selected the troncation parameter. The first
#' way is to choose it yourself with `selectlags = "Fixed"`. In this case the
#' parameter `lags` will be the number of lags to be included in the model.
#'
#' The second way is to let informations criterion (AIC or BIC) choose it. In
#' this case, `lags` will be the maximum number of lag to be included and each
#' possibility will be tested. The selected troncation parameter will be the
#' one that minimize the choosen criterion. To allow comparison between models,
#' the number of used observations will be the same for each model.
#'
#' 
#' @param y Numeric vector of the time series (not specially a ts object).
#' @param y_name Character indicating the name of the serie to be tested.
#' @param lags Numeric indicating the (maximum) number of lags.
#' See details for more informations.
#' @param selectlags Character indicating the method to select the troncation
#' parameter. "Fixed" for a manual choosing. "BIC" for the bayesian criterion.
#' "AIC" (the default) for the Akaike critetion.
#' @param signif Numeric indicating the level of significance. Can be, 0.01,
#' 0.05 (the default) or 0.1.
#' @param message Logical indicating whether a message indicatinb the test result
#' should be printed out. True by default.
#' @param return_res Logical indicating whether the function should return
#' a tibble containing main informations about the test and the results.
#' True by default.
#' 
#' @return If `message = TRUE`, Informations about the test rsultst will be
#' printed out (model selected, number of lags, H0, results). If
#' `return_res = TRUE` a tibble will be return. The tibble will contained in
#' one row : the name of the serie, the model selected, the number of lags used,
#' a logical indicating whether the serie has a unit root or not and the t-stat
#' of the test.
#'
#' @examples
#' # Create 2 series
#' set.seed(123)
#' serie1 <- rnorm(100, mean = 3, sd = 2)
#' serie2 <- rnorm(100, mean = 0, sd = 4) + 3
#' 
#' # Test for one serie : lags choosen by the user
#' result_test_1 <- adf_test_auto(serie1, "Serie 1", lags = 3,
#'                                selectlags = "Fixed", signif = 0.01,
#'                                message = TRUE, return_res = TRUE)
#' 
#' # Test with 2 series and lag choosent by BIC criterion.
#' result <-
#'   purrr::map2(
#'     list(serie1, serie2),
#'     c("Serie 1", "Serie 2"),
#'     \(serie, name) adf_test_auto(serie, y_name = name, lags = 15, selectlags = "BIC")
#'   )  |>
#'   purrr::list_rbind()
#'
#' @export
adf_test_auto <- function (y, y_name, lags = 20, selectlags = c("AIC", "BIC", "Fixed"), signif = 0.05,
                           message = TRUE, return_res = TRUE){
  # Verify parameters of the functions
  selectlags <- match.arg(selectlags)

  if (ncol(as.matrix(y)) > 1) 
    stop("\ny is not a vector or univariate time series.\n")

  # Remove NA if they exists and warn the user
  if (any(is.na(y))){
    y <- y[-which(is.na(y))]
  }
  warnings("NA's had been removed.")

  y <- as.vector(y)

  lag <- as.integer(lags)
  if (lag < 0) {
    stop("\nLags must be set to an non negative integer value.\n")
  }
  
  if (!is.numeric(signif)){
    class_signif <- class(signif)
    stop(glue::glue("signif must be numeric, not {class_signif}"))
  }

  if(!signif %in% c(0.01, 0.05, 0.1)){
    stop(glue::glue("signif must be one of : 0.01 ; 0.05 ; 0.1 ; not {signif}"))
  }

  if (!is.logical(message)){
    class_message <- class(message)
    stop(glue::glue("message must be numeric, not {class_message}"))
  }

  if (!is.logical(return_res)){
    class_return_res <- class(return_res)
    stop(glue::glue("return_res must be numeric, not {class_return_res}"))
  }

  if (!is.character(y_name)){
    class_y_name <- class(y_name)
    stop(glue::glue("y_name must be character, not {class_y_name}"))
  }

  # Prepare the data
  lags <- lags + 1 # Number of lags + 1 bc the first column is not lagged
  y_diff <- diff(y) # First difference of the serie
  n <- length(y_diff) # Number of observations of the first difference
  x <- stats::embed(y_diff, lags) # Matrix with lagged first difference series
  y_diff <- x[, 1] # Only the first difference not lagged
  y_lag_1 <- y[lags:n] # Lag of the serie in level (start at lags bc firsts obs are lost)
  tt <- lags:n # linear trend

  # 1) Model with trend and drift
  test_model <-
    model_test(
      y_name = y_name,
      y_diff = y_diff,
      tt = tt,
      y_lag_1 = y_lag_1,
      x = x,
      lags = lags,
      selectlags = selectlags,
      num_model = "3",
      signif = signif,
      n= n
    )

  # 2) Model with drift and no trend is trend not significant
  if (test_model$signif_tested_param == FALSE) {
    test_model <-
      model_test(
        y_name = y_name,
        y_diff = y_diff,
        tt = tt,
        y_lag_1 = y_lag_1,
        x = x,
        lags = lags,
        selectlags = selectlags,
        num_model = "2",
        signif = signif,
        n = n
      )

    # 3) Model with no drift or trend if drift and trend not significant
    if (test_model$signif_tested_param == FALSE){
      test_model <-
        model_test(
          y_name = y_name,
          y_diff = y_diff,
          tt = tt,
          y_lag_1 = y_lag_1,
          x = x,
          lags = lags,
          selectlags = selectlags,
          num_model = "1",
          signif = signif,
          n = n
        )
    }
  }
  
  # Displya informations on the test results
  if (message == TRUE){
    message_test_adf <- glue::glue("# Augmented Dickey-Fuller test for {y_name} serie #")
    ligne_diez <- strrep("#", nchar(message_test_adf))

    type_model <-
      switch(
        test_model$num_model,
        "1" = "Without drift and trend",
        "2" = "With drift and without trend",
        "3" = "With drift and trend"
      )
    
    message(glue::glue("\n\n\n{ligne_diez}\n{message_test_adf}\n{ligne_diez}\n\n"))
    message(glue::glue("model selected : {type_model}\n\n"))
    message(glue::glue("Lags number : {test_model$nb_lags_model}\n\n"))
    message(glue::glue("H0 : {y_name} serie has at least one unit root\n\n"))      
    message(glue::glue("t stat : {test_model$t_stat_phi}\n\n"))

    if (test_model$logical_ur == TRUE){
      cli::cli_alert_danger(test_model$message_result)
    } else{
      cli::cli_alert_success(test_model$message_result)
    }
    message("\n")
  }

  # Return if needed tibble with test results
  df_res_ur <-
    dplyr::tibble(
      serie_name = y_name,
      nb_lags = test_model$nb_lags_model,
      model = test_model$num_model,
      has_ur = test_model$logical_ur,
      phi_stat = test_model$t_stat_phi,
      phi_crit = test_model$crit_value_phi
    )

  if(return_res == TRUE){
    return(df_res_ur)
  }
}





