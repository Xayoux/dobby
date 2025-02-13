% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/granger_causality.R
\name{granger_causality}
\alias{granger_causality}
\title{Perform a causality granger causality test}
\usage{
granger_causality(
  df,
  endog_var,
  exog_var = NULL,
  selectlag = c("Moving", "Fixed"),
  ic = c("AIC", "HQ", "SC", "FPE"),
  lag_max = 20,
  boot = TRUE,
  boot.runs = 100,
  type = c("Moving", "None"),
  window = integer(0),
  nb_min_obs = 520,
  set_seed = TRUE,
  seed = 1234567
)
}
\arguments{
\item{df}{Dataframe containing all variables and a "date" variable.}

\item{endog_var}{Character vector containing the names of endogenous
variables.}

\item{exog_var}{Character vector containing the names of exogeneous
variables.}

\item{selectlag}{Character indicating if the number of lags for the VAR
should be reestimated each period ("Moving") or not ("Fixed"). By default
"Moving". Usefull only if \code{type == "Moving"}.}

\item{ic}{Name of the information criterion used to determine the optimal
number of lags. Either "AIC" (the default), "HQ", "SC" or "FPE". See
\code{\link[vars]{VARselect}}() for more informations.}

\item{lag_max}{Numeric indicating the maximum number of lags to be tested
for the VAR.}

\item{boot}{Logical indicating if critical values should be computed by
bootstrap or not. TRUE by default. See \code{\link[vars]{causality}}() for
more informations.}

\item{boot.runs}{Numeric indicating the number of runs for the bootstrap
if needed.}

\item{type}{Character indicating if the granger tests should be condcuted
in a "Moving way" (either reccursive if \code{window = 0} or rolling if
\code{window > 0}) or juste conducted on the full sample. "Moving" by default.}

\item{window}{Numeric indicating the size of the window. 0 for reccursive test.
It must be equal or greater than nb_min_obs.}

\item{nb_min_obs}{Numeric indicating the minimal number of observation the
window must have before the VAR is computed.}

\item{set_seed}{Logical indicating whether a seed should be set to always
obtain same random numbers. TRUE by default. The seed generation is from the
"L'Ecuyer-CMRG" algorithm. See \code{\link[base]{set.seed}}() and
\code{\link[furrr]{furrr_options}}().}

\item{seed}{Numeric indicating the seed number to be used.}
}
\value{
A tibble with 5 columns :
\itemize{
\item date : Latest date of the window (only if "Moving")
\item cause : the name of the variable causing the others
\item p_value : the p_value of the test for the causing variable at this date.
(Only if "Moving")
\item p_value_full : the p_value of the test for the causing variable on the full
sample.
\item information : the name of each endogeneous variable present in the VAR
separated by " / ". Usefull when looping.
}
}
\description{
The test is conducted on the full sample, but can also be conducted on
a rolling window or reccursive way. Reminder that the null hypothesis of
a Granger test is  : the variable X do not Granger cause the variable Y.
}
\examples{
n <- 100

X <- numeric(n)
X[1] <- rnorm(1)  # Valeur initiale
for (t in 2:n) {
  X[t] <- 0.7 * X[t-1] + rnorm(1)
}

# Generation of Y_t influenced par X_t-1
Y <- numeric(n)
Y[1] <- rnorm(1)
for (t in 2:n) {
  Y[t] <- 0.5 * Y[t-1] + 0.3 * X[t-1] + rnorm(1)
}

# Generation of Z_t independant of other variables
Z <- numeric(n)
Z[1] <- rnorm(1)  # Valeur initiale
for (t in 2:n) {
  Z[t] <- 0.7 * Z[t-1] + rnorm(1)
}

data <-
  dplyr::tibble(
    date = 1:n,
    X = X,
    Y = Y,
    Z = Z
  )

# Granger causality on full sample with X and Y. No bootstrap.
granger_causality(
  df = data,
  endog_var = c("X", "Y"),
  ic = "AIC",
  boot = FALSE,
  type = "None"
)

# Granger causality reccursive with X, Y, Z. Bootstrap.
granger_causality(
  df = data,
  endog_var = c("X", "Z", "Y"),
  ic = "AIC",
  boot = TRUE,
  type = "Moving",
  nb_min_obs = 50,
  lag_max = 10
)

# Granger causality with XY and XZ. moving window of 50
purrr::map(
  list(c("X", "Y"), c("X", "Z")),
  \(var_names) granger_causality(
    df = data,
    endog_var = var_names,
    ic = "HQ",
    boot = FALSE,
    type = "Moving",
    window = 50,
    nb_min_obs = 50,
    lag_max = 10
  )
)

# Granger causality between X and Y with Z as an exogen variable
granger_causality(
  df = data,
  endog_var = c("X", "Y"),
  exog_var = c("Z"),
  ic = "AIC",
  boot = FALSE,
  type = "Moving",
  window = 50,
  nb_min_obs = 50,
  lag_max = 10
)

}
\seealso{
\itemize{
\item \code{\link[vars]{VARselect}}() For more informations on lags selection.
\item \code{\link[vars]{causality}}() For more informations on Granger test.
\item \code{\link[base]{set.seed}}() For more informations on random numbers
\item \code{\link[furrr]{furrr_options}}()
}
}
