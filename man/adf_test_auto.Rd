% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adf_test_auto.R
\name{adf_test_auto}
\alias{adf_test_auto}
\title{Perform the Augmented Dickey-Fuller procedure}
\usage{
adf_test_auto(
  y,
  y_name,
  lags = 20,
  selectlags = c("AIC", "BIC", "Fixed"),
  signif = 0.05,
  message = TRUE,
  return_res = TRUE
)
}
\arguments{
\item{y}{Numeric vector of the time series (not specially a ts object).}

\item{y_name}{Character indicating the name of the serie to be tested.}

\item{lags}{Numeric indicating the (maximum) number of lags.
See details for more informations.}

\item{selectlags}{Character indicating the method to select the troncation
parameter. "Fixed" for a manual choosing. "BIC" for the bayesian criterion.
"AIC" (the default) for the Akaike critetion.}

\item{signif}{Numeric indicating the level of significance. Can be, 0.01,
0.05 (the default) or 0.1.}

\item{message}{Logical indicating whether a message indicatinb the test result
should be printed out. True by default.}

\item{return_res}{Logical indicating whether the function should return
a tibble containing main informations about the test and the results.
True by default.}
}
\value{
If \code{message = TRUE}, Informations about the test rsultst will be
printed out (model selected, number of lags, H0, results). If
\code{return_res = TRUE} a tibble will be return. The tibble will contained in
one row : the name of the serie, the model selected, the number of lags used,
a logical indicating whether the serie has a unit root or not and the t-stat
of the test.
}
\description{
The Augmented Dickey-Fuller unit root test is perform
automatically.
}
\section{Augmented Dickey-Fuller test}{
The Augmented Dickey-Fuller test is a statistical test used for testing the
presence of a unit root inside a univariate time serie. It is a sequential
procedure. The procedure is as follow :
\enumerate{
\item Estimate the model (3) with drift and trend :
}

\deqn{\Delta Y_t = \alpha + \delta t + \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}

With :
\itemize{
\item \eqn{\Delta} beeing the first difference operator : \eqn{\Delta Y_t = Y_t - Y_{t-1}}
\item \eqn{Y_t} is the time serie at time t
\item \eqn{\alpha} is the drift parameter
\item \eqn{\delta} is the linear trend parameter
\item \eqn{t} is the linear trend of the serie
\item \eqn{\phi} is the parameter to be tested for the unit root presence
\item \eqn{p} is the troncation parameter (i,e., the number of lags to be included
to obtain white noise residus)
\item \eqn{\phi_j} are the parameter for the augmented part of the test
\item \eqn{\epsilon_t} is the residus at time t. The residus serie must be white
noise for the test to be valid.
}

The choice of p is crucial to have a valid test. See next section for more
information on the selection of p.

The first step is to test the significance of the linear trend. The t-stat
of \eqn{\delta} is computed and compared to the values tabulated by
Dickey-Fuller. The null hypothesis H0 : \eqn{\delta = 0}.

If the null is rejected, then we test for the significance of the
\eqn{\phi parameter}. For that we computer the t-stat and compare it to the
tabulated values. The Null hypothesis H0 is \eqn{\phi = 0} indicating that
there is at least of unit root in the serie. If the null is rejected
(i,e., \eqn{\text{t-stat}_\phi \leq \text{DF}_\phi}) then we can conclude that
there is no unit root in the serie. If the null can't be rejected
(i,e., \eqn{\text{t-stat}_\phi > \text{DF}_\phi}) then we can conclude that
the serie has at least one unit root.

If the linear trend is not significative, then :
\enumerate{
\item Estimate the model (2) with drift and without trend
}

\deqn{\Delta Y_t = \alpha + \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}

Repeat the same exercice but the parameter to be tested first is now
\eqn{\alpha}.

if \eqn{\alpha} is not significative then :
\enumerate{
\item Estimate the model without drift and trend :
}

\deqn{\Delta Y_t = \phi Y_{t-1} + \sum_{j=1}^p \phi_j \Delta Y_{t-j} + \epsilon_t}

Now, just test the significance of \eqn{\phi} and conclude.
}

\section{Select the troncation parameter}{
This function implement 2 ways of selected the troncation parameter. The first
way is to choose it yourself with \code{selectlags = "Fixed"}. In this case the
parameter \code{lags} will be the number of lags to be included in the model.

The second way is to let informations criterion (AIC or BIC) choose it. In
this case, \code{lags} will be the maximum number of lag to be included and each
possibility will be tested. The selected troncation parameter will be the
one that minimize the choosen criterion. To allow comparison between models,
the number of used observations will be the same for each model.
}

\examples{
# Create 2 series
set.seed(123)
serie1 <- rnorm(100, mean = 3, sd = 2)
serie2 <- rnorm(100, mean = 0, sd = 4) + 3

# Test for one serie : lags choosen by the user
result_test_1 <- adf_test_auto(serie1, "Serie 1", lags = 3,
                               selectlags = "Fixed", signif = 0.01,
                               message = TRUE, return_res = TRUE)

# Test with 2 series and lag choosent by BIC criterion.
result <-
  purrr::map2(
    list(serie1, serie2),
    c("Serie 1", "Serie 2"),
    \(serie, name) adf_test_auto(serie, y_name = name, lags = 15, selectlags = "BIC")
  )  |>
  purrr::list_rbind()

}
