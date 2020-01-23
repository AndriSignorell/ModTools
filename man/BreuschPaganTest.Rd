\name{BreuschPaganTest}
\alias{BreuschPaganTest}
\encoding{latin1}
\title{Breusch-Pagan Test}
\description{Performs the Breusch-Pagan test against heteroskedasticity.}
\usage{
BreuschPaganTest(formula, varformula = NULL, studentize = TRUE, data = list())
}
\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
 \item{varformula}{a formula describing only the potential explanatory variables
   for the variance (no dependent variable needed). By default the same
   explanatory variables are taken as in the main regression model.}
 \item{studentize}{logical. If set to \code{TRUE} Koenker's studentized
   version of the test statistic will be used.}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which \code{BreuschPaganTest} is
called from.}
}
\details{The Breusch-Pagan test fits a linear regression model to the residuals
of a linear regression model
(by default the same explanatory variables are taken as in the main regression
model) and rejects if too much of the variance
is explained by the additional explanatory variables.

Under \eqn{H_0} the test statistic of the Breusch-Pagan test follows a
chi-squared distribution with \code{parameter} (the number of regressors without
the constant in the model) degrees of freedom.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link[lmtest]{bondyield}}, \code{\link[lmtest]{currencysubstitution}},
\code{\link[lmtest]{growthofmoney}}, \code{\link[lmtest]{moneydemand}},
 \code{\link[lmtest]{unemployment}}, \code{\link[lmtest]{wages}}.
}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{
T.S. Breusch & A.R. Pagan (1979),
A Simple Test for Heteroscedasticity and Random Coefficient Variation.
\emph{Econometrica} \bold{47}, 1287--1294

R. Koenker (1981), A Note on Studentizing a Test for Heteroscedasticity.
\emph{Journal of Econometrics} \bold{17}, 107--112.

W. Kraemer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica
}

\seealso{\code{\link{lm}}, \code{\link[car]{ncvTest}}}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\examples{
## generate a regressor
x <- rep(c(-1,1), 50)

## generate heteroskedastic and homoskedastic disturbances
err1 <- rnorm(100, sd=rep(c(1,2), 50))
err2 <- rnorm(100)

## generate a linear relationship
y1 <- 1 + x + err1
y2 <- 1 + x + err2

## perform Breusch-Pagan test
BreuschPaganTest(y1 ~ x)
BreuschPaganTest(y2 ~ x)
}
\keyword{htest}
