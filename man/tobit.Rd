\name{Tobit}
\alias{Tobit}

\title{Tobit Regression}

\description{
Fitting and testing Tobit regression models for censored data.
}

\usage{
Tobit(formula, left = 0, right = Inf, dist = "gaussian",
      subset = NULL, data = list(), \dots)
}

\arguments{
  \item{formula}{a symbolic description of a regression model of type
    \code{y ~ x1 + x2 + \dots}.}
  \item{left}{left limit for the censored dependent variable \code{y}.
    If set to \code{-Inf}, \code{y} is assumed not to be left-censored.}
  \item{right}{right limit for the censored dependent variable \code{y}.
    If set to \code{Inf}, the default, \code{y} is assumed not to be right-censored.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.
    This is passed to \code{\link[survival]{survreg}}, see the respective man page for
    more details.}
  \item{subset}{a specification of the rows to be used.}
  \item{data}{a data frame containing the variables in the model.}
  \item{\dots}{further arguments passed to \code{\link[survival]{survreg}}.}
}

\details{
  The function \code{Tobit} is an alias for the \pkg{AER} function \code{\link[AER]{tobit}} (Achim Zeileis <Achim.Zeileis@R-project.org>).
  All details can be found there.
}
\author{Andri Signorell}

\value{
  An object of class \code{"Tobit"} inheriting from class \code{"survreg"}.
}

\examples{# still to do
}

\keyword{regression}
