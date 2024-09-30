\name{RobSummary}
\alias{RobSummary}
\title{Robust Summary for Linear Models
}
\description{
For poisson models with mild violation of the distribution assumption that the variance equals the mean, Cameron and Trivedi (2009) recommended using robust standard errors for the parameter estimates. The function uses the function \code{vcovHC} from the package \pkg{sandwich} to obtain the robust standard errors and calculate the p-values accordingly.
It returns a matrix containing the usual results in the model summary, comprising the parameter estimates,  their robust standard errors, p-values,  extended with the 95\% confidence interval.
}
\usage{
RobSummary(mod, conf.level = 0.95, type = "HC0")
}
\arguments{
  \item{mod}{the model for which robust standard errors should be calculated
}
  \item{conf.level}{the confidence level, default is 95\%.
}
  \item{type}{a character string specifying the estimation type. Details in \code{\link[sandwich]{vcovHC}()}. }

}
\details{Further details in \url{https://stats.oarc.ucla.edu/r/dae/poisson-regression/}
}
\value{

a \emph{p x 6} matrix with columns for the estimated coefficient, its standard error, t- or z-statistic, the corresponding (two-sided) p-value, the lower and upper confidence interval.
}
\references{Cameron, A. C. and Trivedi, P. K. (2009) Microeconometrics Using Stata. College Station, TX: Stata Press.

}
\author{Andri Signorell <andri@signorell.net>
}


\seealso{\code{\link{summary.lm}}, \code{\link{summary.glm}}
}
\examples{
r.lm <- lm(Fertility ~ ., swiss)
RobSummary(r.lm)
}
