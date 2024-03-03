\name{RobSummary}
\alias{RobSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust Summary for Linear Models
%%  ~~function to do ... ~~
}
\description{
For poisson models with mild violation of the distribution assumption that the variance equals the mean, Cameron and Trivedi (2009) recommended using robust standard errors for the parameter estimates. The function uses the function \code{vcovHC} from the package \pkg{sandwich} to obtain the robust standard errors and calculate the p-values accordingly.
It returns a matrix containing the usual results in the model summary, comprising the parameter estimates,  their robust standard errors, p-values,  extended with the 95\% confidence interval.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RobSummary(mod, conf.level = 0.95, type = "HC0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{the model for which robust standard errors should be calculated
%%     ~~Describe \code{mod} here~~
}
  \item{conf.level}{the confidence level, default is 95\%.
%%     ~~Describe \code{conf.level} here~~
}
  \item{type}{a character string specifying the estimation type. Details in \code{\link[sandwich]{vcovHC}()}. }
%%     ~~Describe \code{type} here~~

}
\details{Further details in \url{https://stats.oarc.ucla.edu/r/dae/poisson-regression/}
%%  ~~ If necessary, more details than the description above ~~
}
\value{

a \emph{p x 6} matrix with columns for the estimated coefficient, its standard error, t- or z-statistic, the corresponding (two-sided) p-value, the lower and upper confidence interval.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Cameron, A. C. and Trivedi, P. K. (2009) Microeconometrics Using Stata. College Station, TX: Stata Press.

%% ~put references to the literature/web site here ~
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{summary.lm}}, \code{\link{summary.glm}}
%% ~~objects to See Also as , ~~~
}
\examples{
r.lm <- lm(Fertility ~ ., swiss)
RobSummary(r.lm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
