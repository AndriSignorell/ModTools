\name{RefLevel}
\alias{RefLevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Used Reference Levels in a Linear Model
%%  ~~function to do ... ~~
}
\description{Returns all the reference levels in the factors used in a linear model. It is customer friendly to report also the reference level in lm summaries, which normally are suppressed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RefLevel(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{lm object, linear model with factors as predictors.
%%     ~~Describe \code{x} here~~
}
}
\details{For reporting tables of linear models we might want to include an information about the used reference levels, which remain uncommented in the default \code{lm} result output. \code{RefLevel()} allows to add a footnote or integrate the reference levels in the coefficient table.
%%  ~~ If necessary, more details than the description above ~~
}
\note{
It's not clear how general the used algorithm is for more exotic models. \code{\link{dummy.coef}} could in such cases be an alternative.
}
\value{a named vector containing the reference levels of all factors}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

\seealso{\code{dummy.coef}, \code{\link{Response}}, \code{\link{relevel}}, \code{\link{lm}}
}

\examples{
RefLevel(lm(breaks ~ wool + tension, data = warpbreaks))
}



