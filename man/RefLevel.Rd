\name{RefLevel}
\alias{RefLevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Used Reference Levels in a Linear Model
%%  ~~function to do ... ~~
}
\description{
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
\details{For reporting tables of linear models we might want to include an information about the used reference levels, which remain uncommented in the default lm result output. \code{RefLevel()} allows to add a footnote or integrate the reference levesl in the coefficient table.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a named vector containing the reference levels of all factors
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

\seealso{\code{\link{relevel}}, \code{\link{lm}}
}
\examples{
RefLevel(lm(breaks ~ wool + tension, data = warpbreaks))
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
