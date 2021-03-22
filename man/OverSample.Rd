\name{Over-/Undersample}
\alias{OverSample}
\alias{UnderSample}

\title{Oversample and Undersample
%%  ~~function to do ... ~~
}
\description{For classification purposes we might want to have balanced datasets. If the response variable has not a prevalence of 50\%, we can sample records for getting as much response A cases as response B. This is called oversample. Undersample means to sample the (lower) number of cases A from the records of case B.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
OverSample(x, vname)
UnderSample(x, vname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data frame containing predictors and response
%%     ~~Describe \code{x} here~~
}
  \item{vname}{the name of the response variable to be used to over/undersample
%%     ~~Describe \code{vname} here~~
}
}
\value{a data frame with balanced response variable
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}
\seealso{ \code{\link{BestCut}}
}
\examples{
OverSample(d.pima2, "diabetes")

UnderSample(d.pima2, "diabetes")
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
