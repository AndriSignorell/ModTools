\name{Response}
\alias{Response}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract the Response from Several Models
%%  ~~function to do ... ~~
}
\description{Time after time, in the course of our daily work, we experience that the response variable is hidden very deeply in the object. This again leads to superfluous consultation of the documentation.
\code{Reponse()} relieves us of this work.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Response(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the model to use
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{more arguments
%%     ~~Describe \code{\dots} here~~
}
}
\details{The function implements the extraction of the response variables for all the models listed in the package's help text.
%%  ~~ If necessary, more details than the description above ~~
}
\value{the response of model x
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{model.frame}}, \code{\link{model.response}}, \code{\link{RefLevel}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.rpart <- FitMod(diabetes ~ ., d.pima, fitfn="rpart")
Response(r.rpart)

# up to the attribute "response" this is the same
identical(StripAttr(Response(r.rpart), "response"),
          model.response(model.frame(r.rpart)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
