\name{ROC}
\alias{ROC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a ROC curve
%%  ~~function to do ... ~~
}
\description{This is a wrapper to the main function of the pROC package.
It builds a ROC curve and returns a "roc" object, a list of class "roc".
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ROC(x, resp = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a model object, or the predicted probabilities, when resp is not \code{NULL}. }
  \item{resp}{the response}
  \item{\dots}{all arguments are passed to \code{roc()}.
%%     ~~Describe \code{\dots} here~~
}
}
\value{Some values
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{\code{\link[pROC]{pROC}}}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
r.glm <- FitMod(diabetes ~ ., data = d.pima, fitfn="logit")
ROC(r.glm)

# plot ROC curves for a list of models
r.rp <- FitMod(diabetes ~ ., data = d.pima, fitfn="rpart")

# combine models to a list
mlst <- list(r.glm, r.rp)

# do the plot
for(i in seq_along(mlst))
  if(i==1){
    plot(ROC(mlst[[i]], grid=TRUE, col=c(hred, hblue)[i]))
  } else {
    lines(ROC(mlst[[i]], col=c(hred, hblue)[i]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
