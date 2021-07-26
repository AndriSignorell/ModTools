\name{ROC}
\alias{ROC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a ROC curve
%%  ~~function to do ... ~~
}
\description{This is a wrapper to the main function \code{\link[pROC]{pROC}} of the \pkg{pROC} package (by Xavier Robin et al.). It builds a ROC curve and returns a \code{"roc"} object, a list of class \code{"roc"}.
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

\details{
Partial ROC is calculated following Peterson et al.
(2008; \doi{10.1016/j.ecolmodel.2007.11.008}). This function is a modification
of the PartialROC funcion, available at \url{https://github.com/narayanibarve/ENMGadgets}.
}

\value{
A data.frame containing the AUC values and AUC ratios calculated for each iteration.
}
\seealso{\code{\link[pROC]{pROC}}}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}


\references{
Peterson, A.T. et al. (2008) Rethinking receiver operating characteristic analysis applications in ecological niche modeling. Ecol. Modell., 213, 63-72.
}%% ~Make other sections like Warning with \section{Warning }{....} ~

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
