\name{ROC}
\alias{ROC}
\title{Build a ROC curve
}
\description{This is a wrapper to the main function \code{\link[pROC]{pROC}} of the \pkg{pROC} package (by Xavier Robin et al.). It builds a ROC curve and returns a \code{"roc"} object, a list of class \code{"roc"}.
}
\usage{
ROC(x, resp = NULL, ...)
}
\arguments{
  \item{x}{a model object, or the predicted probabilities, when resp is not \code{NULL}. }
  \item{resp}{the response}
  \item{\dots}{all arguments are passed to \code{roc()}.
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
\author{Andri Signorell <andri@signorell.net>}


\references{
Peterson, A.T. et al. (2008) Rethinking receiver operating characteristic analysis applications in ecological niche modeling. Ecol. Modell., 213, 63-72.
}

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
\keyword{ model }
