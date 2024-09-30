\name{TModC}
\alias{TModC}
\alias{plot.TModC}
\title{Compare Classification Models
}
\description{For the comparison of several classification models, the AUC values and BrierScore values of the models are determined and tabulated. Both the absolute values and the relative values are reported, each related to the model with the highest corresponding value.
}
\usage{
TModC(..., newdata = NULL, reference = NULL, ord = NULL)

\method{plot}{TModC}(x, col = NULL, args.legend = NULL,...)
}
\arguments{
  \item{\dots}{the models to be compared}
  \item{x}{\code{TModC} object to plot}

  \item{newdata}{the data to use for predicting. If not provided, the \code{model.frame} will be used.}

  \item{reference}{the reference values}
  \item{ord}{character defining the order of the results table, can be any of \code{"auc"}, \code{"bs"},
  \code{"auc_p"}, \code{"bs_p"}, \code{"bs_rnk"}, \code{"auc_rnk"}, \code{"ensemble"} (using the mean of \code{"auc_p"} and \code{"bs_p"} for the ranking).}
  \item{col}{the color for the lines in the ROC plot}
  \item{args.legend}{the legend to be placed in the ROC plot}

}
\value{a matrix with the columns
  \item{auc }{absolute value of area under the ROC curve (AUC) }
  \item{auc_p}{percentage of the auc based on the best observerd AUC}
  \item{auc_rnk}{the rank of the auc}
  \item{bs }{absolute value of the Brier score}
  \item{bs_p}{percentage of the Brier score based on the best observed BS}
  \item{bs_rnk}{the rank of the BS}
  \item{auc_grnk}{character representation of the AUC rank}
  \item{bs_grnk}{character representation of the BS rank}
}

\author{Andri Signorell <andri@signorell.net>
}
\seealso{\code{\link[DescTools]{TMod}}, \code{\link[DescTools]{BrierScore}}, \code{\link[DescTools]{AUC}}, \code{\link{ROC}}
}
\examples{
d.pim <- SplitTrainTest(d.pima, p = 0.2)
mdiab <- formula(diabetes ~ pregnant + glucose + pressure + triceps +
                            insulin + mass + pedigree + age)

r.glm <- FitMod(mdiab, data=d.pim$train, fitfn="logit")
r.rp <- FitMod(mdiab, data=d.pim$train, fitfn="rpart")
mods <- list(glm=r.glm, rp=r.rp)

# the table with the measures
(tm <- TModC(mods, ord="auc"))

# plotting the ROC curves
plot(tm, col=c("darkmagenta", "dodgerblue"))
}
\keyword{models}
