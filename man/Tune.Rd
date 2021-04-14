\name{Tune}
\alias{Tune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tune Classificators
%%  ~~function to do ... ~~
}
\description{Some classifiers benefit more from adjusted parameters to a particular dataset than others. However, it is often not clear from the beginning how the parameters have to be determined.  What often only remains is a grid search when several parameters have to be found in combination. The present function uses a grid search approch for the decisive arguments (typically for a neural network, a random forest or a classification tree). However it's not restricted to these models, any model fulfilling weak interface standards could be provided.
}
\usage{
Tune(x, ..., testset = NULL, keepmod = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the model to be tuned, best (but not necessarily) trained with \code{\link{FitMod}}.
}
  \item{\dots}{a list of parameters, containing the values to be used for a grid search.}

\item{testset}{a testset containing all variables required in the model to be used for calculating independently the accuracy (normally a subset of the original dataset).}
  \item{keepmod}{logical, defining if all fitted models should be returned in the result set. Default is \code{TRUE}. (Keep an eye on your RAM!)}
}
\details{The function creates a n-dimensional grid according to the given parameters and calculates the model with the combinations of all the parameters. The accuracy for the models are calculated insample and on a test set, if one has been provided.

It makes sense to avoid overfitting to provide a test set to also be evaluated.
A matrix with all combination of the values for the given parameters, sorted by accuracy, either by the accuracy achieved in the test set or the insample accuracy is returned.
}
\value{a matrix with all supplied parameters and a column \code{"acc"} and \code{"test_acc"} (if a test set has been provided) }

\author{Andri Signorell <andri@signorell.net>}


\examples{
d.pim <- SplitTrainTest(d.pima, p = 0.2)
mdiab <- formula(diabetes ~ pregnant + glucose + pressure + triceps
                 + insulin + mass + pedigree + age)

# tune a neural network for size and decay
r.nn <- FitMod(mdiab, data=d.pim$train, fitfn="nnet")
(tu <- Tune(r.nn, size=12:17, decay = 10^(-4:-1), testset=d.pim$test))

# tune a random forest
r.rf <- FitMod(mdiab, data=d.pim$train, fitfn="randomForest")
(tu <- Tune(r.rf, mtry=seq(2, 20, 2), testset=d.pim$test))

# tune a SVM model
r.svm <- FitMod(mdiab, data=d.pim$train, fitfn="svm")
tu <- Tune(r.svm,
           kernel=c("radial", "sigmoid"),
           cost=c(0.1,1,10,100,1000),
           gamma=c(0.5,1,2,3,4), testset=d.pim$test)

# let's get some more quality measures
tu$modpar$Sens <- sapply(tu$mods, Sens)     # Sensitivity
tu$modpar$Spec <- sapply(tu$mods, Spec)     # Specificity
Sort(tu$modpar, ord="test_acc", decreasing=TRUE)
}
