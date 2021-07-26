\name{LogitBoost}
\alias{LogitBoost}
\alias{LogitBoost.default}
\alias{LogitBoost.formula}
\title{LogitBoost Classification Algorithm}

\description{Train logitboost classification algorithm using decision
  stumps (one node decision trees) as weak learners.  }

\usage{
LogitBoost(x, ...)

\method{LogitBoost}{formula}(formula, data, ..., subset, na.action)

\method{LogitBoost}{default}(x, y, nIter=ncol(x), ...)

}


\arguments{

  \item{formula}{a formula expression as for regression models, of the form \code{response ~ predictors}. The response should be a factor or a matrix with K columns, which will be interpreted as counts for each of K classes. See the documentation of \code{\link{formula}()} for other details.
%%     ~~Describe \code{formula} here~~
}
  \item{data}{an optional data frame in which to interpret the variables occurring in formula.
%%     ~~Describe \code{data} here~~
}
  \item{\dots}{additional arguments for nnet
%%     ~~Describe \code{\dots} here~~
}
  \item{subset}{expression saying which subset of the rows of the data should be used in the fit. All observations are included by default.
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{a function to filter missing data.
%%     ~~Describe \code{na.action} here~~
}

  \item{x}{A matrix or data frame with training data. Rows contain samples
    and columns contain features}
  \item{y}{Class labels for the training data samples.
    A response vector with one label for each row/component of \code{xlearn}.
    Can be either a factor, string or a numeric vector.}
  \item{nIter}{An integer, describing the number of iterations for
     which boosting should be run, or number of decision stumps that will be
     used.}
}

\details{
  The function was adapted from logitboost.R function written by Marcel
  Dettling. See references and "See Also" section. The code was modified in
  order to make it much faster for very large data sets. The speed-up was
  achieved by implementing a internal version of decision stump classifier
  instead of using calls to \code{\link[rpart]{rpart}}. That way, some of the most time
  consuming operations were precomputed once, instead of performing them at
  each iteration. Another difference is that training and testing phases of the
  classification process were split into separate functions.
}

\value{
  An object of class "LogitBoost" including components:
  \item{Stump}{List of decision stumps (one node decision trees) used:
    \itemize{
      \item column 1: feature numbers or each stump, or which column each stump
       operates on
      \item column 2: threshold to be used for that column
      \item column 3: bigger/smaller info: 1 means that if values in the column
       are above threshold than corresponding samples will be labeled as
       \code{lablist[1]}. Value "-1" means the opposite.
    }
    If there are more than two classes, than several "Stumps" will be
    \code{cbind}'ed
   }
  \item{lablist}{names of each class}
}

\references{
  Dettling and Buhlmann (2002), \emph{Boosting for Tumor Classification of Gene
    Expression Data}.

  \url{http://www.cs.princeton.edu/~schapire/boost.html}
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}}


\examples{
# basic interface
r.lb <- LogitBoost(Species ~ ., data=iris, nIter=20)
pred <- predict(r.lb)
prob <- predict(r.lb, type="prob")
d.res <- data.frame(pred, prob)
d.res[1:10, ]

# accuracy increases with nIter (at least for train set)
table(predict(r.lb, iris, type="class", nIter= 2), iris$Species)
table(predict(r.lb, iris, type="class", nIter=10), iris$Species)
table(predict(r.lb, iris, type="class"),           iris$Species)

# example of spliting the data into train and test set
d.set <- SplitTrainTest(iris)
r.lb <- LogitBoost(Species ~ ., data=d.set$train, nIter=10)
table(predict(r.lb, d.set$test, type="class", nIter=2), d.set$test$Species)
table(predict(r.lb, d.set$test, type="class"),          d.set$test$Species)
}

\keyword{classif}
