\name{SplitTrainTest}
\alias{SplitTrainTest}
\title{Split DataFrame in Train an Test Sample
}
\description{For modeling we usually split our data frame in a train sample, where we train our model on, and a test sample, where we test, how good it works. This function splits a given data frame in two parts, one being the training sample and the other the test sample in form of a list with two elements.
}
\usage{
SplitTrainTest(x, p = 0.1, seed = NULL, logical = FALSE)
}
\arguments{
  \item{x}{data.frame
}
  \item{p}{proportion for test sample. Default is 10\%.
}
  \item{seed}{initialization for random number generator.
}
\item{logical}{logical, defining if a logical vector should be returned or the list with train and test data. Default is \code{FALSE}. }
}
\details{In order to obtain reasonable models, we should ensure two points. The dataset must be large enough to yield statistically meaningful results and it should be representative of the data set as a whole. Assuming that our test set meets the preceding two conditions, our goal is to create a model that generalizes well to new data. We are aiming for a model that equally well predicts training and test data. We should never train on test data. If we are seeing surprisingly good results on the evaluation metrics, it might be a sign that we're accidentally training on the test set.
}
\value{If \code{logical} is \code{FALSE} a list with two data frames, \code{train} and \code{test}, of the same structure as the given data in \code{x}
\cr if \code{logical} is \code{TRUE} a logical vector containing \code{nrow} elements of \code{TRUE} and \code{FALSE}
}
\author{Andri Signorell <andri@signorell.net>
}

\examples{
SplitTrainTest(d.pima)
}
\keyword{model }
