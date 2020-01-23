\name{SplitTrainTest}
\alias{SplitTrainTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Split DataFrame in Train an Test Sample
%%  ~~function to do ... ~~
}
\description{For modeling we usually split our data frame in a train sample, where we train our model on and a test sample, where we test, how good it works. This function splits a given data frame in two parts, one being the training sample and the other the test sample in form of a list with two elements.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SplitTrainTest(x, p = 0.1, seed = NULL, logical = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data.frame
%%     ~~Describe \code{x} here~~
}
  \item{p}{proportion for test sample. Default is 10\%.
%%     ~~Describe \code{p} here~~
}
  \item{seed}{initialization for random number generator.
}
\item{logical}{logical, defining if a logical vector should be returned or the list with train and test data. Default is \code{FALSE}. }
}
\details{In order to obtain reasonable models, we should ensure two points. The dataset must be large enough to yield statistically meaningful results and it should be representative of the data set as a whole. Assuming that our test set meets the preceding two conditions, our goal is to create a model that generalizes well to new data. We are aiming for a model that equally well predicts training and test data. We should never train on test data. If we are seeing surprisingly good results on the evaluation metrics, it might be a sign that we're accidentally training on the test set.
%%  ~~ If necessary, more details than the description above ~~
}
\value{depending on \code{logical} either \cr list with two entries named \code{train} and \code{test}
or \cr logical vector having \code{nrow} elements
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
}

\examples{
SplitTrainTest(d.pima)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{model }% use one of  RShowDoc("KEYWORDS")
