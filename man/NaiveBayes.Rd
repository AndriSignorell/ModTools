\name{NaiveBayes}
\alias{NaiveBayes}
\alias{NaiveBayes.default}
\alias{NaiveBayes.formula}
\title{Naive Bayes Classifier}
\usage{
\method{NaiveBayes}{default}(x, y, prior = NULL, laplace = 0,
  usekernel = FALSE, ...)

\method{NaiveBayes}{formula}(formula, data, prior = NULL, laplace = 0,
  usekernel = FALSE, ..., subset, na.action)
}

\arguments{
\item{x}{numeric matrix or dataframe with categorical (character/factor/logical) or metric (numeric) predictors.}

\item{y}{class vector (character/factor/logical).}

\item{formula}{an object of class \code{"formula"} (or one that can be coerced to "formula") of the form \code{class ~ predictors} (class has to be a factor/character/logical).}

\item{data}{numeric matrix or dataframe with categorical (character/factor/logical) or metric (numeric) predictors.}

\item{laplace}{value used for Laplace smoothing. Defaults to 0 (no Laplace smoothing).}

\item{usekernel}{logical; if \code{TRUE}, \code{\link[stats]{density}} is used to estimate the densities of metric predictors.}

\item{prior}{vector with prior probabilities of the classes. If unspecified, the class proportions for the training set are used. If present, the probabilities should be specified in the order of the factor levels.}


\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}

\item{na.action}{a function which indicates what should happen when the data contain \code{NAs}. By default (\code{na.pass}), missing values are not removed from the data and are then omited while constructing tables. Alternatively, \code{na.omit} can be used to exclude rows with at least one missing value before constructing tables.}

\item{...}{other graphical parameters to \code{\link[stats]{density}} (\code{na.rm} defaults to \code{TRUE}).}

}
\value{
\code{NaiveBayes} returns an object of class \code{"NaiveBayes"} which is a list with following components:
    \item{data}{list with two components: \code{x} (dataframe with predictors) and \code{y} (class variable).}
    \item{tables}{list of tables. For each categorical predictor a table with class-conditional probabilities and for each metric predictor a table with means and standard deviations or \code{\link[stats]{density}} objects for each class.}
    \item{prior}{numeric vector with prior probabilities.}
    \item{laplace}{value of Laplace smoothing.}
    \item{levels}{character vector with values of the class variable.}
    \item{usekernel}{logical; \code{TRUE}, if the density was used for estimating densities of numeric variables.}
    \item{call}{the call that produced this object.}
}
\description{
\code{NaiveBayes} is used to fit Naive Bayes model in which predictors are assumed to be independent within each class label.
}
\details{
Metric predictors are handled by assuming that they follow Gaussian distribution, given the class label. Alternatively, kernel density estimation can be used to estimate their class-conditional distributions. Missing values are not included into constructing tables. Logical variables are treated as categorical variables.
}
\examples{
data(iris)
nb <- NaiveBayes(Species ~ ., data = iris)
plot(nb)

nb_kernel <- NaiveBayes(x = iris[-5], y = iris[ ,5], usekernel = TRUE)
plot(nb_kernel)


\dontrun{
vars <- 10
rows <- 500000
y <- sample(c("a", "b"), rows, TRUE)

# Only categorical variables
X1 <- as.data.frame(matrix(sample(letters[5:9], vars * rows, TRUE),
                           ncol = vars))
nb_cat <- FitMod(x = X1, y = y, fitfn="naiveBayes")
nb_cat
system.time(pred2 <- predict(nb_cat, X1))

# Only numeric variables
X2 <- as.data.frame(matrix(rnorm(vars * rows), ncol = vars))
nb_num <- FitMod(x = X2, y = y, fitfn="naiveBayes")
nb_num
system.time(pred2 <- predict(nb_num, X2))
}
}
\author{
Michal Majka, \email{michalmajka@hotmail.com}
}
\seealso{
\code{\link{predict.NaiveBayes}}, \code{\link{plot.NaiveBayes}}
}

