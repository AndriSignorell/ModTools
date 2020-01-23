\name{predict.NaiveBayes}
\alias{predict.NaiveBayes}
\title{Predict Method for NaiveBayes Objects}
\usage{
\method{predict}{NaiveBayes}(object, newdata = NULL, type = c("class","prob"),
        na.action = na.pass, threshold = 0.001, ...)
}
\arguments{
\item{object}{object of class inheriting from \code{"NaiveBayes"}.}

\item{newdata}{numeric matrix or dataframe with categorical (character/factor/logical) or metric (numeric) predictors.}

\item{type}{if "class", new data points are classified according to the highest posterior probabilities. If "prob", the posterior probabilities for each class are returned.}

\item{threshold}{value by which zero probabilities corresponding to metric variables are replaced (zero probabilities corresponding to categorical variables can be handled with Laplace smoothing).}

\item{na.action}{a function to determine what should be done with missing values in \code{newdata}. The default is to pass them down the tree using surrogates in the way selected when the model was built. Other possibilities are \code{na.omit} and \code{na.fail}.}

\item{...}{not used.}
}
\value{
\code{predict.NaiveBayes} returns a factor with class corresponding to the maximal conditional posterior probability or a matrix with conditional posterior probabilities for each class.
}
\description{
Classification based on Naive Bayes models.
}
\details{
Computes conditional posterior probabilities for each class using the Bayes' rule. If no new data is provided, the data from the object is used. Logical variables are treated as categorical variables. Predictors with missing values are not included into the computation of posterior probabilities.
}
\examples{
ind_iris <- sample(1:nrow(iris), size = round(0.3 * nrow(iris)))
iris_train <- iris[-ind_iris, ]
iris_test <- iris[ind_iris, ]

nb_iris <- NaiveBayes(Species ~ ., iris_train)
predict(nb_iris, iris_test)
head(predict(nb_iris, iris_test, type = "prob"))


\dontrun{
vars <- 10
rows <- 500000
y <- sample(c("a", "b"), rows, TRUE)

# Only categorical variables
X1 <- as.data.frame(matrix(sample(letters[5:9], vars * rows, TRUE),
                           ncol = vars))
nb_cat <- NaiveBayes(x = X1, y = y)
nb_cat
system.time(pred2 <- predict(nb_cat, X1))

# Only numeric variables
X2 <- as.data.frame(matrix(rnorm(vars * rows), ncol = vars))
nb_num <- NaiveBayes(x = X2, y = y)
nb_num
system.time(pred2 <- predict(nb_num, X2))
}
}
\author{
Michal Majka, \email{michalmajka@hotmail.com}
}
\seealso{
\code{\link{NaiveBayes}}, \code{\link{plot.NaiveBayes}}, \code{\link[naivebayes]{tables}}
}

