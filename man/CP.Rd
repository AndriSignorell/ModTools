\name{CP}
\alias{CP}
\alias{print.CP}
\alias{plot.CP}
\title{Complexity Parameter of an rpart Model
}
\description{Extracts, prints and plots the complexity table of an \code{rpart} model.
}
\usage{
CP(x, ...)

\method{print}{CP}(x, digits = getOption("digits") - 2L, ...)
\method{plot}{CP}(x, minline = TRUE, lty = 3, col = 1,
     upper = c("size", "splits", "none"), ...)
}
\arguments{
  \item{x}{fitted model object of class \code{"rpart"}. This is assumed to be the result of some function that produces an object with the same named components as that returned by the \code{rpart} function.
}
\item{digits}{the number of digits of numbers to print.}
\item{minline}{whether a horizontal line is drawn 1SE above the minimum of the curve.}
\item{lty}{line type for this line}
\item{col}{colour for this line}
\item{upper}{what is plotted on the top axis: the size of the tree (the number of leaves) ("\code{size}"), the number of splits ("\code{splits}") or nothing ("\code{none}").}

  \item{\dots}{further arguments passed to \code{print} and \code{plot}
}
}
\details{The complexity parameter table is hidden deep in the entrails of the \code{rpart} result object, it is convenient to have a function to extract it.
}
\value{A list containing the following components:
\item{cp }{the complexity table}
\item{x }{the \code{rpart} object}

}
\author{Andri Signorell <andri@signorell.net>
}
\seealso{\code{\link[rpart]{printcp}}, \code{\link[rpart]{plotcp}}
}
\examples{
r.rp <- FitMod(diabetes ~ ., d.pima, fitfn="rpart")

CP(r.rp)
plot(CP(r.rp))
}
\keyword{models }
