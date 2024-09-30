\name{Response}
\alias{Response}
\title{Extract the Response from Several Models
}
\description{Time after time, in the course of our daily work, we experience that the response variable is hidden very deeply in the object. This again leads to superfluous consultation of the documentation.
\code{Reponse()} relieves us of this work.
}
\usage{
Response(x, ...)
}
\arguments{
  \item{x}{the model to use
}
  \item{\dots}{more arguments
}
}
\details{The function implements the extraction of the response variables for all the models listed in the package's help text.
}
\value{the response of model x
}
\author{Andri Signorell <andri@signorell.net>
}


\seealso{\code{\link{model.frame}}, \code{\link{model.response}}, \code{\link{RefLevel}}
}
\examples{
r.rpart <- FitMod(diabetes ~ ., d.pima, fitfn="rpart")
Response(r.rpart)

# up to the attribute "response" this is the same
identical(StripAttr(Response(r.rpart), "response"),
          model.response(model.frame(r.rpart)))
}
\keyword{ models }
