\name{Over-/Undersample}
\alias{OverSample}
\alias{UnderSample}

\title{Oversample and Undersample
}
\description{For classification purposes we might want to have balanced datasets. If the response variable has not a prevalence of 50\%, we can sample records for getting as much response A cases as response B. This is called oversample. Undersample means to sample the (lower) number of cases A from the records of case B.
}
\usage{
OverSample(x, vname)
UnderSample(x, vname)
}
\arguments{
  \item{x}{a data frame containing predictors and response
}
  \item{vname}{the name of the response variable to be used to over/undersample
}
}
\value{a data frame with balanced response variable
}
\author{Andri Signorell <andri@signorell.net>
}
\seealso{ \code{\link{BestCut}}
}
\examples{
OverSample(d.pima2, "diabetes")

UnderSample(d.pima2, "diabetes")
}
