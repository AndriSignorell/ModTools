\name{d.pima}
\alias{d.pima}
\alias{d.pima2}
\title{Diabetes survey on Pima Indians}
\description{The National Institute of Diabetes and
Digestive and Kidney Diseases conducted a study on 768 adult female
Pima Indians living near Phoenix.
}
\usage{
data(d.pima)
data(d.pima2)
}
\format{The dataset contains the following variables
  \describe{
    \item{\code{pregnant}}{ Number of times pregnant}
    \item{\code{glucose}}{ Plasma glucose concentration at 2 hours in an
      oral glucose tolerance test}
    \item{\code{pressure}}{ Diastolic blood pressure (mm Hg)}
    \item{\code{triceps}}{ Triceps skin fold thickness (mm)}
    \item{\code{insulin}}{ 2-Hour serum insulin (mu U/ml)}
    \item{\code{mass}}{ Body mass index (weight in kg/(height in metres
      squared))}
    \item{\code{pedigree}}{ Diabetes pedigree function}
    \item{\code{age}}{ Age (years)}
    \item{\code{diabetes}}{ test whether the patient shows signs of
      diabetes (coded \code{neg} if negative, \code{pos} if positive)}
  }
}
\details{\code{d.pima2} is the same dataset as \code{d.pima} with the only change, that invalid 0-values are replaced by \code{NA}s.}
\note{This dataset has been borrowed from Julian Faraway's package:\cr \emph{faraway}: Functions and datasets for books by Julian Faraway, 2015}
\source{ The data may also
be obtained from the package \code{MASS}.
}
\keyword{datasets}
