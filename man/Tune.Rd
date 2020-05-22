\name{Tune}
\alias{Tune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tune Classificators
%%  ~~function to do ... ~~
}
\description{Several classificators all need their individual tuning parameters. This function uses a grid search approch for the decisive arguments or randomForest, rpart tree, neural networks, and so on.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Tune(x, scale = FALSE, sort = TRUE, type = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the model to be tuned
%%     ~~Describe \code{x} here~~
}
  \item{scale}{the scale
%%     ~~Describe \code{scale} here~~
}
  \item{sort}{should the result be sorted?
%%     ~~Describe \code{sort} here~~
}
  \item{type}{what kind of model type.
%%     ~~Describe \code{type} here~~
}
  \item{\dots}{ar passed on
%%     ~~Describe \code{\dots} here~~
}
}
\details{The details.
%%  ~~ If necessary, more details than the description above ~~
}
\value{the optimized parameters
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}
