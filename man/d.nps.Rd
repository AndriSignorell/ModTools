\name{d.nps}
\alias{d.nps}
\docType{data}
\title{Net Promoter Score
%%   ~~ data name/kind ... ~~
}
\description{The net promoter score is a fairly simple and yet meaningful customer loyalty and satisfaction measurement. It is taken from asking customers how likely they are to recommend a product or service to others (on a scale of 1-10, sometimes also from 0-10). \cr This data set is taken as a random sample from a survey on customer satisfaction conducted by 10 Swiss insurance companies and contains the likelihood of recommendation (the response variable we're normally interested in), but also other potential determinants as customer satisfaction, image assessment, and level of premiums.
}
\usage{data("d.nps")}
\format{
  A data frame with 8'963 observations on the following 15 variables.
  \describe{
    \item{\code{inr}}{a numeric vector, the interview number}
    \item{\code{resp_x}}{an ordered factor with levels \code{very unlikely (1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{very likely (10)}, the reponse to the question: \emph{"How likely is it, that you would recommend . to your best friend?"}}
    \item{\code{resp_n}}{the numeric value of the ordinal response \code{resp_x}}
    \item{\code{rcat_x}}{an ordered factor with levels \code{unlikely} < \code{undecided} < \code{likely}, the response, aggregated to 3 categories (unlikely=1:6, undecided=7:8, likely=9:10).}
    \item{\code{age_n}}{a numeric vector, the age in years}
    \item{\code{sex_c}}{the sex of the person, a factor with levels \code{male} \code{female}}
     \item{\code{region_x}}{the region of the person, a factor with levels \code{LEM} \code{MIT} \code{NW} \code{ZH} \code{OST} \code{ZS} \code{TI}}
    \item{\code{hins_x}}{code for the insurance company, a factor with levels \code{UHE} \code{FAS} \code{HAT} \code{WCO} \code{ECS} \code{RGR} \code{PKP} \code{ZSA} \code{BSW} \code{FVI}}
    \item{\code{image_x}}{an ordered factor with levels \code{disagree at all (1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{fully agree (10)}, \emph{"The company . has a good image."}}
    \item{\code{custsat_x}}{an ordered factor with levels \code{not satisfied at all (1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{fully satisfied (10)}, \emph{"How satisfied are you with company .?"}}
    \item{\code{lowprem_x}}{an ordered factor with levels \code{disagree at all (1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{fully agree (10)}, \emph{"The company . has low premiums."}}
    \item{\code{stabprem_x}}{an ordered factor with levels \code{disagree at all (1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{fully agree (10)}, \emph{"The company . has a stable premiums."}}

        \item{\code{durcon_x}}{an ordered factor with levels \code{[0-1)} < \code{[1-2)} < \code{[2-5)} < \code{[5-10)} < \code{[10-Inf)} < \code{don}'\code{t know}, duration of contract}
            \item{\code{income_x}}{an ordered factor with levels \code{[0-4500)} < \code{[4500-7000)} < \code{[7000-9000)} < \code{[9000-12000)} < \code{[12000-15000)} < \code{[15000-Inf)} < \code{no comment}, income per month in categories.}
    \item{\code{polit_x}}{an ordered factor with levels \code{(0) left} < \code{(1)} < \code{(2)} < \code{(3)} < \code{(4)} < \code{(5)} < \code{(6)} < \code{(7)} < \code{(8)} < \code{(9)} < \code{(10) right} < \code{don}'\code{t know}
    the relative political position from left to right.}
  }
}
\details{
The net promoter score describes the relationship between people who are rather positive about a brand and recommend it to others ("promoters") and people who do not ("detractors") and is calculated by subtracting the detractors (in percent of all respondents) from the promoters (in percent of all respondents). Detractors are people rating from 1-6, promoters from 9-10, the people between are called "passive".

An example: 54\% of the respondents state that they recommend a company or brand (scale 10+9). 34\% say that they are unlikely to do so (scale 0 to 6). In this case, the NPS of the company or brand would be 54\% - 34\% = 20\%. The higher the Net Promoter Score, the greater the attractiveness of a company or brand and the more loyal its customers are to it.
}
\source{own data
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\keyword{datasets}
