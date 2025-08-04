

# print.multinom
# wie https://stats.oarc.ucla.edu/stata/output/multinomial-logistic-regression/
#
# library(haven)
# d.ice_cream <- as.data.frame(read_dta("C:/Users/HK1S0/Downloads/mlogit.dta"), stringsAsFactors=TRUE)
# d.ice_cream$ice_cream <- relevel(
#   factor(d.ice_cream$ice_cream,
#          labels=names(attr(d.ice_cream$ice_cream, "labels"))),
#   ref = "vanilla")
#
# r.mult <- multinom(ice_cream ~ video + puzzle + female,
#                    data = d.ice_cream, model=TRUE)
#
# d.ice_cream$ice_cream
#


summary.FitMod.multinom <- function(x, conf.level=NULL, ...) {

  if(is.null(conf.level)) conf.level <- 0.95
  alpha <- 1-(1-conf.level) / 2

  # class(x) <- class(x)[class(x)!="regr"]
  r.summary <- summary(x)

  est <- t(r.summary$coefficients)
  est <- reshape(data.frame(est, id=row.names(est)),
                 varying=1:ncol(est), idvar="id",
                 times=colnames(est), v.names="estimate", direction="long")

  se <- t(r.summary$standard.errors)
  se <- reshape(data.frame(se), varying=1:ncol(se),
                times=colnames(se), v.names="se", direction="long")[, "se"]

  z <- abs(est[, "estimate"]/se)

  d.res <- data.frame(
    est[,c("time","id","estimate")],
    stderr = se,
    z = z,
    "ci" = est[, "estimate"] - qnorm(alpha) * se,
    "uci" = est[, "estimate"] + qnorm(alpha) * se,
    "pval" = 2*(1-pnorm(q=z))
  )

  colnames(d.res)[6] <- gettextf("%s-lci", Fm(conf.level, fmt="%", digits=0))
  rownames(d.res) <- NULL

  res <- list(coefficients = d.res, call = x$call,
              nobs=NROW(x$residuals),
              na.action=length(x$na.action),
              # Probably not available:
              # BrierScore = BrierScore(x),
              PseudoR2 = PseudoR2(x, which="all"),
              response=c(attr(Response(x), "response"),
                         levels(Response(x))[1]),
              results = "Coefficients"
          )

  return(res)

}



print.FitMod.multinom  <- function(x, digits=3, pdigits=3,
                                   conf.level=0.95, ...){

  if(inherits(x, "SummaryFitMod"))
    xx <- x
  else
    xx <- summary.FitMod.multinom(x, conf.level=conf.level, ...)

  out <- cbind(
    " "=xx$coefficients$id,
    Fm(xx$coefficients[, c("estimate","95%-lci","uci")], digits=3),
    pval=Fm(xx$coefficients$pval, fmt="p", digits=3, eps=1e-3),
    " "=Fm(xx$coefficients$pval, fmt="*")
  )

  g <- which(!duplicated(xx$coefficients$time))
  out <- out[sort(c(seq(nrow(out)), g)), ]
  out[, 1] <- paste("  ", out[, 1], sep="")
  ti <- g + c(0, seq(length(g)-1))

  words <- c(xx$coefficients$time, out[, 1])
  out[ti, 1] <- StrAlign(c(words[which.max(nchar(words))],
                           unique(xx$coefficients$time)), sep="\\l")[-1]
  out[ti, -1] <- ""

  cat("\nCall:\n", paste(deparse(xx$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat(gettextf("\n%s:\n", xx$results))
  cat(gettextf("(%s == %s is the base outcome)\n\n",
               xx$response[1], xx$response[2]))

  print(out, quote=FALSE, print.gap=2, row.names=FALSE)

  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  cat(gettextf("\nObs (NAs): %s (%s)", xx$nobs, xx$na.action))

  cat("\tPseudo R\u00B2 (McFadden):", Fm(xx$PseudoR2["McFadden"], digits=3))
  cat("   ", "AIC:", Fm(xx$PseudoR2["AIC"], digits=3))
  cat("\n\n")

}



OddsRatio.multinom <- function(x, conf.level=NULL, digits=3, ...) {

  xx <- summary.FitMod.multinom(x, conf.level=conf.level, digits=digits, ...)

  xx$coefficients[c("estimate", "95%-lci", "uci")] <-
    exp(xx$coefficients[c("estimate", "95%-lci", "uci")])
  xx["results"] <- "OddsRatio"

  class(xx) <- c("SummaryFitMod", "FitMod", class(x))
  return(xx)

}




