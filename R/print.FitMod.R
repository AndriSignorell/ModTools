

print.FitMod <- function(x, digits = 3, pdigits=3,
                         conf.level = 0.95, ...){

  if(identical(class(x)[-1], "lm")){
    print.FitMod.lm(x, digits, pdigits, conf.level, ...)

  } else   if(inherits(x, "glm")){
    print.FitMod.lm(x, digits, pdigits, conf.level, ...)

  } else if(inherits(x, "multinom")){
    print.FitMod.multinom(x, digits, pdigits, conf.level, ...)

  } else {
    # default
    class(x) <- class(x)[class(x) %nin% "FitMod"]
    print(x, ...)
  }

}



print.FitMod.lm <- function(x, digits = 3, pdigits=3,
                            conf.level = 0.95, ...) {

  n <- nobs(x)
  ci <- confint(x, level=conf.level)
  ref <- RefLevel(x)

  isGLM <- inherits(x, "glm")

  # anova_p <- drop1(x, test="F",
  #                  scope= formula(gettextf("~ %s",
  #                                          paste(names(ref),
  #                                                collapse="+"))))[names(ref), "Pr(>F)"]

  if(isGLM)
    anova_p <- drop1(x, test="Chisq")[names(ref), "Pr(>Chi)"]
  else
    anova_p <- drop1(x, test="F")[names(ref), "Pr(>F)"]


  xx <- summary(x)

  cat("\nCall:\n", paste(deparse(xx$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  resid <- xx$residuals
  df <- xx$df
  rdf <- df[2L]

  if (length(xx$aliased) == 0L) {
    cat("\nNo Coefficients\n")

  } else {
    if (nsingular <- df[3L] - df[1L])
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
          sep = "")
    else cat("\nCoefficients:\n")

    coefs <- xx$coefficients
    if (any(aliased <- xx$aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                              colnames(coefs)))
      coefs[!aliased, ] <- xx$coefficients
    }


    out <- cbind(Format(cbind(coefs[, 1], ci), digits=digits),
                 gsub(" ", "", Format(coefs[, 4], fmt="p", eps = 10^-pdigits, digits=pdigits)), Format(coefs[, 4], fmt="*"))
    colnames(out) <- c("estimate",
                       gettextf(c("%s-lci", "uci"), Format(conf.level, fmt="%",
                                                           digits=max(0, Ndec(as.character(signif(conf.level)))-2))),
                       "p-val", "")

    for(i in seq(ref)){
      # rnr <- grep(gettextf("^%s", names(ref)[i]), rownames(out))[1]
      rnr <- grep(gettextf("^%s", gsub("[^a-zA-Z0-9_]", " ", names(ref)[i])),
                  gsub("[^a-zA-Z0-9_]", " ", rownames(out)))[1]

      p <- anova_p[i]

      out <- Append(out, rbind(c(rep(".", 3),
                                 gettextf("%s", gsub(" ", "", Format(p, fmt="p", eps=10^-pdigits, digits=pdigits))),
                                 Format(p, fmt="*") )),
                    after = rnr-1, rows = TRUE)
      rownames(out)[rnr] <- paste(names(ref)[i], "(ref: ", ref[i], ")", sep="")
    }

    for(i in seq(ref)){
      rnr <- grep(gettextf("^%s", names(ref)[i]), rownames(out))
      rownames(out)[rnr] <- gsub(names(ref)[i], paste(names(ref)[i], " ", sep=""), rownames(out)[rnr])
    }


    print(out, quote=FALSE, right=TRUE, print.gap=2)

    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  }


  cat(gettextf("\nObs (NAs): %s (%s)", nobs(x), length(x$na.action)))
  if (!is.null(xx$fstatistic)) {
    cat("\tR\u00B2/R\u00B2adj:", paste(formatC(c(xx$r.squared, xx$adj.r.squared), digits = digits), collapse="/"))
  } else if(inherits(x, "glm")){
    cat("\tPseudo R\u00B2 (McFadden):", Fm(PseudoR2(x)["McFadden"], digits=3))
  }

  cat("   ", "AIC:", AIC(x))
  cat("\n\n")

  invisible(xx)

}

