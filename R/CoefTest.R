


CoefTest <- function(x, vcov. = NULL, df = NULL, ...) {
  UseMethod("CoefTest")
}


CoefTest.default <- function(x, vcov. = NULL, df = NULL, ...) {

  # ## use S4 methods if loaded
  # coef0 <- if("stats4" %in% loadedNamespaces()) stats4::coef else coef
  # vcov0 <- if("stats4" %in% loadedNamespaces()) stats4::vcov else vcov
  coef0 <- coef
  vcov0 <- vcov

  ## extract coefficients and standard errors
  est <- coef0(x)
  if(is.null(vcov.)) se <- vcov0(x) else {
    if(is.function(vcov.)) se <- vcov.(x)
    else se <- vcov.
  }
  se <- sqrt(diag(se))

  ## match using names and compute t/z statistics
  if(!is.null(names(est)) && !is.null(names(se))) {
    anames <- names(est)[names(est) %in% names(se)]
    est <- est[anames]
    se <- se[anames]
  }
  tval <- as.vector(est)/se

  ## apply central limit theorem
  if(is.null(df)) {
    df <- try(df.residual(x), silent = TRUE)
    if(inherits(df, "try-error")) df <- NULL
  }
  if(is.null(df)) df <- 0

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "CoefTest"
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  ##  dQuote(class(x)[1]), "object", sQuote(deparse(substitute(x))))
  return(rval)
}

CoefTest.glm <- function(x, vcov. = NULL, df = Inf, ...)
  CoefTest.default(x, vcov. = vcov., df = df, ...)

CoefTest.mlm <- function(x, vcov. = NULL, df = NULL, ...) {

  ## obtain vcov
  v <- if(is.null(vcov.)) vcov(x) else if(is.function(vcov.)) vcov.(x) else vcov.

  ## nasty hack: replace coefficients so that their names match the vcov() method
  x$coefficients <- structure(as.vector(x$coefficients), .Names = colnames(vcov(x)))

  ## call default method
  CoefTest.default(x, vcov. = v, df = df, ...)
}

CoefTest.survreg <- function(x, vcov. = NULL, df = Inf, ...) {

  if(is.null(vcov.)) v <- vcov(x) else {
    if(is.function(vcov.)) v <- vcov.(x)
    else v <- vcov.
  }
  if(length(x$coefficients) < NROW(x$var)) {
    x$coefficients <- c(x$coefficients, "Log(scale)" = log(x$scale))
  }
  CoefTest.default(x, vcov. = v, df = df, ...)
}


CoefTest.multinom <- function(x, vcov. = NULL, df = NULL, ...)
{
  ## extract coefficients
  est <- coef(x)
  if(!is.null(dim(est))) {
    est <- structure(as.vector(t(est)),
                     names = as.vector(t(outer(rownames(est), colnames(est), paste, sep = ":"))))
  }

  ## process vcov.
  if(is.null(vcov.)) vc <- vcov(x) else {
    if(is.function(vcov.)) vc <- vcov.(x)
    else vc <- vcov.
  }
  se <- sqrt(diag(vc))
  tval <- as.vector(est)/se

  ## process degrees of freedom
  if(is.null(df)) df <- Inf

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  return(rval)
}

CoefTest.polr <- function(x, vcov. = NULL, df = NULL, ...)
{
  ## extract coefficients
  est <- c(x$coefficients, x$zeta)

  ## process vcov.
  if(is.null(vcov.)) vc <- vcov(x) else {
    if(is.function(vcov.)) vc <- vcov.(x)
    else vc <- vcov.
  }
  se <- sqrt(diag(vc))
  tval <- as.vector(est)/se

  ## process degrees of freedom
  if(is.null(df)) df <- Inf

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  return(rval)
}

lrtest.fitdistr <- function(object, ..., name = NULL)
{
  if(is.null(name)) name <- function(x) if(is.null(names(x$estimate))) {
    paste(round(x$estimate, digits = max(getOption("digits") - 3, 2)), collapse = ", ")
  } else {
    paste(names(x$estimate), "=", round(x$estimate, digits = max(getOption("digits") - 3, 2)), collapse = ", ")
  }
  lmtest::lrtest.default(object, ..., name = name)
}



# CoefTest.breakpointsfull <- function(x, vcov. = NULL, df = NULL, ...) {
#
#   stopifnot(require("strucchange"))
#
#   est <- coef(x, ...)
#
#   if(is.null(df)) {
#     df <- df.residual(x, ...)
#     df <- as.vector(rep(df, rep(NCOL(est), length(df))))
#   }
#
#   rnames <- as.vector(t(outer(rownames(est), colnames(est), paste)))
#   est <- as.vector(t(est))
#
#   se <- vcov(x, vcov. = vcov., ...)
#
#   se <- as.vector(sapply(seq(along = se), function(x) sqrt(diag(se[[x]]))))
#   tval <- est/se
#
#   if(any(is.finite(df) && df > 0)) {
#     pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
#     cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
#     mthd <- "t"
#   } else {
#     pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
#     cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
#     mthd <- "z"
#   }
#   rval <- cbind(est, se, tval, pval)
#   colnames(rval) <- cnames
#   rownames(rval) <- rnames
#   class(rval) <- "CoefTest"
#   attr(rval, "method") <- paste(mthd, "test of coefficients")
#   ##  dQuote(class(x)[1]), "object", sQuote(deparse(substitute(x))))
#   return(rval)
# }


print.CoefTest <- function(x, ...) {
  mthd <- attr(x, "method")
  if(is.null(mthd)) mthd <- "Test of coefficients"
  cat(paste("\n", mthd,":\n\n", sep = ""))
  printCoefmat(x, ...)
  cat("\n")
  invisible(x)
}


