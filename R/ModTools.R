# ****************************************************************************
#
# Projekt:	ModTools.r
#
# Zweck:    Modeling Tools for regression and classification
#
# Autor:    Andri Signorell
# Version:	0.5 (in development)
#
# Depends:  DescTools
#
# Datum:
#           21.10.2018  komplettes Redesign (aber ohne Check-Errors...)
#           04.12.2016  mit neuen DescTools
#           02.11.2013  version 0.3
#           20.09.2013 	created
#
# ****************************************************************************

# TODO (Andri):
#   Vergleich f?r lines in den Plot zeichnen etc.

# bereits in DescTools (weil basics...)
# statistics:
#   Brier, c-statistic, pseudoR2 (McFadden etc.)
#   DurbinWatson, HosmerLemeshow, BreuschGodfrey
#   findCorrelation, Conf, VIF, MAE etc., ModTable

# modelling:
#   linear regression, rlm, glm, tobit, cox ph, quantreg
# mixed models, gam, mars, lasso,
# discriminant analysis,
#
# rpart, randomForest, svm, NaiveBayes, nnet, C5.0
# knn, logitBoost, adaBoost
# (http://www.inf.fu-berlin.de/inst/ag-ki/adaboost4.pdf)
#
# unsupervised:
#   SOM, ART Adaptive Resonance Theory
# cluster, association rules,
# multidimensional scaling
#
# tuning:
#   finding the appropriate parameters
#
# preprocess
#
# 1. Check for near zero variance predictors and flag as near zero if:
#   * the percentage of unique values is less than 20
# * the ratio of the most frequent to the second most frequent
# value is greater than 20,
# 2. Check for susceptibility to multicollinearity
# * Calculate correlation matrix
# * Find variables with correlation 0.9 or more and delete them
#
#
# transformations, retransformation:
#   rcs - restricted cubic splines
# smearingEstimator
#
# variable selection:
#   lasso, stepAIC
#
# imputation:
#   ... ?
# knnImputation
#
# evaluation, verfication of models:
#   CrossValidation
# variable importance

# ROC and all what s around
#
# tests:
# plots:
#   regression plots, termplots, residuals etc.
# tree plots
# interactive exploring of a tree
#
# reporting:
#   tables for regression, OR
# tables for comparing models
#
# ChemometricsWithR:
#    GA Genetic Algorithms for variable selection in classification

# imports from libraries:
# randomForest:
#   "tuneRF", "varImpPlot", "importance" from randomForest wrapped into Tune, VarImp, plot(VarImp)
#
# , "tune"


Tobit <- AER::tobit


SplitTrainTest <- function(x, p=0.1, seed=NULL, logical=FALSE){

  if(!is.null(seed)) set.seed(seed)
  n <- ifelse(is.atomic(x), length(x), nrow(x))
  idx <- 1:n
  itest <- sample(idx, trunc(length(idx) * p))
  if(logical)
    res <- Unwhich(itest, n)
  else {
    if(is.atomic(x))
      res <- list(train=x[idx[-itest]], test=x[idx[itest]])
    else
      res <- list(train=x[idx[-itest], ], test=x[idx[itest], ])

  }

  return(res)

}




.print.multinom <- function (x, digits = x$digits, ...) {

  object <- x

  vc <- vcov(object)
  r <- length(object$vcoefnames)
  se <- sqrt(diag(vc))
  if (length(object$lev) == 2L) {
    coef <- object$wts[1L + (1L:r)]
    stderr <- se
    names(coef) <- names(stderr) <- object$vcoefnames

  } else {
    coef <- matrix(object$wts, nrow = object$n[3L], byrow = TRUE)[-1L,
                                                                  1L + (1L:r), drop = FALSE]
    stderr <- matrix(se, nrow = object$n[3L] - 1L, byrow = TRUE)
    if (length(l <- object$lab) || length(l <- object$lev))
      dimnames(coef) <- dimnames(stderr) <- list(l[-1L],
                                                 object$vcoefnames)
  }

  object$is.binomial <- (length(object$lev) == 2L)
  object$digits <- digits
  object$coefficients <- coef
  object$standard.errors <- stderr
  object$Wald.ratios <- coef/stderr
  # if (correlation)
  #   object$correlation <- vc/outer(se, se)
  # class(object) <- "summary.multinom"


  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl, control = NULL)
  }

  cat("\n")
  print(Format(x$drop1[,-3], digits=c(0,1,3), fmt=c("","","p*"), na.form = ""), print.gap = 3)
  cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("\nCoefficients:\n Level: ", x$lev[1], " (reference)\n\n")

  ci <- confint(x)
  b <- t(coef(x))
  lst <- list()
  for(i in seq(dim(ci)[3])){
    lst[[i]] <- cbind(b[, i], ci[,,i])
    # rownames(lst[[i]]) <- paste(colnames(b)[i], rownames(lst[[i]]) , sep=": ")
  }

  beta <- do.call(rbind, rbind(NA, lst, NA))

  betaor <- cbind(beta, exp(beta))
  colnames(betaor) <- c("coef", "lci", "uci", "or", "or.lci", "or.uci")
  rownames(betaor) <- paste0("  ", rownames(betaor))
  rownames(betaor)[seq(1, nrow(betaor), by=nrow(lst[[1]])+2)] <- paste0("Level: ", colnames(b))
  betaor[seq(2, nrow(betaor), by=nrow(lst[[1]])+2), 4:6] <- 0
  print(Format(betaor, digits=3, na.form = "", zero.form = "."), print.gap=2)

  cat("Residual Deviance: ", format(x$deviance), ",  ",
      "AIC: ", format(x$AIC), ",  nobs: ", nrow(x$residuals), "\n\n", sep="")
  if (!is.null(correl <- x$correlation)) {
    p <- dim(correl)[2L]
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1L, -p], quote = FALSE, ...)
    }
  }
  invisible(x)
}


.print.polr <- function(x, digits = x$digits, ...){

  if (!is.null(cl <- x$call)) {
    cat("\nCall:\n")
    dput(cl, control = NULL)
  }

  cat("\n")

  xdrop <- data.frame(x$drop1[, 3], x$drop1[,1], x$drop1[,4], x$drop1[,4])
  rownames(xdrop) <- rownames(x$drop1)
  colnames(xdrop) <- c("Df", "AIC",  "pval", "")

  print(Format(xdrop, fmt=c("", "", "p", "*"),
               digits=c(0,1,4,NULL)),
        print.gap = 3)

  cat("---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")

  cat("\nCoefficients:\n")

  tval <- x$coefficients / head(sqrt(diag(vcov(x))), -length(x$zeta))
  pval <- pnorm(abs(tval), lower.tail = FALSE) * 2

  res <- data.frame(coef=x$coefficients, x$ci,
                    exp(x$coefficients), exp(x$ci),
                    Format(pval, fmt="p"), Format(pval, fmt="*"))
  colnames(res) <- c("coef","lci","uci", "or", "or.lci", "or.uci", "p.val", "")

  res[, 1:6] <- Format(res[, 1:6], digits=digits)
  res <- as.matrix(res)

  # insert reference levels for factors
  rr <- RefLevel(x)$simple
  i <- 2
  for(i in seq_along(rr)){
    ff <- names(rr[i])
    rl <- rr[i]
    j <- head(grep(ff, rownames(res)), 1)
    res <- Append(res, "-", after=j-1, rows=TRUE)
    rownames(res)[grep(ff, rownames(res))] <-
      gsub(ff, "  ", rownames(res)[grep(ff, rownames(res))])
    rownames(res)[j] <- gettextf("%s (Ref: %s)", ff, rl)
    res[j, 8] <- Format(xdrop[grep(ff, rownames(xdrop)), 4], fmt="*")

  }
  res <- as.data.frame(res)
  colnames(res)[8] <- ""
  print(res, print.gap = 2)

  cat("\nIntercepts:\n")
  print(x$zeta)

  cat("\n")

}



.predict.nnet <- function (object, newdata, type = c("class", "prob"), cutoff=0.5, ...) {

  type <- match.arg(type)

  if (!requireNamespace("nnet", quietly = TRUE))
    stop("package 'nnet' must be installed")
  predict.nnet <- get("predict.nnet", asNamespace("nnet"),
                          inherits = FALSE)

  Y <- predict.nnet(object, newdata, type="raw", ...)

  if(type=="prob" && dim(Y)[2]==1){
    # make compatible to the others
    Y <- cbind(1-Y, Y)
    colnames(Y) <- object$lev
  }


  switch(type, class = {
    if (length(object$lev) > 2L)  Y <- factor(max.col(Y),
                                              levels = seq_along(object$lev), labels = object$lev)
    if (length(object$lev) == 2L) Y <- factor((Y > cutoff) + 1,
                                              levels = 1L:2L, labels = object$lev)
    if (length(object$lev) == 0L) Y <- factor(max.col(Y),
                                              levels = seq_along(object$lab), labels = object$lab)
  }, prob = {
  })


  return(Y)

}


.predict.multinom <- function (object, newdata, type = c("class", "prob"), ...) {

  type <- match.arg(type)

  if (!requireNamespace("nnet", quietly = TRUE))
    stop("package 'nnet' must be installed")
  predict.multinom <- get("predict.multinom", asNamespace("nnet"),
                          inherits = FALSE)
  res <- predict.multinom(object, newdata, type, ...)

  # extend here with the complementary probability

  if(type=="prob" && (is.vector(res) || dim(res)[2]==1)){
    res <- cbind(1-apply(cbind(res), 1, sum), cbind(res))
    if(is.vector(res))
      names(res) <- object$lev
    else
      colnames(res) <- object$lev
  }
  res

}




.predict.glm <- function(object, newdata, type = c("class", "prob", "link", "terms", "response"), cutoff = 0.5, ...) {

  type <- match.arg(type)

  object$lev <- levels(model.response(model.frame(object)))

  switch(type[1],
         "prob"={
           Y <- stats::predict.glm(object, newdata, type="response", ... )
           if(is.list(Y)){
             Y[["fit"]] <- cbind(1-Y[["fit"]], Y[["fit"]])
             colnames(Y[["fit"]]) <- object$lev
           } else {
             Y <- cbind(1-Y, Y)
             colnames(Y) <- object$lev
           }
         },

         "class"={
           Y <- (stats::predict.glm(object, newdata, type="response", ...) > cutoff) + 1
         },
         "link"={
           Y <- stats::predict.glm(object, newdata, type="link", ... )
         },
         "terms"={
           Y <- stats::predict.glm(object, newdata, type="terms", ... )
         },
         "response"={
           Y <- stats::predict.glm(object, newdata, type="response", ... )
         }
  )


  switch(type[1], class = {
    if (length(object$lev) > 2L)  Y <- factor(max.col(Y),
                                              levels = seq_along(object$lev), labels = object$lev)
    if (length(object$lev) == 2L) Y <- factor(Y,
                                              levels = 1L:2L, labels = object$lev)
    if (length(object$lev) == 0L) Y <- factor(max.col(Y),
                                              levels = seq_along(object$lab), labels = object$lab)
  }, prob = {
  })

  return(Y)

}




.predict.lda <- function(object, newdata, type = c("class", "prob", "term"),
                        prior = object$prior, dimen, method = c("plug-in", "predictive", "debiased"),
                        ...) {

  type <- match.arg(type)

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("package 'MASS' must be installed")
  predict.lda <- get("predict.lda", asNamespace("MASS"),
                      inherits = FALSE)

  res <- predict.lda(object, newdata, prior=prior, dimen=dimen, method=method, ...)

  switch(type,
         "prob"={
           res <- res$posterior
         },

         "class"={
           res <- res$class
         },
         "term"={
           res <- res$x
         }
  )

  return(res)

}





# .predict.naiveBayes <- function(object, newdata, type = c("class", "prob", "term")
#                                 , threshold = 0.001, eps = 0, ...) {
#
#   type <- match.arg(type)
#
#   if (!requireNamespace("e1071", quietly = TRUE))
#     stop("package 'e1071' must be installed")
#   predict.naiveBayes <- get("predict.naiveBayes", asNamespace("e1071"),
#                      inherits = FALSE)
#
#   if(missing(newdata))
#     newdata <- object$data$x
#   else
#     newdata <- as.data.frame(newdata)
#
#
#   if(type=="prob")
#     type <- "raw"
#   res <- predict.naiveBayes(object, newdata, type=type, threshold=threshold, eps=eps, ...)
#
#   return(res)
#
# }



.predict.qda <- function(object, newdata, type = c("class", "prob", "term"),
                        prior = object$prior, dimen, method = c("plug-in", "predictive", "debiased", "looCV"),
                        ...) {

  type <- match.arg(type)

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("package 'MASS' must be installed")
  predict.qda <- get("predict.qda", asNamespace("MASS"), inherits = FALSE)

  res <- predict.qda(object, newdata, prior=prior, dimen=dimen, method=method, ...)

  switch(type,
         "prob"={
           res <- res$posterior
         },

         "class"={
           res <- res$class
         },
         "term"={
           res <- res$x
         }
  )

  return(res)

}


.predict.svm <- function(object, newdata, type = c("class", "prob", "term"),
                         decision.values = TRUE, probability = TRUE ,
                         na.action = na.omit, ...) {

  type <- match.arg(type)

  if (!requireNamespace("e1071", quietly = TRUE))
    stop("package 'e1071' must be installed")
  predict.svm <- get("predict.svm", asNamespace("e1071"),
                     inherits = FALSE)

  if(missing(newdata)) {
    # newdata <- eval(object$call$data) # [, as.character(attr(object$terms, "predvars"))]
    # select all predictors: all.vars(attr(r.lm$model, "terms"))
    newdata <- eval(object$call$data)[, all.vars(object$terms)]
  }

  res <- predict.svm(object, newdata, decision.values=decision.values, probability=probability,
                     na.action=na.action, ...)

  switch(type,
         "prob"={
           res <- attr(res, "probabilities")
           if(object$nclasses == 2){
             res <- res[, c(2, 1)]
           }
         },

         "class"={
           res <- res
         },
         "term"={
           res <- attr(res, "decision.values")
         }
  )

  return(res)

}



.predict.C5.0 <- function(object, newdata, type = c("class","prob"),
                          na.action=na.omit, ...) {

  type <- match.arg(type)

  if(missing(newdata)) {
    newdata <- eval(object$call$data)[, all.vars(object$Terms)]
  }

  # res <- predict.svm(object, newdata, decision.values=decision.values, probability=probability,
  #                    na.action=na.action, ...)

  C50::predict.C5.0(object, newdata, type = type, na.action=na.action, ...)
}





OddsRatio.FitMod <- function (x, conf.level = NULL, ...){
  class(x) <- class(x)[class(x) != "FitMod"]
  NextMethod(x, conf.level=conf.level, ...)
}





ROC <- function (x, resp = NULL, ...) {

  if(is.null(resp))
    roc(predictor = predict(x, type="prob")[, 2], response = Response(x), plot=FALSE, ...)
  else
    roc(predictor=x, response = resp, plot=FALSE, ...)

}


BestCut <- function(x, method=c("youden", "closest.topleft")){
  pROC::coords(roc=x, x="best", best.method=method,
               transpose=TRUE)
}


confint.ROC <- function(object, parm, level = 0.95, ...) {
  pROC::ci.coords(roc=object, conf.level = level, ...)
}



FitMod <- function(formula, data, ..., subset, na.action=na.pass, fitfn=NULL){

  cl <- match.call()

  if(is.null(fitfn)){
    # guess fitting function if not provided, based on type of response variable

    # this would return the response variable
    # resp <- eval(cl$data, parent.frame())[, all.vars(cl$formula)[1]]
    # see also:  https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object

    # evaluate left hand side in order to trap transformations
    resp <- eval(parse(text=cl$formula[2]), envir = eval(cl$data, parent.frame()))

    if(IsDichotomous(resp))
      cl[["fitfn"]] <- fitfn <- "logistic"

    else if(inherits(resp, "factor"))
      if(inherits(resp, "ordered"))
        cl[["fitfn"]] <- fitfn <- "polr"
      else
        cl[["fitfn"]] <- fitfn <- "multinom"

    else if(inherits(resp, "integer"))
      cl[["fitfn"]] <- fitfn <- "poisson"

    else if(inherits(resp, "numeric"))
      cl[["fitfn"]] <- fitfn <- "lm"

    else {
      cl[["fitfn"]] <- fitfn <- "randomForest"
      warning("fitting function could not be reasonably guessed, consider providing fitfn")
    }

  }


  if(fitfn == "lm"){
    fun <- get("lm", asNamespace("stats"), inherits = FALSE)

  } else if(fitfn == "poisson"){
    fun <- get("glm", asNamespace("stats"), inherits = FALSE)
    cl[["family"]] <- "poisson"

  } else if(fitfn == "quasipoisson"){
    fun <- get("glm", asNamespace("stats"), inherits = FALSE)
    cl[["family"]] <- "quasipoisson"

  } else if(fitfn == "gamma"){
    fun <- get("glm", asNamespace("stats"), inherits = FALSE)
    if(is.null(cl[["link"]]))
      cl[["family"]] <- "Gamma"
    else {
      cl[["family"]] <- eval(parse(text = gettextf("function() Gamma(link=%s)", cl[["link"]])), parent.frame())
      cl[["link"]] <- NULL
    }

  } else if(fitfn == "negbin"){
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("package 'MASS' must be installed")
    if(is.null(cl[["theta"]]))
      # no theta available, use glm.nb and fit theta
      fun <- get("glm.nb", asNamespace("MASS"), inherits = FALSE)
    else {
      # theta is given, use usual glm and define family
      fun <- glm
      cl[["family"]] <- eval(parse(text = gettextf("function() negative.binomial(theta=%s)", cl[["theta"]])), parent.frame())
      cl[["theta"]] <- NULL
    }

  } else if(fitfn == "polr"){
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("package 'MASS' must be installed")
    fun <- get("polr", asNamespace("MASS"), inherits = FALSE)
    if(is.null(cl[["Hess"]])) cl[["Hess"]] <- TRUE
    if(is.null(cl[["model"]])) cl[["model"]] <- TRUE

  } else if(fitfn == "lmrob"){
    if (!requireNamespace("robustbase", quietly = TRUE))
      stop("package 'robustbase' must be installed")
    fun <- get("lmrob", asNamespace("robustbase"), inherits = FALSE)

  } else if(fitfn == "tobit"){
    fun <- Tobit

  } else if(fitfn == "zeroinfl"){
    if (!requireNamespace("pscl", quietly = TRUE))
      stop("package 'pscl' must be installed")
    fun <- get("zeroinfl", asNamespace("pscl"), inherits = FALSE)

  } else if(fitfn == "multinom"){
    if (!requireNamespace("nnet", quietly = TRUE))
      stop("package 'nnet' must be installed")
    fun <- get("multinom", asNamespace("nnet"), inherits = FALSE)
    if(is.null(cl[["maxit"]])) cl[["maxit"]] <- 500   # we set a higher default for the iterations here...
    if(is.null(cl[["model"]])) cl[["model"]] <- TRUE
    if(is.null(cl[["trace"]])) cl[["trace"]] <- FALSE

  } else if(fitfn == "randomForest"){
    if (!requireNamespace("randomForest", quietly = TRUE))
      stop("package 'randomForest' must be installed")
    fun <- get("randomForest", asNamespace("randomForest"), inherits = FALSE)
    if(is.null(cl[["na.action"]])) cl[["na.action"]] <- na.omit

  } else if(fitfn == "rpart"){
    if (!requireNamespace("rpart", quietly = TRUE))
      stop("package 'rpart' must be installed")
    fun <- get("rpart", asNamespace("rpart"), inherits = FALSE)
    if(is.null(cl[["model"]])) cl[["model"]] <- TRUE
    if(is.null(cl[["y"]])) cl[["y"]] <- TRUE

  } else if(fitfn == "nnet"){
  if (!requireNamespace("nnet", quietly = TRUE))
    stop("package 'nnet' must be installed")
    fun <- get("nnet", asNamespace("nnet"), inherits = FALSE)
    if(is.null(cl[["maxit"]])) cl[["maxit"]] <- 300   # we set a default for the iterations here...
    if(is.null(cl[["trace"]])) cl[["trace"]] <- FALSE
    if(is.null(cl[["size"]])) cl[["size"]] <- 5

  } else if(fitfn == "C5.0"){
    if (!requireNamespace("C50", quietly = TRUE))
      stop("package 'C50' must be installed")
    fun <- get("C5.0", asNamespace("C50"), inherits = FALSE)

  } else if(fitfn == "lda"){
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("package 'MASS' must be installed")
    fun <- get("lda", asNamespace("MASS"), inherits = FALSE)

  } else if(fitfn == "qda"){
    if (!requireNamespace("MASS", quietly = TRUE))
      stop("package 'MASS' must be installed")
    fun <- get("qda", asNamespace("MASS"), inherits = FALSE)

  } else if(fitfn == "svm"){
    if (!requireNamespace("e1071", quietly = TRUE))
      stop("package 'e1071' must be installed")
    fun <- get("svm", asNamespace("e1071"), inherits = FALSE)
    if(is.null(cl[["probability"]]))
      cl[["probability"]] <- TRUE

  } else if(fitfn == "naive_bayes"){
    if (!requireNamespace("naivebayes", quietly = TRUE))
      stop("package 'naivebayes' must be installed")
    fun <- get("naive_bayes", asNamespace("naivebayes"), inherits = FALSE)

  } else if(fitfn == "lb"){
    fun <- get("LogitBoost", asNamespace("ModTools"), inherits = FALSE)

  } else if(fitfn == "logit"){
    fun <- get("glm", asNamespace("stats"), inherits = FALSE)
    cl[["family"]] <- "binomial"
  }
  cl[["fitfn"]] <- NULL

  # cl[[1]] <- as.name(fun)
  cl[[1]] <- fun

  # evaluate model with given fit function
  res <- eval(cl, parent.frame())
  class(res) <- c("FitMod", class(res))

  if(fitfn == "logit")
    res[["call"]][[1]] <- as.name("glm")

  else if(fitfn %in% c("multinom", "rpart")){
    res[["call"]][[1]] <- as.name(fitfn)

    # extend result object with p-values
    capture.output(sm <- summary(res))
    z <- sm$coefficients / sm$standard.errors

    # 2-tailed Wald z tests to test significance of coefficients
    res[["pval"]] <- (1 - pnorm(abs(z), 0, 1)) * 2

  }


  if(inherits(res, "nnet"))
    # ... and display a warning, if it did not converge
    if(res$convergence == 1)
      warning("nnet() did not (yet) converge, consider increase maxit (and set trace=TRUE)!")

  if(inherits(res, "multinom")){
    res[["drop1"]] <- .drop1.multinom(res)
  }

  if(inherits(res, "polr")){
    ## res[["drop1"]] <- drop1(res, test="Chisq")
    res[["drop1"]] <- .drop1.polr(res)
    if(is.null(cl[["profiling"]]) || !cl[["profiling"]])
      res[["ci"]] <- confint.default(res)
    else
      res[["ci"]] <- confint(res)
  }

  if(inherits(res, "rpart")){
    # add a variables actually used in tree construction
    frame <- res$frame
    leaves <- frame$var == "<leaf>"
    used <- unique(frame$var[!leaves])
    res[["used"]] <- sort(as.character(used))

  }

  if(fitfn %in% c("poisson","gamma", "quasipoisson"))
    res[["call"]][[1]] <- as.name("glm")

  else if(fitfn == "negbin") {
    if (is.null(cl[["family"]]))
      res[["call"]][[1]] <- as.name("glm.nb")
    else
      res[["call"]][[1]] <- as.name("glm")

  } else if(fitfn %in% c("lm","lmrob","polr"))
    res[["call"]][[1]] <- as.name(fitfn)


  res$fitfn <- fitfn

  return(res)

}


predict.FitMod <- function(object, ...){

  if(inherits(object, "multinom")){
    .predict.multinom(object, ...)

  } else if(inherits(object, "nnet")){
    .predict.nnet(object, ...)


  } else if(inherits(object, "glm")){
    .predict.glm(object, ...)

  } else if(inherits(object, "rpart")){
    .predict.rpart(object, ...)

  } else if(inherits(object, "lda")){
    .predict.lda(object, ...)

  } else if(inherits(object, "qda")){
    .predict.qda(object, ...)

  } else if(inherits(object, "svm")){
    .predict.svm(object, ...)

  } else if(inherits(object, "C5.0")){
    .predict.C5.0(object, ...)

  } else {
    NextMethod(object, ...)
  }
}


summary.FitMod <- function(object, ...){
  NextMethod(object, ...)
}


drop1.FitMod <- function(object, ...){
  if(inherits(object, "multinom"))
    .drop1.multinom(object, ...)
  else
    NextMethod(object, ...)
}



print.FitMod <- function(x, ...){
  if(inherits(x, "multinom"))
    .print.multinom(x, ...)

  # else if(inherits(x, "nnet"))
  #   .print.nnet(x, ...)

  else if(inherits(x, "polr"))
     .print.polr(x, ...)

  else
    NextMethod(x, ...)
}

rpart:::plot.rpart


plot.FitMod <- function(x, ...){

  if(inherits(x, "rpart")){
    plot.rpart(x, ...)

    # evt. other defaults:
    # plot.rpart <- function (x = stop("no 'x' arg"), type = 0, extra = 0, under = FALSE,
    #           clip.right.labs = TRUE, fallen.leaves = FALSE, branch = if (fallen.leaves) 1 else 0.2,
    #           uniform = TRUE, digits = 2, varlen = -8, faclen = 3, cex = NULL,
    #           tweak = 1, compress = TRUE, ycompress = uniform, snip = FALSE,
    #           ...){


  } else if(inherits(x, "nnet")){
    # print.nnet <- nnet:::print.nnet
    # print.nnet <- get("print.nnet", asNamespace("NeuralNetTools"),
    #                      inherits = FALSE)
    plotnet(x, ...)

  } else {
    NextMethod(x, ...)
  }
}

as.FitMod <- function(x){
  structure(x, class=c("FitMod", class(x)))
}


update.FitMod <- function(object, ...){
  # update kills the model class, so restore here
  oclass <- class(object)
  class(object) <- class(object)[class(object) != "FitMod"]
  res <- update(object, ...)
  class(res) <- oclass
  return (res)
}




Response <- function(x, ...){

  # response name
  # deparse(attr(terms(model), "variables")[[2]])


  if(inherits(x, "C5.0")){
    x$terms <- eval(x$call$formula)
    res <- model.response(model.frame(x))

  } else if(inherits(x, "rpart") & is.null(x$model)){
    res <- factor(attr(x, "ylevels")[x$y])

  } else if(inherits(x, "naive_bayes")){
    x$terms <- eval(x$call$formula)
    res <- model.response(model.frame(x))

  } else {
    res <- model.response(model.frame(x))
  }

  # get name of response
  attr(res, "response") <- attr(attr(x$terms,"dataClasses"),"names")[attr(x$terms,"response")]

  return(res)

}


ResponseName <- function(x) {

  # https://stat.ethz.ch/pipermail/r-help/2007-November/145649.html
  f <- formula(x)
  if(length(f) < 3) stop("no response")

  resp <- f[[2]]

  if(!is.name(resp)) stop("response is not a name")

  as.character(resp)
}



Conf.FitMod <- function(x, ...){
  if(inherits(x, c("C5.0", "nnet", "naive_bayes", "lda", "LogitBoost")))
    Conf(x=predict(x, type="class"), ref=Response(x), ... )
  else
    NextMethod(x, ...)
}



Cstat.FitMod <- function(x, pos=NULL, ...) {
  pred <- predict(x, type = "prob")
  if(is.null(pos))
    pos <- 2
  if(ncol(pred)==2)
    Cstat(pred[, pos], Response(x), ...)
  else
    NA
}


# test  -------------------------


BreuschPaganTest <- function(formula, varformula = NULL, studentize = TRUE,
                             data = list()) {

  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
    Z <- if(is.null(varformula)) X
    else model.matrix(varformula, data = data)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if(is.null(varformula)) X
    else model.matrix(varformula, data = data)
  }

  ## only use complete cases that are in both models
  if(!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% row.names(X))))) {
    allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
    X <- X[allnames,]
    Z <- Z[allnames,]
    y <- y[allnames]
  }

  k <- ncol(X)
  n <- nrow(X)

  resi <- lm.fit(X,y)$residuals
  sigma2 <- sum(resi^2)/n

  if(studentize)
  {
    w <- resi^2 - sigma2
    fv <- lm.fit(Z,w)$fitted
    bp <- n * sum(fv^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
  }
  else
  {
    f <- resi^2/sigma2 -1
    fv <- lm.fit(Z,f)$fitted
    bp <- 0.5 * sum(fv^2)
    method <- "Breusch-Pagan test"
  }

  names(bp) <- "BP"
  df <- ncol(Z)-1
  names(df) <- "df";
  RVAL <- list(statistic = bp,
               parameter = df,
               method = method,
               p.value= pchisq(bp,df,lower.tail=FALSE),
               data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}




TModC <- function(..., newdata=NULL, reference=NULL, ord=NULL){

  graphrank <- function(i)
    paste(strrep(".", i-1), "x", strrep(".", max(i) - i), sep="")


  mods <- list(...)
  if(length(mods)==1)
    mods <- unlist(mods, recursive = FALSE)

  if(!is.null(newdata)){
    rocs <- list()
    for(i in seq_along(mods))
      rocs[[i]] <-  ROC(predict(mods[[i]], newdata=newdata, type="prob")[, 2],
                        resp = reference, plotit=FALSE)

  } else {
    rocs <- lapply(mods, ROC)
  }


  xname <- names(mods)
  if(is.null(xname))
    xname <- paste("model", seq_along(mods))


  auc <- sapply(rocs, function(x) x$auc)
  bs <- sapply(mods, BrierScore)

  res <- data.frame(auc=sapply(rocs, function(x) x$auc),
                    auc_p=auc/max(auc),
                    auc_rnk=rank(-auc),
                    bs=sapply(mods, BrierScore),
                    bs_p=((1-bs)/max(1-bs)),
                    bs_rnk=rank(bs))
  res$auc_grnk <- graphrank(res$auc_rnk)
  res$bs_grnk <- graphrank(res$bs_rnk)


  ord <- match.arg(ord, choices = c("ensemble", colnames(res)))

  if(ord=="bs")
    res <- Sort(res, ord="bs", decreasing=FALSE)

  else if(ord=="ensemble"){
    id <- order(apply(res[, c("auc_p","bs_p")], 1, mean),
                decreasing = TRUE)
    res <- res[id, ]

  } else {
    res <- Sort(res, ord=ord, decreasing=TRUE)

  }

  res <- list(tab=res, rocs=rocs, xname=xname)

  class(res) <- c("TModC", class(res))
  res
}



print.TModC <- function(x, ...){

  print.data.frame(Format(x$tab, digits=c(3,3,0,3,3,0,0,0)),
                   print.gap=2, ...)

}


plot.TModC <- function(x, col=NULL, args.legend=NULL, ...){

  if(is.null(col))
    cols <- Pal("Tibco")
  else
    cols <- col


  sapply(seq_along(x$rocs),
         function(i) plot(x$rocs[[i]], add=(i!=1), col=cols[i]))

  grid()

  args.legend1 <- list(x = "bottomright", inset = 0, legend = x$xname,
                       fill = cols, bg = "white", cex = 1.0)
  if (!is.null(args.legend)) {
    args.legend1[names(args.legend)] <- args.legend
  }
  add.legend <- TRUE
  if (!is.null(args.legend))
    if (identical(args.legend, NA)) {
      add.legend <- FALSE
    }
  if (add.legend)
    DoCall("legend", args.legend1)


}




Conf.polr <- function(x, ...) {

  # extract response from the model

  Conf(x=predict(x),
       ref=model.extract(model.frame(x), "response") , ... )
}

.drop1.multinom <-  function (object, scope, test = c("Chisq","none"), ...) {

  if (!inherits(object, "multinom"))
    stop("Not a multinom fit")
  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, attr(object$terms, "term.labels"),
                   nomatch = FALSE)))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>",
                                                           scope), c("Df", "AIC")))
  ans[1, ] <- c(object$edf, object$AIC)
  if (test[1]=="Chisq") ans <- cbind(ans,Chisq=NA, p.value=NA)
  env <- environment(formula(object))
  for(i in seq(ns)) {
    tt <- scope[i]

    nfit <- update(object, as.formula(paste("~ . -", tt)), evaluate = FALSE)

    # this is additionally to automatically handle missing values
    # it is not always clear, what impact missings have, so handle with care
    # maybe place a warning... see help drop1.multinom
    # !!!!
    nfit$data <- str2lang(gettextf("na.omit(%s[, c(%s)])",
                        deparse(nfit$data), paste(sQuote(all.vars(formula(object))), collapse=",")))

    nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)

    if (nfit$edf == object$edf)
      nfit$AIC <- NA
    ans[i+1, ] <- c(nfit$edf, nfit$AIC,
                    if (test[1]=="Chisq") unlist(anova(object, nfit)[2,6:7]))
  }
  as.data.frame(ans)
}


OddsRatio.polr <- function(x, conf.level=NULL, digits=3, ...) {

  if(is.null(conf.level)) conf.level <- 0.95

  # class(x) <- class(x)[class(x)!="regr"]
  r.summary <- summary(x, Wald.ratios = TRUE)

  coe <- SetNames(r.summary$coefficients[,1, drop=FALSE], colnames="or")
  se <- r.summary$coefficients[,2]

  d.res <- r.summary

  d.print <- data.frame(
    "or"= Format(exp(coe[, "or"]), digits=digits),
    "or.lci" = Format(exp(coe[, "or"] + qnorm(0.025) * se), digits=digits),
    "or.uci" = Format(exp(coe[, "or"] - qnorm(0.025) * se), digits=digits),
    "pval" = Format(2*(1-pnorm(q = abs(coe[, "or"]/se), mean=0, sd=1)), fmt="p", digits=3),
    "sig" = Format(2*(1-pnorm(q = abs(coe[, "or"]/se), mean=0, sd=1)), fmt="*"),
    stringsAsFactors = FALSE
  )

  colnames(d.print)[4:5] <- c("Pr(>|z|)","")
  rownames(d.print) <- paste(coe$time, coe$id, sep=":")

  res <- list(or = d.print, call = x$call,
              BrierScore = NA, # BrierScore(x),
              PseudoR2 = PseudoR2(x, which="all"), res=d.res)

  class(res) <- "OddsRatio"
  return(res)

}


.drop1.polr <- function (mod, ...) {

  # based on car:::Anova.II.polr

  term.names <-
    function (model, ...)
    {
      term.names <- labels(terms(model))
      if (any(names(coefficients(model)) == "(Intercept)"))  # (has.intercept(model))
        c("(Intercept)", term.names)
      else term.names
    }


  df.terms.polr <-  function (model, term, ...) {

    if (!missing(term) && 1 == length(term)) {
      assign <- attr(model.matrix(model), "assign")
      which.term <- which(term == labels(terms(model)))
      if (0 == length(which.term))
        stop(paste(term, "is not in the model."))
      sum(assign == which.term)
    }
    else {
      terms <- if (missing(term))
        labels(terms(model))
      else term
      result <- numeric(0)
      for (term in terms) result <- c(result, Recall(model,
                                                     term))
      names(result) <- terms
      result
    }
  }

  relatives <- function (term, names, factors) {
    is.relative <- function(term1, term2) {
      all(!(factors[, term1] & (!factors[, term2])))
    }
    if (length(names) == 1)
      return(NULL)
    which.term <- which(term == names)
    (1:length(names))[-which.term][sapply(names[-which.term],
                                          function(term2) is.relative(term, term2))]
  }

  responseName <- function(model, ...)
    deparse(attr(terms(model), "variables")[[2]])


  if (!requireNamespace("MASS"))
    stop("MASS package is missing")
  which.nms <- function(name) which(asgn == which(names == name))

  fac <- attr(terms(mod), "factors")
  names <- term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  wt <- model.weights(model.frame(mod))
  asgn <- attr(X, "assign")
  aic <- p <- LR <- rep(0, n.terms)
  df <- df.terms.polr(mod)

  for (term in 1:n.terms) {
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
                                         which.nms)))
    mod.1 <- if (n.terms > 1)
      MASS::polr(y ~ X[, -c(1, exclude.1)], weights = wt)
    else MASS::polr(y ~ 1, weights = wt)

    dev.1 <- deviance(mod.1)
    mod.2 <- if (length(rels) == 0)
      mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      MASS::polr(y ~ X[, -c(1, exclude.2)], weights = wt)
    }
    dev.2 <- deviance(mod.2)
    LR[term] <- dev.1 - dev.2
    p[term] <- pchisq(LR[term], df[term], lower.tail = FALSE)
    aic[term] <- AIC(mod.1)
  }

  result <- data.frame(aic, LR, df, p)
  row.names(result) <- names
  names(result) <- c("AIC", "LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
                               paste("Response:", responseName(mod)))
  result

}



LowVar <- function(x, na.rm=FALSE, uvalp = 0.2, freq = 20, method=c("both", "unique", "freq")) {

  # The functionality is realized in two main steps:
  # 1. Check for near zero variance predictors and flag as near zero if:
  #   (a) the percentage of unique values is less than 20
  #   (b) the ratio of the most frequent to the second most frequent value is greater than 20,
  # 2. Check for susceptibility to multicollinearity
  #   (a) Calculate correlation matrix
  #   (b) Find variables with correlation 0.9 or more and delete them



  if (na.rm) x <- na.omit(x)
  method <- match.arg(method)

  lw1 <- (1 - sum(duplicated(x)) / length(x)) < uvalp

  m1 <- DescTools::Mode(x, na.rm)
  if(identical(m1, structure(NA_real_, freq = 1L)))
    # all values are unique, lw2 = FALSE
    lw2 <- FALSE

  else {
    m2 <- DescTools::Mode(x[x != m1], na.rm)
    lw2 <- (attr(m1, "freq") / attr(m2, "freq")) > freq

  }

  if(method == "unique")
    res <- lw1

  else if(method == "freq")
    res <- lw2

  else if(method == "both")
    res <- lw1 | lw2

  return(res)

}



RefLevel <- function(x){

  refCat <- function(model, var) {
    cs <- attr(model.matrix(model), "contrasts")[[var]]
    if (is.character(cs)) {
      if (cs == "contr.treatment")
        ref <- 1
      else stop("No treatment contrast")
    }
    else {
      zeroes <- !cs
      ones <- cs == 1
      stopifnot(all(zeroes | ones))
      cos <- colSums(ones)
      stopifnot(all(cos == 1))
      ros <- rowSums(ones)
      stopifnot(sum(!ros) == 1 && sum(ros) != ncol(cs))
      ref <- which(!ros)
    }
    return(levels(model$model[[var]])[ref])
  }

  # find factor predvals
  fpred <- names(grep("factor", attr(x[["terms"]], "dataClasses"), value=TRUE))
  resp <- all.vars(formula(x))[1]
  fpred <- fpred[fpred != resp]

  # list( complicated=
  # sapply(fpred , function(z) refCat(x, z))
  #
  # # or simply?
  # , simple=sapply(x$xlevels, "[", 1))

  sapply(fpred , function(z) refCat(x, z))

}


.InsertRefLevel <- function(ref, coeftable, pval){

  # call: .InsertRefLevel(ref=RefLevel(r.lm)$simple,
  #                       coefficients(r.lm))

  # insert reference levels for factors in a coefficient table
  # rr <- RefLevel(x)$simple

  for(i in seq_along(ref)){
    ff <- names(ref[i])
    rl <- ref[i]
    j <- head(grep(ff, rownames(coeftable)), 1)
    coeftable <- Append(coeftable, NA, after=j-1, rows=TRUE)
    rownames(coeftable)[grep(ff, rownames(coeftable))] <-
      gsub(ff, "  ", rownames(coeftable)[grep(ff, rownames(coeftable))])
    rownames(coeftable)[j] <- gettextf("%s (Ref: %s)", ff, rl)

    if(!is.null(pval))
      coeftable[j, 8] <- Format(pval[grep(ff, rownames(pval)), ], fmt="*")

  }

  as.data.frame(coeftable)

}



InsertRow <- function(x, val, after=nrow(x)) {

  # insert a row in a data.frame

  x[seq(after+1, nrow(x)+1), ] <- x[seq(after, nrow(x)), ]
  x[after, ] <- val

  x

}


OverSample <- function(x, vname){
  p <- sort(prop.table(table(x[, vname])))[1]
  xbig <- x[x[,vname] != names(p),]
  return(rbind(xbig, Sample(x[x[,vname] == names(p),], size=nrow(xbig), replace=TRUE)))
}


UnderSample <- function(x, vname){
  p <- sort(prop.table(table(x[, vname])))[1]
  xsmall <- x[x[,vname] == names(p),]
  return(rbind(xsmall, Sample(x[x[,vname] != names(p),], size=nrow(xsmall), replace=FALSE)))
}




RobSummary <- function(mod, conf.level=0.95, type="HC0"){

  sterr <- sqrt(diag(vcovHC(mod, type=type)))
  alpha <- 1 - (1-conf.level)/2

  if(inherits(mod, "glm")) {
    res <- cbind(est= coef(mod),
                 lci = coef(mod) - qnorm(alpha) * sterr,
                 uci = coef(mod) + qnorm(alpha) * sterr,
                 rse = sterr,
                 zval = coef(mod)/sterr,
                 pval = 2 * pnorm(abs(coef(mod)/sterr), lower.tail=FALSE))

  } else if(inherits(mod, "lm")) {

    res <- cbind(est= coef(mod),
                 lci = coef(mod) - qt(alpha, df = mod$df) * sterr,
                 uci = coef(mod) + qt(alpha, df = mod$df) * sterr,
                 rse = sterr,
                 zval = coef(mod)/sterr,
                 pval = 2 * pt(abs(coef(mod)/sterr), df=mod$df, lower.tail=FALSE)
    )

  }

  return(res)

}




PredictCI <- function(mod, newdata, conf.level=0.95){

  preds <- predict(mod, newdata=newdata, type="link", se.fit=TRUE)
  critval <- -qnorm((1-conf.level)/2)
  res <- cbind(
    fit <- preds$fit
    , lwr <- preds$fit - critval * preds$se.fit
    , upr <- preds$fit + critval * preds$se.fit
  )

  res <- SetNames(mod$family$linkinv(res),
                  colnames=c("prob", "lci", "uci"))

  return(res)
}


























