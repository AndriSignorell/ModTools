

# Original:
# library(BCA)
# lift.chart


## Lift charts and scoring
PlotLift <- function(modelList, data, targLevel, trueResp, type="cumulative",
                       sub="") {
  if(type != "cumulative" & type != "incremental") {
    stop("An improper lift chart type is specified.")
  }
  # not needed at all here... 17.5.2021
  # set.seed(1)
  data <- data[order(runif(nrow(data))), ]
  yvar1 <- rep(NA, length(modelList))
  modAvail <- rep(NA, length(modelList))
  probVar <- NULL
  for(i in 1:length(modelList)) {
    mod <- eval(parse(text=modelList[i]))
    modtype <- class(mod)[1]
    if(modtype != "FitMod" & modtype != "glm" & modtype != "rpart" & modtype != "nnet.formula") {
      stop("Models can only be estimated using glm, rpart, or nnet.")
    }
    yvar1[i] <- as.character(mod$call$formula)[2]
    xvars <- unlist(strsplit(as.character(mod$call$formula)[3]," + ",
                             fixed=TRUE))
    if(!all(xvars %in% names(data))) {
      probVar <- c(probVar, xvars[!(xvars %in% names(data))])
      modAvail[i] <- FALSE
    }
    else {
      modAvail[i] <- TRUE
    }
  }
  if(any(yvar1!=yvar1[1])) {
    stop("Not all the models have the same dependent variable")
  }
  yvar2 <- data[[yvar1[1]]]
  if(!is.factor(yvar2)) {
    stop("The y variable must be a two-level factor.")
  }
  if(length(levels(yvar2)) != 2) {
    stop("The y variable must be a two-level factor.")
  }
  if(!any(as.character(yvar2) == targLevel)) {
    stop(paste('None of the levels of the response variable is "', targLevel,
               '".', sep=""))
  }
  yvar <- as.numeric(yvar2 == targLevel)
  sampResp <- sum(yvar)/length(yvar)

  # print(sampResp)

  if(length(probVar)>0) {
    probVar <- unique(probVar)
    probModel <- modelList[!modAvail]
    warnString <- paste("The models", paste(probModel, collapse=", "),
                        "are not in the lift chart because the variables",
                        paste(probVar, collapse=", "), "are not available.")
    #warning(warnString, call.=FALSE, immediate.=TRUE)
    warning(warnString)
    modelList <- modelList[modAvail]
    if(length(modelList)==0) {
      # Message(message=gettextRcmdr(paste(
      #   "All models are missing at least one of the variables: ",
      #   paste(probVar, collapse=", "), ".", sep="")), type="error")
      return()
    }
  }
  sampWt <- (sampResp*(1 - trueResp))/(trueResp*(1 - sampResp))
  nmodels <- length(modelList) # new
  colr <- rep(palette(), ceiling(nmodels/length(palette()))) # new
  if(type == "cumulative") {
    plot(seq(0.1, 1, 0.1), seq(0.1, 1, 0.1), main=
           "Weighted Cumulative Response Captured", sub=sub, xlab="Sample Proportion",
         ylab="Percent of Total Response Captured", type="l", lwd=2, xaxs="i",
         yaxs="i")
    for(i in 1:nmodels) {
      model1 <- eval(parse(text=modelList[i]))
      modtype <- class(model1)[1]
      model <- eval(model1$call)
      if(inherits(model1, c("nnet.formula", "glm"))) {
        if(levels(yvar2)[1] == targLevel) {
          var1 <- yvar[order(predict(model, newdata = data),
                             decreasing = FALSE)]
        }
        else {
          var1 <- yvar[order(predict(model, newdata = data),
                             decreasing = TRUE)]
        }
      }
      else {
        var1 <- yvar[order(as.vector(predict(
          model, newdata = data, method="class")[ , targLevel]), decreasing = TRUE)]
      }
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- sampWt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1),
                     include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      lines(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[i], lwd = 2)
      points(seq(0.1, 1, 0.1), cumsum(var2)/sum(var2), col = colr[i], pch=i)
    }
    legend("bottomright", legend=modelList, col=colr[1:nmodels],
           pch=1:length(modelList), lty=1, lwd=2)
  }
  else {
    resp.matrix <- matrix(NA, nrow=10, ncol=length(modelList))
    for(i in 1:nmodels) {
      model1 <- eval(parse(text=modelList[i]))
      modtype <- class(model1)[1]
      model <- eval(model1$call)
      if(inherits(model1, c("nnet.formula", "glm"))) {
        if(levels(yvar2)[1] == targLevel) {
          var1 <- yvar[order(predict(model, newdata = data),
                             decreasing = FALSE)]
        }
        else {
          var1 <- yvar[order(predict(model, newdata = data),
                             decreasing = TRUE)]
        }
      }
      else {
        var1 <- yvar[order(as.vector(predict(
          model, newdata = data)[ , targLevel]), decreasing = TRUE)]
      }
      var.ind1 <- rep(1, length(var1))
      var.ind1[var1 == 0] <- sampWt
      var.ind <- cut(cumsum(var.ind1)/sum(var.ind1), seq(0, 1, 0.1),
                     include.lowest = TRUE)
      var2 <- as.vector(by(var1, var.ind, sum))
      var3 <- as.vector(by(var.ind1, var.ind, sum))
      resp.matrix[, i] <- var2/var3
    }
    max.resp <- max(resp.matrix)
    plot(seq(0.1, 1, 0.1), seq(0, max.resp, length=10), type="n", main=
           "Weighted Incremental Response Rate", sub=sub, xlab="Sample Percentile",
         ylab="Resposne Rate", lwd=2, xaxs="i", yaxs="i")
    lines(seq(0.1, 1, 0.1), rep(trueResp, 10), lwd=2)
    for(j in 1:nmodels) {
      lines(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], lwd=2)
      points(seq(0.1, 1, 0.1), as.vector(resp.matrix[, j]), col=colr[j], pch=j)
    }
    print(length(modelList))
    legend("topright",legend=modelList,col=colr[1:nmodels],
           pch=1:nmodels, lty=1, lwd=2)
  }

  invisible(sampResp)

}

