

# Variable importance for lm, glm, nnet, C5.0, rpart, randomForest, svm, logitBoost, naiveBayes  -----------

VarImp <- function(x, scale = FALSE, sort=TRUE, ...){
  UseMethod("VarImp")
}


plot.VarImp <- function(x, sort=TRUE, maxrows=NULL, cex=1.1,
                        main="Variable importance", ...){
  dotchart(rev(x), labels=rev(rownames(x)), cex=cex, main=main, ...)
}




VarImp.FitMod <- function(x, scale = FALSE, sort=TRUE, type=NULL, ...){

  if(inherits(x, "glm")){
    res <- .VarImp.glm(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "rpart")){
    res <- .VarImp.rpart(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "randomForest")){
    res <- .VarImp.randomForest(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "nnet")){
    res <- .VarImp.nnet(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "C5.0")){
    res <- .VarImp.C5.0(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "svm")){
    res <- .VarImp.svm(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "lm")){
     res <- VarImp.lm(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "nnet")){
     res <- .VarImp.gbm(x, ...)

  } else {
    res <- NextMethod(x, scale=scale, sort=sort, type=type, ...)
  }


  res <- as.matrix(res)
  class(res) <- c("VarImp", class(res))
  return(res)

}





print.VarImp <- function(x, digits = 3, ...) {

  # print as data.frame if something was changed
  # print(round(data.frame(x[]), digits=digits))
  if(!identical(x, NA))
    print(Format(data.frame(x), digits=digits))
  else
    print(NA_real_)
}



VarImp.lm <- function(x, scale = FALSE, sort=TRUE, type = NULL, ...)  {

#   example:
#   r.lmx <- regr(Fertility ~ ., swiss)
#   VarImp(r.lmx, type="level")
#   VarImp(r.lmx, type=c("lmg","level"))
#
#   VarImp(r.lmx, type="lmg")
#   VarImp(r.lmx, type=c("lmg","pratt"))



  if(is.null(type)) type <- "lmg"
  # variable importances following library(relaimpo), calc.relimp
  crtype <- type[type %in% c("lmg", "pmvd", "first", "last", "betasq", "pratt")]
  if(length(crtype) > 0){
    res <- calc.relimp(x, type = crtype)
    res <- sapply(crtype, function(x) slot(res,x))
  }

  # implementation uses $allcoef from regr! so lm would not work.... ********
  if(any(type == "level")) {
    # get the level importance, Achen 1982

    ## *********** CHECKME ********* CHECKME *********** CHECKME ************
    # transformations are not handled correct here!! Correct!!!!!!!!!!!!!

    varnames <- names(x$allcoef)[ names(x$allcoef) == names(x$allvars)]
    vx <- numeric(length(varnames))
    levimp <- numeric(length(varnames))
    for(i in 1:length(varnames)) {
      xx <- x$allvars[varnames[i]][,]
      if(is.factor(xx))  {
        vx[i] <- sum(prop.table(table(xx)) * unlist(x$allcoef[varnames[i]]))
      } else {
        vx[i] <- mean(xx, na.rm=TRUE) * unlist(x$allcoef[varnames[i]])
      }
    }
    # percentage of mean(response)
    res.lvl <- matrix(vx / mean(x$response[,], na.rm=TRUE), dimnames=list(varnames,"level"))

    # bind level importance with potential relimpo results
    if(exists("res", inherits=FALSE))  res <- cbind(res, "level"=res.lvl[rownames(res),])
    else res <- res.lvl
  }

  if(sort)
    res <- Sort(res, decreasing=TRUE)

  return(res)

}




.VarImp.rpart <- function(x, scale = TRUE, sort=TRUE, ...){

  res <- as.matrix(x$variable.importance)
  if(scale)
    res <- res/sum(res)
  colnames(res) <- "varimp"
  if(sort) res <- Sort(res, decreasing=TRUE)

  res

}

.VarImp.randomForest <- function(x, scale = TRUE, sort = TRUE, ...){
  res <- as.matrix(importance(x, ...))
  colnames(res) <- "varimp"
  if(sort) res <- Sort(res, decreasing=TRUE)

}


.VarImp.C5.0 <- function(x, scale = TRUE, sort = TRUE, ...) {
  C5imp(object=x, ...)
}


GarsonWeights <- function(object) {

  beta <- coef(object)
  abeta <- abs(beta)
  nms <- names(beta)
  i2h <- array(NA, dim = object$n[2:1])
  h2o <- array(NA, dim = object$n[2:3])

  for (hidden in 1:object$n[2]) {
    for (input in 1:object$n[1]) {
      label <- paste("i", input, "->h", hidden,"$", sep = "")
      i2h[hidden, input] <- abeta[grep(label, nms, fixed = FALSE)]
    }
  }
  for(hidden in 1:object$n[2]){
    for(output in 1:object$n[3]){
      label <- paste("h", hidden, "->o",
                     ifelse(object$n[3] == 1, "", output),
                     sep = "")
      h2o[hidden,output] <- abeta[grep(label, nms, fixed = TRUE)]
    }
  }

  if(FALSE)
  {
    ## Test case from Gevrey, M., Dimopoulos, I., & Lek,
    ## S. (2003). Review and comparison of methods to study the
    ## contribution of variables in artificial neural network
    ## models. ecological modelling, 160(3), 249-264.
    i2h <- matrix(c(-1.67624,  3.29022,  1.32466,
                    -0.51874, -0.22921, -0.25526,
                    -4.01764,  2.12486, -0.08168,
                    -1.75691, -1.44702,  0.58286),
                  ncol = 3, byrow = TRUE)
    h2o <- matrix(c(4.57857, -0.48815, -5.73901, -2.65221),
                  ncol = 1)
  }

  ##  From Gevrey et al. (2003): "For each hidden neuron i, multiply
  ##  the absolute value of the hidden-output layer connection
  ##  weight by the absolute value of the hidden-input layer
  ##  connection weight. Do this for each input variable j. The
  ##  following products Pij are obtained"


  ## We'll do this one response at a time. Gevrey et al. (2003) do
  ## not discuss multiple outputs, but the results are the same (at
  ## least in the case of classification).

  imp <- matrix(NA, nrow = object$n[1], ncol = object$n[3])


  for(output in 1:object$n[3])
  {
    Pij <- i2h * NA
    for(hidden in 1:object$n[2]) Pij[hidden,] <- i2h[hidden,] * h2o[hidden,output]

    ## "For each hidden neuron, divide Pij by the sum for all the
    ## input variables to obtain Qij. For example for Hidden 1, Q11 =
    ## P11/(P11+P12+P13).

    Qij <- Pij * NA
    for(hidden in 1:object$n[2]) Qij[hidden,] <- Pij[hidden,] / sum(Pij[hidden,])


    ## "For each input neuron, sum the product Sj formed from the
    ## previous computations of Qij. For example, S1 =
    ## Q11+Q21+Q31+Q41."

    Sj <- apply(Qij, 2, sum)

    ## "Divide Sj by the sum for all the input variables. Expressed as
    ## a percentage, this gives the relative importance or
    ## distribution of all output weights attributable to the given
    ## input variable. For example, for the input neuron 1, the
    ## relative importance is equal to (S1/100)/(S1+S2+S3)"

    imp[,output] <- Sj/sum(Sj)*100
    rm(Pij, Qij, Sj)
  }

  colnames(imp) <- if(!is.null(colnames(object$residuals))) colnames(object$residuals) else paste("Y", 1:object$n[3], sep = "")
  rownames(imp) <- if(!is.null(object$coefnames)) object$coefnames else  paste("X", 1:object$n[1], sep = "")

  class(imp) <- c("VarImp", class(imp))
  imp
}




.VarImp.nnet <- function(x, scale = FALSE, sort=TRUE, ...) {

  object <- x

  imp <- GarsonWeights(object)
  if(ncol(imp) > 1) {
    imp <- cbind(apply(imp, 1, mean), imp)
    colnames(imp)[1] <- "Overall"

  } else {
    imp <- Sort(as.data.frame(imp))
    colnames(imp) <- "Overall"
  }

  if(!is.null(object$xNames))
    rownames(imp) <- object$xNames

  # structure(imp, class="VarImp")
  if(sort)
    Sort(imp, decreasing=TRUE)
  else
    imp

}


# VarImp.gam <- function(x, scale = FALSE, ...){
#
#   object <- x
#
#   if(any(names(object) %in% c("edf", "mgcv.conv", "gcv.ubre")))
#   {
#     library(mgcv)
#     tmp <- mgcv:::anova.gam(object)
#     smoothed <- data.frame(Overall = tmp$s.table[, 4])
#
#     if(nrow(tmp$p.table) > 1)
#     {
#       linear <- data.frame(Overall = tmp$p.table[, 4])
#       out <- rbind(linear, smoothed)
#     } else out <- smoothed
#     out$Overall[!is.na(out$Overall)] <- -log10(out$Overall[!is.na(out$Overall)])
#
#     out$Overall[is.na(out$Overall)] <- 0
#
#     nms <- strsplit(rownames(out), "[()]")
#     nms <- unlist(lapply(nms, function(x) x[length(x)]))
#     rownames(out) <- nms
#     out <- subset(out, rownames(out) != "Intercept")
#   } else {
#     library(gam)
#     trms <- attr(object$terms, "term.labels")
#
#     vars <- all.vars(object$terms)[-1]
#     out <- data.frame(Overall = rep(NA, length(vars)))
#     rownames(out) <- vars
#     for(i in seq(along = trms))
#     {
#       reduced <- update(object, as.formula(paste("~.-", trms[i])))
#       out[i,1] <- -log10(gam:::anova.gam(object, reduced)[2, "P(>|Chi|)"])
#     }
#   }
#
#   structure(out, class="VarImp")
#   # out
#
# }


# still to do:
.VarImp.gbm <- function(x, scale = FALSE, sort = TRUE, ...){  }

.VarImp.svm <- function(x, scale = FALSE, sort = TRUE, ...){
  return(NA)
}



.VarImp.glm <- function(x, scale = FALSE, sort = TRUE, ...){

  # this is simply the z-value of the model,
  # z-values are uniformly distributed random values
  # therefore it's a completely useless definition of variable importance

  # Look for alternatives!!!  *******************

  object <- x

  sObj <- summary(object)$coef
  sObj <- sObj[rownames(sObj) != "(Intercept)",]
  sObj <- abs(sObj[, "z value", drop = FALSE])
  colnames(sObj) <- "Overall"
  sObj <- as.data.frame(sObj)
  sObj
}



.VarImp.multinom <- function(x, scale = FALSE, sort = TRUE, ...) {

  out <- abs(coef(x))
  if(is.vector(out))
  {
    out <- data.frame(Overall = out)
    rownames(out) <- names(coef(x))
  } else {
    out <- as.data.frame(apply(out, 2, sum))
    names(out)[1] <- "Overall"
  }
  subset(out, rownames(out) != "(Intercept)")

}




VarImp.default <- function(x, scale = FALSE, sort = TRUE, ...){
  warning(gettextf("no VarImp definition found for %s (class %s)", deparse(substitute(x)),
               paste(class(x), collapse = ", ")))
  return(NA)
}





