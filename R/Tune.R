
Tune <- function(x, scale = FALSE, sort=TRUE, type=NULL, ...){
  UseMethod("Tune")
}





# Tune.FitMod <- function(x, scale = FALSE, sort=TRUE, type=NULL, ...){
#
#   if(inherits(x, "lm")){
#     VarImp.lm(x, scale=scale, sort=sort, type=type, ...)
#
#   } else if(inherits(x, "nnet")){
#     .predict.nnet(x, ...)
#
#
#   } else if(inherits(x, "glm")){
#     .VarImp.glm(x, scale=scale, sort=sort, type=type, ...)
#
#   } else if(inherits(x, "rpart")){
#     .VarImp.rpart(x, scale=scale, sort=sort, type=type, ...)
#
#   } else if(inherits(x, "rf")){
#     tuneRF(x, scale=scale, sort=sort, type=type, ...)
#
#   } else if(inherits(x, "nnet")){
#     .VarImp.nnet(x, scale=scale, sort=sort, type=type, ...)
#
#   } else if(inherits(x, "C5.0")){
#     .VarImp.C5.0(x, scale=scale, sort=sort, type=type, ...)
#
#   } else {
#     NextMethod(x, sort=sort, type=type, ...)
#   }
#
# }


Tune.nnet <- function(x, ..., testset=NULL){

  mpar <- list(...)

  tpar <- do.call(expand.grid, mpar)

  lst <- list()
  for(i in 1:nrow(tpar)){
    for(j in 1:length(mpar))
      x$call[names(mpar[j])] <- tpar[i, names(mpar[j])]
    lst[[i]] <- eval(x$call)
    class(lst[[i]]) <- c("FitMod", class(lst[[i]]))
  }

  tpar$acc <- sapply(lst, function(x) Conf(x)$acc)
  if(!is.null(testset)){
    xresp <- attr(Response(x), "response")
    tpar$test_acc <- sapply(lst, function(x) Conf(predict(x, newdata=testset),
                                                  ref=testset[, xresp])$acc)

  }
  return(list(modpar=Sort(tpar, ncol(tpar), decreasing = TRUE), mods=lst))

}



Tune.default <- function(x, ...){

  cat(gettextf("No tune procedure found for %s.\n", class(x)))
}





print.Tune <- function(x, digits = 2, ...) {
  # print as data.frame if something was changed
  print(round(data.frame(x[]), digits=digits))
}




