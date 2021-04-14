
Tune <- function(x, ..., testset=NULL, keepmod=TRUE){
  UseMethod("Tune")
}



Tune.default <- function(x, ..., testset=NULL, keepmod=TRUE){

  mpar <- list(...)

  tpar <- do.call(expand.grid, mpar)
  tpar[sapply(mpar, is.character)] <- lapply(tpar[sapply(mpar, is.character)], as.character)

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

  if(keepmod)
    res <- list(modpar=tpar, mods=lst)
  else
    res <- list(modpar=tpar, mods=NULL)

  class(res) <- "Tune"
  return(res)


}



print.Tune <- function(x, digits = 4, ...) {
  # print as data.frame if something was changed
  xp <- data.frame(x$modpar)
  xp[sapply(xp, is.numeric)] <- Format(xp[sapply(xp, is.numeric)], digits=digits)
  xp <- Sort(xp, ncol(xp), decreasing = TRUE)
  print(xp)
}




