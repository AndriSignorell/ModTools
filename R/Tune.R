
Tune <- function(x, scale = FALSE, sort=TRUE, type=NULL, ...){
  UseMethod("Tune")
}





Tune.FitMod <- function(x, scale = FALSE, sort=TRUE, type=NULL, ...){

  if(inherits(x, "lm")){
    VarImp.lm(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "nnet")){
    .predict.nnet(x, ...)


  } else if(inherits(x, "glm")){
    .VarImp.glm(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "rpart")){
    .VarImp.rpart(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "rf")){
    tuneRF(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "nnet")){
    .VarImp.nnet(x, scale=scale, sort=sort, type=type, ...)

  } else if(inherits(x, "C5.0")){
    .VarImp.C5.0(x, scale=scale, sort=sort, type=type, ...)

  } else {
    NextMethod(x, sort=sort, type=type, ...)
  }

}





print.Tune <- function(x, digits = 2, ...) {
  # print as data.frame if something was changed
  print(round(data.frame(x[]), digits=digits))
}




