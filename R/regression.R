

# Same for regression



# Regr <- function(formula, data, ..., subset, na.action=na.pass, fitfn=NULL){
#
#   fitfn_null <- FALSE
#
#   if(is.null(fitfn)) {
#     fitfn <- "lm"
#     fitfn_null <- TRUE
#   }
#
#   cl <- match.call()
#
#   if(fitfn == "lm"){
#     fun <- get("lm", asNamespace("stats"), inherits = FALSE)
#
#   } else if(fitfn == "poisson"){
#     fun <- get("glm", asNamespace("stats"), inherits = FALSE)
#     cl[["family"]] <- "poisson"
#
#   } else if(fitfn == "quasipoisson"){
#     fun <- get("glm", asNamespace("stats"), inherits = FALSE)
#     cl[["family"]] <- "quasipoisson"
#
#   } else if(fitfn == "gamma"){
#     fun <- get("glm", asNamespace("stats"), inherits = FALSE)
#     if(is.null(cl[["link"]]))
#       cl[["family"]] <- "Gamma"
#     else {
#       cl[["family"]] <- eval(parse(text = gettextf("function() Gamma(link=%s)", cl[["link"]])), parent.frame())
#       cl[["link"]] <- NULL
#     }
#
#   } else if(fitfn == "negbin"){
#     if (!requireNamespace("MASS", quietly = TRUE))
#       stop("package 'MASS' must be installed")
#     if(is.null(cl[["theta"]]))
#       # no theta available, use glm.nb and fit theta
#       fun <- get("glm.nb", asNamespace("MASS"), inherits = FALSE)
#     else {
#       # theta is given, use usual glm and define family
#       fun <- glm
#       cl[["family"]] <- eval(parse(text = gettextf("function() negative.binomial(theta=%s)", cl[["theta"]])), parent.frame())
#       cl[["theta"]] <- NULL
#     }
#
#   } else if(fitfn == "polr"){
#     if (!requireNamespace("MASS", quietly = TRUE))
#       stop("package 'MASS' must be installed")
#     fun <- get("polr", asNamespace("MASS"), inherits = FALSE)
#
#   } else if(fitfn == "lmrob"){
#     if (!requireNamespace("robustbase", quietly = TRUE))
#       stop("package 'robustbase' must be installed")
#     fun <- get("lmrob", asNamespace("robustbase"), inherits = FALSE)
#
#   } else if(fitfn == "tobit"){
#     fun <- Tobit
#
#   } else if(fitfn == "zeroinfl"){
#     if (!requireNamespace("pscl", quietly = TRUE))
#       stop("package 'pscl' must be installed")
#     fun <- get("zeroinfl", asNamespace("pscl"), inherits = FALSE)
#
#   }
#
#   if(!fitfn_null)
#     cl[["fitfn"]] <- NULL
#
#   # cl[[1]] <- as.name(fun)
#   cl[[1]] <- fun
#
#   res <- eval(cl, parent.frame())
#
#   if(fitfn %in% c("poisson","gamma", "quasipoisson"))
#     res[["call"]][[1]] <- as.name("glm")
#
#   else if(fitfn == "negbin") {
#     if (is.null(cl[["family"]]))
#       res[["call"]][[1]] <- as.name("glm.nb")
#     else
#       res[["call"]][[1]] <- as.name("glm")
#
#   }
#
#   else if(fitfn %in% c("lm","lmrob","polr"))
#     res[["call"]][[1]] <- as.name(fitfn)
#
#
#   class(res) <- c("Regr", class(res))
#   return(res)
#
# }
#
#
#
# .drop1.multinom <-  function (object, scope, test = c("Chisq","none"), ...)
#   {
#     if (!inherits(object, "multinom"))
#       stop("Not a multinom fit")
#     if (missing(scope))
#       scope <- drop.scope(object)
#     else {
#       if (!is.character(scope))
#         scope <- attr(terms(update.formula(object, scope)),
#                       "term.labels")
#       if (!all(match(scope, attr(object$terms, "term.labels"),
#                      nomatch = FALSE)))
#         stop("scope is not a subset of term labels")
#     }
#     ns <- length(scope)
#     ans <- matrix(nrow = ns + 1, ncol = 2, dimnames = list(c("<none>",
#                                                              scope), c("Df", "AIC")))
#     ans[1, ] <- c(object$edf, object$AIC)
#     if (test[1]=="Chisq") ans <- cbind(ans,Chisq=NA, p.value=NA)
#     env <- environment(formula(object))
#     for(i in seq(ns)) {
#       tt <- scope[i]
#       ##        cat("trying -", tt, "\n")
#       nfit <- update(object, as.formula(paste("~ . -", tt)),
#                      evaluate = FALSE)
#       nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)
#       ##-        nobject <- update(object, paste("~ . -", tt))
#       if (nfit$edf == object$edf)
#         nfit$AIC <- NA
#       ans[i+1, ] <- c(nfit$edf, nfit$AIC,
#                       if (test[1]=="Chisq") unlist(anova(object,nfit)[2,6:7]))
#     }
#     as.data.frame(ans)
#   }
