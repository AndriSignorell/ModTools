\name{ModTools-package}
\alias{ModTools-package}
\alias{ModTools}
\docType{package}
\title{Regression and Classification Tools

}
\description{
There is a rich selection of R packages implementing algorithms for classification and regression tasks out there. The authors legitimately take the liberty to tailor the function interfaces according to their own taste  and needs. For us other users, however, this often results in struggling with user interfaces, some of which are - to put it mildly - rather weird and almost always different in terms of arguments and result structures.
\bold{ModTools} pursues the goal of offering uniform handling for at least the most important regression and classification models.\cr
The function \code{FitMod()} is designed as simple and consistent interface to these original functions while maintaining the flexibility to pass on all possible arguments. \code{print}, \code{plot}, \code{summary} and \code{predict} operations can so be carried out consistently. The results will again be reshaped to a reasonable standard.

For all the functions of this package Google styleguides are used as naming rules (in absence of convincing alternatives). The 'CamelStyle' has been consequently applied to functions borrowed from contributed R packages as well.

As always: Feedback, feature requests, bugreports and other suggestions are welcome!
}

\section{Warning}{This package is still under development. You should be aware that everything in the package might be subject to change. Backward compatibility is not yet guaranteed. Functions may be deleted or renamed and new syntax may be inconsistent with earlier versions. By release of version 1.0 the "deprecated-defunct process" will be installed.
}


\details{
ModTools comprises interfaces to the following models: \cr

\tabular{ll}{
\cr\cr
\bold{Regression}:\tab\cr
 \code{\link{lm}()}  \tab 	 Linear model OLS\cr
 \code{\link[robustbase]{lmrob}()}  \tab 	 Robust linear model (package \bold{robustbase})\cr
 \code{poisson()}  \tab 	GLM model with family \code{poisson}\cr
 \code{negbin()}  \tab 	GLM model with family \code{negative.binomial} (package  \bold{MASS})\cr
 \code{gamma()}  \tab 	Gamma model \code{\link{glm}}, family \code{gamma}\cr
 \code{tobit()}  \tab 	Tobit model for censored responses (package  \bold{AER})\cr
\cr\cr\cr
\bold{Classification}:	\tab\cr
\code{\link{lda}()}  \tab 	 Linear discriminant analysis (package  \bold{MASS})\cr
\code{\link{qda}()}  \tab 	 Quadratic discriminant analysis (package  \bold{MASS})\cr
\code{logit()}  \tab 	 Logistic Regression model (\code{\link{glm}}, family \code{binomial(logit)})\cr
\code{\link{multinom}()}  \tab 	 Multinomial Regression model (package  \bold{nnet})\cr
 \code{\link[MASS]{polr}()}  \tab 	 Proportional odds model (package \bold{MASS})\cr
\code{\link{rpart}()}  \tab 	 Regression and classification trees (package  \bold{rpart})\cr
\code{\link[nnet]{nnet}()}  \tab 	 Neuronal networks (package  \bold{nnet})\cr
\code{\link[randomForest]{randomForest}()}   \tab 	 Random forests (package  \bold{randomForest})\cr
\code{\link[C50]{C5.0}()}   \tab 	 C5.0 tree (package  \bold{C50})\cr
\code{\link[e1071]{svm}()}   \tab 	 Support vector machines (package  \bold{e1071})\cr
\code{\link[e1071]{naiveBayes}()}   \tab 	 Naive Bayes classificator (package  \bold{e1071})\cr
\code{\link{LogitBoost}()}  \tab 	 Logit boost (using decision
  stumps as weak learners)\cr
\cr\cr\cr \bold{Other interfaced functions}:\tab\cr
\code{\link[class]{knn}()}  \tab 	 k nearest neighbours (package  \bold{class})\cr
\code{\link[class]{SOM}()}  \tab 	 Self organizing maps (package  \bold{class})\cr
% \link{BestTree}  \tab 	 Variable importance for many models\cr
% \link{PlotHistLeaf}  \tab 	 Tuning for many models\cr
% \link{PlotTree}  \tab 	 Confusion matrix for classification tasks\cr
% \link{LeafRates}  \tab 	 Confusion matrix for classification tasks\cr
% \link{Purity}  \tab 	 Confusion matrix for classification tasks\cr
% \link{Snip}  \tab 	 Variable importance for many models\cr
% \link{Flop}  \tab 	 Variable importance for many models\cr
% \link{Improve}  \tab 	 Variable importance for many models\cr
% \link{Desc.rpart}  \tab 	 Variable importance for many models\cr
% \link{Splits}  \tab 	 Variable importance for many models\cr
% \link{Rules}  \tab 	 Variable importance for many models\cr
% \link{TreeDepth}  \tab 	 Variable importance for many models\cr
\cr\cr\cr \bold{Preprocess}:\tab\cr
\code{\link{SplitTrainTest}()}  \tab 	Splits a data frame or index vector in train an test sample\cr
\cr\cr\cr \bold{Prediction and Validation}:\tab\cr
\code{\link{predict}()}  \tab 	 Consistent predict for \code{FitMod} models\cr
\code{\link{VarImp}()}  \tab 	 Variable importance for most \code{FitMod} models\cr
\code{\link{ROC}()}  \tab ROC curves for all dichotomous classification \code{FitMod} models	 \cr
\code{\link{TModC}()}  \tab Aggregated results for multiple \code{FitMod} classification models \cr
%\code{\link{Tune}}  \tab Tuning approaches for \code{FitMod} models	 \cr
%\code{\link{CV}}  \tab Crossvalidating approaches for \code{FitMod} models	 \cr
\cr\cr\cr
\bold{Tests}:	\tab\cr
\code{\link{BreuschPaganTest}()}  \tab 	 Breusch-Pagan test against heteroskedasticity.\cr
%\code{\link{CoefTest}}  \tab 	 Coefficient test\cr
}

We might want to add support later down the road for:
mixed models, gam, mars, lasso\cr
ART Adaptive Resonance Theory cluster, association rules, multidimensional scaling

preprocess:
1. Check for near zero variance predictors and flag as near zero if:
* the percentage of unique values is less than 20
* the ratio of the most frequent to the second most frequent value is greater than 20,

transformations, retransformation: rcs - restricted cubic splines smearingEstimator
imputation: ... ? knnImputation)

}

\author{
Andri Signorell\cr
Helsana Versicherungen AG, Health Sciences, Zurich\cr
HWZ University of Applied Sciences in Business Administration Zurich.\cr

Includes R source code and/or documentation previously published by (in alphabetical order): \cr

Bernhard Compton, Ulrike Groemping, Max Kuhn, Werner A. Stahel
Michal Majka, <michalmajka@hotmail.com>
Jarek Tuszynski (SAIC) <jaroslaw.w.tuszynski@saic.com>
Achim Zeileis <Achim.Zeileis@R-project.org>

Special thanks go to Beat Bruengger for their valuable contributions and testing.

The good things come from all these guys, any problems are likely due to my tweaking.
Thank you all! \cr


Maintainer: Andri Signorell <andri@signorell.net>

}

\keyword{ package }


\examples{

r.swiss <- FitMod(Fertility ~ ., swiss, fitfn="lm")
r.swiss
# PlotTA(r.swiss)
# PlotQQNorm(r.swiss)


## Count models

data(housing, package="MASS")

# poisson count
r.pois <- FitMod(Freq ~ Infl*Type*Cont + Sat, family=poisson, data=housing, fitfn="poisson")

# negative binomial count
r.nb <- FitMod(Freq ~ Infl*Type*Cont + Sat, data=housing, fitfn="negbin")
summary(r.nb)

r.log <- FitMod(log(Freq) ~ Infl*Type*Cont + Sat, data=housing, fitfn="lm")
summary(r.log)

r.ols <- FitMod(Freq ~ Infl*Type*Cont + Sat, data=housing, fitfn="lm")
summary(r.ols)

r.gam <- FitMod(Freq ~ Infl*Type*Cont + Sat, data=housing, fitfn="gamma")
summary(r.gam)

r.gami <- FitMod(Freq ~ Infl*Type*Cont + Sat, data=housing, fitfn="gamma", link="identity")
summary(r.gami)

options(digits=3)
TMod(r.pois, r.nb, r.log, r.ols, r.gam, r.gami)



## Ordered Regression

r.polr <- FitMod(Sat ~ Infl + Type + Cont, data=housing, fitfn="polr", weights = Freq)

# multinomial Regression
r.mult <- FitMod(factor(Sat, ordered=FALSE) ~ Infl + Type + Cont, data=housing,
                   weights = housing$Freq, fitfn="multinom")


# Regression tree
r.rp <- FitMod(factor(Sat, ordered=FALSE) ~ Infl + Type + Cont, data=housing,
                 weights = housing$Freq, fitfn="rpart")

# compare predictions
d.p <- expand.grid(Infl=levels(housing$Infl), Type=levels(housing$Type), Cont=levels(housing$Cont))
d.p$polr <- predict(r.polr, newdata=d.p)
# ??
# d.p$ols <- factor(round(predict(r.ols, newdata=d.p)^2), labels=levels(housing$Sat))
d.p$mult <- predict(r.mult, newdata=d.p)
d.p$rp <- predict(r.rp, newdata=d.p, type="class")

d.p


# Classification with 2 classes  ***************

r.pima <- FitMod(diabetes ~ ., d.pima, fitfn="logit")
r.pima
Conf(r.pima)
plot(ROC(r.pima))
OddsRatio(r.pima)


# rpart tree
rp.pima <- FitMod(diabetes ~ ., d.pima, fitfn="rpart")
rp.pima
Conf(rp.pima)
lines(ROC(rp.pima), col=hblue)
# to be improved
plot(rp.pima, col=SetAlpha(c("blue","red"), 0.4), cex=0.7)


# Random Forest
rf.pima <- FitMod(diabetes ~ ., d.pima, method="class", fitfn="randomForest")
rf.pima
Conf(rf.pima)
lines(ROC(r.pima), col=hred)



# more models to compare

d.pim <- SplitTrainTest(d.pima, p = 0.2)
mdiab <- formula(diabetes ~ pregnant + glucose + pressure + triceps
                            + insulin + mass + pedigree + age)

r.glm <- FitMod(mdiab, data=d.pim$train, fitfn="logit")
r.rp <- FitMod(mdiab, data=d.pim$train, fitfn="rpart")
r.rf <- FitMod(mdiab, data=d.pim$train, fitfn="randomForest")
r.svm <- FitMod(mdiab, data=d.pim$train, fitfn="svm")
r.c5 <- FitMod(mdiab, data=d.pim$train, fitfn="C5.0")
r.nn <- FitMod(mdiab, data=d.pim$train, fitfn="nnet")
r.nb <- FitMod(mdiab, data=d.pim$train, fitfn="naiveBayes")
r.lda <- FitMod(mdiab, data=d.pim$train, fitfn="lda")
r.qda <- FitMod(mdiab, data=d.pim$train, fitfn="qda")
r.lb <- FitMod(mdiab, data=d.pim$train, fitfn="lb")

mods <- list(glm=r.glm, rp=r.rp, rf=r.rf, svm=r.svm, c5=r.c5
             , nn=r.nn, nb=r.nb, lda=r.lda, qda=r.qda, lb=r.lb)

# insight in the Regression tree
plot(r.rp, box.palette = as.list(Pal("Helsana", alpha = 0.5)))

# Insample accuracy ...
TModC(mods, ord="a")
# ... is substantially different from the out-of-bag:
TModC(mods, newdata=d.pim$test, reference=d.pim$test$diabetes, ord="b")
# C5 and SVM turn out to be show-offs! They overfit quite ordinary
# whereas randomforest and logit keep their promises. ...

sapply(mods, function(z) VarImp(z))


# Multinomial classification problem with n classes  ***************

d.gl <- SplitTrainTest(d.glass, p = 0.2)
mglass <- formula(Type ~ RI + Na + Mg + Al + Si + K + Ca + Ba + Fe)

r.mult <- FitMod(mglass, data=d.gl$train, maxit=600, fitfn="multinom")
r.rp <- FitMod(mglass, data=d.gl$train, fitfn="rpart")
r.rf <- FitMod(mglass, data=d.gl$train, fitfn="randomForest")
r.svm <- FitMod(mglass, data=d.gl$train, fitfn="svm")
r.c5 <- FitMod(mglass, data=d.gl$train, fitfn="C5.0")
set.seed(1984)
r.nn <- FitMod(mglass, data=d.gl$train, fitfn="nnet")
r.nbay <- FitMod(mglass, data=d.gl$train, fitfn="naiveBayes")
r.lda <- FitMod(mglass, data=d.gl$train, fitfn="lda")
# r.qda <- FitMod(mglass, data=d.glass, fitfn="qda")
r.lb <- FitMod(mglass, data=d.gl$train, fitfn="lb")

mods <- list(glm=r.mult, rp=r.rp, rf=r.rf, svm=r.svm, c5=r.c5,
             nn=r.nn, nbay=r.nbay, lda=r.lda, lb=r.lb)

# confusion matrix and other quality measures can be calculated with Conf()
Conf(r.rf)

# we only extract the general accuracy
sapply(lapply(mods, function(z) Conf(z)), "[[", "acc")

# let's compare r.mult with a model without RI as predictor
Conf(r.mult)
Conf(update(r.mult, . ~ . -RI))
}
