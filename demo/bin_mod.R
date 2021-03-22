## -----------------------------------------------------------------------------------
## Demo file for binomial models; start with 'demo(bin_mod)'
## -----------------------------------------------------------------------------------

set.seed(453)
d.pim <- SplitTrainTest(d.pima, p = 0.2)
mdiab <- formula(diabetes ~ pregnant + glucose + pressure + triceps
                 + insulin + mass + pedigree + age)

r.glm <- FitMod(mdiab, data=d.pim$train, fitfn="logit")
r.rp  <- FitMod(mdiab, data=d.pim$train, fitfn="rpart")
r.rf  <- FitMod(mdiab, data=d.pim$train, fitfn="randomForest")
r.svm <- FitMod(mdiab, data=d.pim$train, fitfn="svm")
r.c5  <- FitMod(mdiab, data=d.pim$train, fitfn="C5.0")
r.nn  <- FitMod(mdiab, data=d.pim$train, fitfn="nnet", decay=0.2, size=7, maxit=1000)
r.nb  <- FitMod(mdiab, data=d.pim$train, fitfn="naive_bayes")
r.lda <- FitMod(mdiab, data=d.pim$train, fitfn="lda")
r.qda <- FitMod(mdiab, data=d.pim$train, fitfn="qda")
r.lb  <- FitMod(mdiab, data=d.pim$train, fitfn="lb")

mods <- list(glm=r.glm, rp=r.rp, rf=r.rf, svm=r.svm, c5=r.c5,
             nn=r.nn, nb=r.nb, lda=r.lda, qda=r.qda, lb=r.lb)

# insight in the Regression tree
plot(r.rp, box.palette = as.list(Pal("Helsana", alpha = 0.5)))

# Insample accuracy ...
TModC(mods, ord="auc")
# ... is substantially different from the out-of-bag:

TModC(mods, newdata=d.pim$test, reference=d.pim$test$diabetes, ord="b")
# C5 and SVM turn out to be show-offs! They overfit quite ordinary
# whereas randomforest and logit keep their promises. ...

sapply(mods, function(z) VarImp(z))


# undebug(Tune)
# zz <- ModTools::Tune(FitMod(mdiab, data=d.pim$train, fitfn="nnet"),
#     decay=seq(0.1, 1.0, 0.05), size=c(5:15))

# zz$modpar


