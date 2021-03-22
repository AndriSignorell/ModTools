## -----------------------------------------------------------------------------------
## Demo file for binomial models; start with 'demo(mult_mod)'
## -----------------------------------------------------------------------------------

r.mult <- FitMod(rcat_x ~ N(image_x) + N(custsat_x) + N(lowprem_x) + age_n +
                   sex_c , data=d.nps, fitfn = "multinom")
r.polr <- FitMod(rcat_x ~ N(image_x) + N(custsat_x) + N(lowprem_x) + age_n +
                   sex_c , data=d.nps, fitfn = "polr")
r.rf <- FitMod(rcat_x ~ (image_x) + (custsat_x) + (lowprem_x) + age_n +
                 sex_c , data=d.nps, fitfn = "randomForest")

TModC(r.polr, r.mult)



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
r.nbay <- FitMod(mglass, data=d.gl$train, fitfn="naive_bayes")
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

