## -----------------------------------------------------------------------------------
## Demo file for binomial models; start with 'demo(num_mod)'
## -----------------------------------------------------------------------------------


r.rlm <- FitMod(rcat_x ~ N(image_x) + N(custsat_x) + N(lowprem_x) + age_n +
                  sex_c + hins_x, d.nps, fitfn = "rlm")

r.lm <- FitMod(resp_n ~ N(image_x) + N(custsat_x) + N(lowprem_x) + age_n +
                 sex_c + hins_x, d.nps, fitfn = "lm")

r.pois <- FitMod(resp_n ~ N(image_x) + N(custsat_x) + N(lowprem_x) + age_n +
                   sex_c + hins_x, d.nps, fitfn = "poisson")


Conf(cut(predict(r.lm), breaks=c(-Inf, 7,9, Inf),
         labels=c("unlikely", "undecided","likely")),
     cut(model.frame(r.lm)$resp_n, breaks=c(-Inf, 7,9, Inf),
         labels=c("unlikely", "undecided","likely")))

