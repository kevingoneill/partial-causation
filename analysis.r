library(brms)
library(ggplot2)

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) {
    writeLines('Usage: ./analisis.r <file-name>')
    quit()
}


judgments <- read.csv(args[1], header=TRUE)
judgments$rating <- judgments$rating / 100
judgments$confidence <- judgments$confidence / 100
judgments$condition <- as.factor(judgments$condition)

time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
mean(time$duration)
sd(time$duration) / sqrt(nrow(time))



mNormal <- brm(bf(rating ~ n * condition + (1 |v| vignette),
                  sigma ~ n * condition +  (1 |v| vignette)),
               prior=c(set_prior('normal(0, 10.0)', class='b'),
                       set_prior('normal(0, 10.0)', class='b', dpar='sigma')),
               data=judgments, file='mNormal', inits="0",
               warmup=1500, iter=3000,
               control=list(adapt_delta=0.99))

summary(mNormal)
pdf("Normal.pdf")
plot(mNormal)
pp_check(mNormal, bw=0.005) + theme_bw()
marginal_effects(mNormal, resp='rating', probs=c(0.05, 0.95))
#marginal_effects(mNormal, resp='confidence', probs=c(0.05, 0.95))
dev.off()



mZOIB <- brm(bf(rating ~ n * condition + (1 |v| vignette),
                phi ~ n * condition +  (1 |v| vignette),
                zoi ~ n * condition +  (1 |v| vignette),
                coi ~ n * condition +  (1 |v| vignette)),
             prior=c(set_prior('normal(0, 10.0)', class='b'),
                     set_prior('normal(0, 10.0)', class='b', dpar='phi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='zoi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='coi')),
             data=judgments, family=zero_one_inflated_beta(),
             file='mZOIB', inits="0", control=list(adapt_delta=0.95))

summary(mZOIB)
pdf("ZOIB.pdf")
plot(mZOIB)
pp_check(mZOIB, bw=0.005) + theme_bw()
marginal_effects(mZOIB, resp='rating', probs=c(0.05, 0.95))
#marginal_effects(mZOIB, resp='confidence', probs=c(0.05, 0.95))
dev.off()
