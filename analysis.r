#!/usr/bin/Rscript
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
judgments$n <- factor(judgments$n)

time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
mean(time$duration)
sd(time$duration) / sqrt(nrow(time))



## binned_pp_check(model, judgments)
##
##  Compute/display pp_checks for a model for different levels of confidence.
##  Uses type='hist' to avoid problems with bandwidth
##
binned_pp_check <- function(model, judgments) {
    for (c in levels(judgments$condition)) {
        for (n in unique(judgments$n)) {
            j <- judgments[judgments$n == n & judgments$condition == c, ]
            
            if (!is.null(nrow(j)) && nrow(j) > 1) {
                print(pp_check(model, bw=0.01, newdata=j) + theme_bw() +
                      ggtitle(sprintf("Ratings (%s cause(s), %s)", n, c)))
            }
        }
    }
}


mNormal <- brm(bf(rating ~ n * condition + (1 + n*condition |v| vignette),
                  sigma ~ n * condition +  (1 + n*condition |v| vignette)),
               prior=c(set_prior('normal(0, 10.0)', class='b'),
                       set_prior('normal(0, 10.0)', class='b', dpar='sigma')),
               data=judgments, file='mNormal_factor', inits="0",
               warmup=1500, iter=3000, cores=4,
               control=list(adapt_delta=0.99))

summary(mNormal)
pdf("Normal_factor.pdf")
plot(mNormal)
pp_check(mNormal, bw=0.01) + theme_bw()
binned_pp_check(mNormal, judgments)
marginal_effects(mNormal, resp='rating', probs=c(0.05, 0.95))
dev.off()



mZOIB <- brm(bf(rating ~ n * condition + (1 + n*condition |v| vignette),
                phi ~ n * condition +  (1 + n*condition |v| vignette),
                zoi ~ n * condition +  (1 + n*condition |v| vignette),
                coi ~ n * condition +  (1 + n*condition |v| vignette)),
             prior=c(set_prior('normal(0, 10.0)', class='b'),
                     set_prior('normal(0, 10.0)', class='b', dpar='phi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='zoi'),
                     set_prior('normal(0, 10.0)', class='b', dpar='coi')),
             data=judgments, family=zero_one_inflated_beta(), cores=4,
             file='mZOIB_factor', inits="0", control=list(adapt_delta=0.95))

summary(mZOIB)
pdf("ZOIB_factor.pdf")
plot(mZOIB)
pp_check(mZOIB, bw=0.01) + theme_bw()
binned_pp_check(mZOIB, judgments)
marginal_effects(mZOIB, resp='rating', probs=c(0.05, 0.95))
dev.off()
