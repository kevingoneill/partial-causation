#!/usr/bin/Rscript
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidybayes)
library(modelr)

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) {
    writeLines('Usage: ./analisis.r <file-name>')
    quit()
}


judgments <- read.csv(args[1], header=TRUE)
judgments$rating <- judgments$rating / 100
judgments$confidence <- judgments$confidence / 100
judgments$n <- as.factor(judgments$n)
judgments$normality <- factor(judgments$normality, levels=c('N', 'A'))

length(unique(judgments$id))


time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
writeLines(sprintf("%f, %f\n",
                   mean(time$duration),
                   sd(time$duration) / sqrt(nrow(time))))





mZOIB <- brm(bf(rating ~ n*structure*normality + (1 |v| vignette),
                phi ~ n*structure*normality + (1 |v| vignette),
                zoi ~ n*structure*normality + (1 |v| vignette),
                coi ~ n*structure*normality + (1 |v| vignette)),
             prior=c(set_prior('normal(0, 10.0)'),
                     set_prior('normal(0, 10.0)', dpar='phi'),
                     set_prior('normal(0, 10.0)', dpar='zoi'),
                     set_prior('normal(0, 10.0)', dpar='coi')),
             data=judgments, family=zero_one_inflated_beta(),
             file='mZOIB', inits="0",
             cores=4, control=list(adapt_delta=0.99))
summary(mZOIB)

pdf("ZOIB.pdf")
plot(mNormal)
pp_check(mNormal, bw=0.01) + theme_bw()
conditional_effects(mNormal, 'structure:normality',
                    conditions=make_conditions(judgments, c('n')))
#conditional_effects(mNormal, 'structure:normality', dpar='sigma',
#                    conditions=make_conditions(judgments, c('n')))
dev.off()

quit()


## Fit a multivariate model (for rating and confidence)
## allowing for unqual variances (sigma),
## with uncorrelated random intercepts/slopes per vignette
mNormal <- brm(bf(rating ~ n*structure*normality + (1 |v| vignette),
                  sigma ~ n*structure*normality + (1 |v| vignette)) +
               bf(confidence ~ n*structure*normality + (1 |v| vignette),
                  sigma ~ n*structure*normality + (1 |v| vignette)),
               prior=c(set_prior('normal(0, 10.0)', resp='rating'),
                       set_prior('normal(0, 10.0)', resp='rating', dpar='sigma'),
                       set_prior('normal(0, 10.0)', resp='confidence'),
                       set_prior('normal(0, 10.0)', resp='confidence', dpar='sigma')),
               data=judgments, file='mNormal_multi_reduced', inits="0",
               cores=4, control=list(adapt_delta=0.99))
#summary(mNormal)


## Plot model fits/diagnostics
#pdf("Normal_multi.pdf")
#plot(mNormal)
#pp_check(mNormal, bw=0.01, resp='rating') + theme_bw()
#pp_check(mNormal, bw=0.01, resp='confidence') + theme_bw()
#conditional_effects(mNormal, 'structure:normality', resp='rating',
#                    conditions=make_conditions(judgments, c('n')))
#conditional_effects(mNormal, 'structure:normality', resp='rating', dpar='sigma',
#                    conditions=make_conditions(judgments, c('n')))
#conditional_effects(mNormal, 'structure:normality', resp='confidence',
#                    conditions=make_conditions(judgments, c('n')))
#conditional_effects(mNormal, 'structure:normality', resp='confidence', dpar='sigma',
#                    conditions=make_conditions(judgments, c('n')))
#dev.off()


## Get posterior estimates over our data
ratings <- judgments %>% data_grid(n, structure, normality) %>%
    add_fitted_draws(mNormal, re_formula=NA, resp='rating', dpar='sigma')
confidence <- judgments %>% data_grid(n, structure, normality) %>%
    add_fitted_draws(mNormal, re_formula=NA, resp='confidence', dpar='sigma')


## Do overall contrasts for levels of N, collapsing structure and normality
ratings %>%
    compare_levels(.value, by=normality, comparison='ordered') %>%
    group_by(normality, structure) %>%
    median_hdi(.value)
quit()
ratings %>%
    compare_levels(.value, by=normality) %>%
    group_by(structure) %>%
    median_hdi(.value)
ratings %>%
    compare_levels(.value, by=normality) %>%
    group_by(structure) %>%
    summarise(median=quantile(.value, probs=c(0.5)),
              lower=quantile(.value, probs=c(0.05)),
              upper=quantile(.value, probs=c(0.95)))
quit()

ratings %>%
    compare_levels(sigma, by=normality) %>%
    median_hdi(sigma)

quit()

ratings %>% compare_levels(.value, by=n, comparison='ordered') %>%
    group_by(n) %>%
    median_hdi(.value)
ratings %>% compare_levels(sigma, by=n, comparison='ordered') %>%
    group_by(n) %>%
    median_hdi(sigma)

ratings %>% compare_levels(.value, by=n, comparison='ordered') %>%
    group_by(n, structure, normality) %>%
    median_hdi(.value)
ratings %>% compare_levels(sigma, by=n, comparison='ordered') %>%
    group_by(n, structure, normality) %>%
    median_hdi(sigma)

confidence %>% compare_levels(.value, by=n, comparison='ordered') %>%
    group_by(n) %>%
    median_hdi(.value)
confidence %>% compare_levels(sigma, by=n, comparison='ordered') %>%
    group_by(n) %>%
    median_hdi(sigma)


## Plot estimates/contrasts for Rating
ratings %>%
    ggplot(aes(x=structure, y=.value,
               group=interaction(normality, structure), fill=normality)) +
    ylab('Causal Rating') +
    ylim(0, 1.0) +
    geom_violin(aes(color=normality), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    facet_wrap(n ~ .) + theme_light()
ratings %>%
    compare_levels(.value, by=normality) %>%
    ggplot(aes(x=n, y=.value, group=interaction(n, structure), fill=structure)) +
    geom_violin(aes(color=structure), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    ylab('Causal Rating Contrasts: Abnormal - Normal') +
    xlab('N') + theme_light()

## Plot estimates/contrasts for sigma(Rating)
ratings %>%
    ggplot(aes(x=structure, y=sigma,
               group=interaction(normality, structure), fill=normality)) +
    ylab('sigma(Causal Rating)') +
    ylim(0, 1.0) +
    geom_violin(aes(color=normality), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    facet_wrap(n ~ .) + theme_light()
ratings %>%
    compare_levels(sigma, by=normality) %>%
    ggplot(aes(x=n, y=sigma, group=interaction(n, structure), fill=structure)) +
    geom_violin(aes(color=structure), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    ylab('sigma(Causal Rating) Contrasts: Abnormal - Normal') +
    xlab('N') + theme_light()

## Plot estimates/contrasts for Confidence
confidence %>%
    ggplot(aes(x=structure, y=.value,
               group=interaction(normality, structure), fill=normality)) +
    ylab('Confidence') +
    ylim(0, 1.0) +
    geom_violin(aes(color=normality), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    facet_wrap(n ~ .) + theme_light()
confidence %>%
    compare_levels(.value, by=normality) %>%
    ggplot(aes(x=n, y=.value, group=interaction(n, structure), fill=structure)) +
    geom_violin(aes(color=structure), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    ylab('Confidence Contrasts: Abnormal - Normal') +
    xlab('N') + theme_light()

## Plot estimates/contrasts for sigma(Confidence)
confidence %>%
    ggplot(aes(x=structure, y=sigma,
               group=interaction(normality, structure), fill=normality)) +
    ylab('sigma(Confidence)') +
    ylim(0, 1.0) +
    geom_violin(aes(color=normality), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    facet_wrap(n ~ .) + theme_light()
confidence %>%
    compare_levels(sigma, by=normality) %>%
    ggplot(aes(x=n, y=sigma, group=interaction(n, structure), fill=structure)) +
    geom_violin(aes(color=structure), position=position_dodge(width=1)) +
    stat_summary(fun.data=median_hdi, geom='pointrange',
                 position=position_dodge(width=1)) +
    ylab('sigma(Confidence) Contrasts: Abnormal - Normal') +
    xlab('N') + theme_light()
