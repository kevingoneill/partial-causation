#!/usr/bin/Rscript
library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tidybayes)
library(modelr)
library(emmeans)
library(bayestestR)
library(wesanderson)
library(ggisoband)
library(ggh4x)
library(patchwork)


## Set a color palette for plotting
PALETTE <- rev(wes_palette("Darjeeling1", n=2))

## Read in and normalize data
judgments <- read.csv('data/processed_data.csv', header=TRUE)
judgments$rating <- judgments$rating / 100
judgments$confidence <- judgments$confidence / 100
judgments$n <- as.factor(judgments$n)
judgments$normality <- factor(judgments$normality, levels=c('N', 'A'))

## Print descriptives
writeLines(sprintf("Number of subjects: %d",
                   length(unique(judgments$id))))
writeLines(sprintf("Number of vignettes: %d",
                   length(unique(judgments$vignette))))
writeLines("Design: 2 (structure) x 2 (normality) * 4 (N)")
length(unique(judgments$id)) /
    (length(unique(judgments$vignette)) * 2 * 2 * 4)
judgments %>% select(sex) %>% table
writeLines(sprintf("Average age: %f (%f)\n", mean(judgments$age, na.rm=TRUE),
                   sd(judgments$age, na.rm=TRUE)))
writeLines(sprintf("Average duration: %f (%f)\n", mean(judgments$duration),
                   sd(judgments$duration)))

## Fit a multivariate model (for rating and confidence)
## allowing for unqual variances (sigma),
## with correlated random intercepts/slopes per vignette
mMulti <- brm(bf(rating ~ n*structure*normality + (1 |v| vignette),
                 sigma ~ n*structure*normality + (1 |v| vignette)) +
              bf(confidence ~ n*structure*normality + (1 |v| vignette),
                 sigma ~ n*structure*normality + (1 |v| vignette)),
              prior=c(set_prior('normal(0, 1.0)'),
                      set_prior('normal(0, 5.0)', resp='rating', dpar='sigma'),
                      set_prior('normal(0, 1.0)', resp='confidence'),
                      set_prior('normal(0, 5.0)', resp='confidence', dpar='sigma')),
              data=judgments, file='mMulti', iter=5000,
              inits="0", sample_prior="yes", save_all_pars=TRUE,
              cores=4, control=list(adapt_delta=0.99))
summary(mMulti, priors=TRUE)


#########################################################################################
##
## Descriptives:
##    Perform model testing on alternate models
##
mReduced <- brm(bf(rating ~ n*structure*normality + (1 |v| vignette)) +
                bf(confidence ~ n*structure*normality + (1 |v| vignette)),
                prior=c(set_prior('normal(0, 1.0)'),
                        set_prior('normal(0, 1.0)', resp='confidence')),
                data=judgments, file='mReduced', iter=5000,
                inits="0", sample_prior="yes", save_all_pars=TRUE,
                cores=4, control=list(adapt_delta=0.99))

mZOIB <- brm(bf(rating ~ n*structure*normality + (1 |v| vignette),
                phi ~ n*structure*normality + (1 |v| vignette),
                zoi ~ n*structure*normality + (1 |v| vignette),
                coi ~ n*structure*normality + (1 |v| vignette)) +
             bf(confidence ~ n*structure*normality + (1 |v| vignette),
                phi ~ n*structure*normality + (1 |v| vignette),
                zoi ~ n*structure*normality + (1 |v| vignette),
                coi ~ n*structure*normality + (1 |v| vignette)),
             prior=c(set_prior('normal(0, 10.0)', resp='rating', class='b'),
                     set_prior('normal(0, 10.0)', resp='rating',
                               dpar='phi', class='b'),
                     set_prior('normal(0, 10.0)', resp='rating',
                               dpar='zoi', class='b'),
                     set_prior('normal(0, 10.0)', resp='rating',
                               dpar='coi', class='b'),
                     set_prior('normal(0, 10.0)', resp='confidence', class='b'),
                     set_prior('normal(0, 10.0)', resp='confidence',
                               dpar='phi', class='b'),
                     set_prior('normal(0, 10.0)', resp='confidence',
                               dpar='zoi', class='b'),
                     set_prior('normal(0, 10.0)', resp='confidence',
                               dpar='coi', class='b')),
             data=judgments, family=zero_one_inflated_beta(),
             file='mZOIB', inits="0", iter=5000, cores=4,
             save_all_pars=TRUE, sample_prior="yes",
             control=list(adapt_delta=0.99, max_treedepth=20))

LOO(mMulti, mReduced, mZOIB)
model_weights(mMulti, mReduced, mZOIB)


## Plot raw causal ratings
judgments %>%
    ggplot(aes(x=n, y=rating,
               group=normality, fill=normality)) +
    ylab('Causal Rating') + ylim(0, 1) + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.85), adjust=0.75,
             n=nrow(judgments), normalize='panels', geom='slab') +
    stat_pointinterval(position=position_dodge(0.85)) +
    scale_alpha_discrete(range = c(0.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("data-ratings.png", width=6, height=4)

## Plot raw confidence ratings
judgments %>%
    ggplot(aes(x=n, y=confidence,
               group=normality, fill=normality)) +
    ylab('Confidence') + ylim(0, 1) + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.85), adjust=0.75,
             n=nrow(judgments), normalize='panels', geom='slab') +
    stat_pointinterval(position=position_dodge(0.85)) +
    scale_alpha_discrete(range = c(0.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("data-confidence.png", width=6, height=4)

## Plot raw causal and confidence ratings
judgments %>%
    ggplot(aes(x=confidence, y=rating,
               color=normality, fill=normality)) +
    geom_density_bands(aes(alpha=stat(ndensity)), bins=10,
                       h=c(0.5, 0.5), show.legend=FALSE) +
    facet_nested(structure+normality ~ n,
                 labeller=labeller(n=c('1'='1 Cause', '2'='2 Causes',
                                       '3'='3 Causes', '4'='4 Causes'),
                                   normality=c('N'='Normal', 'A'='Abnormal'),
                                   structure=c('C'='Conjunctive',
                                               'D'='Disjunctive'))) +
    xlab('Confidence Rating') + ylab('Causal Rating') +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75),
                       labels=c('0', '.25', '.5', '.75'),
                       expand=expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75),
                       labels=c('0', '.25', '.5', '.75'),
                       expand=expansion(mult = c(0, 0))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          panel.spacing=unit(0, 'lines'),
          legend.position='none')
ggsave('data-bivariate.png', width=5, height=5)


## extract model predictions
constrain <- function (x) {pmax(0, pmin(1, x))}
predictions <- full_join(predicted_draws(mMulti, judgments, resp='rating',
                                         prediction='.prediction.rating'),
                         predicted_draws(mMulti, judgments, resp='confidence',
                                         prediction='.prediction.confidence'))

## Plot model predictions for causal ratings
predictions %>%
    ggplot(aes(x=n, y=constrain(.prediction.rating),
               group=normality, fill=normality)) +
    ylab('Causal Rating') + coord_cartesian(ylim=c(0, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             normalize='panels', geom='slab') +
    stat_pointinterval(aes(y=.prediction.rating), position=position_dodge(0.75)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ .,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("predictions-ratings.png", width=6, height=4)

## Plot model predictions for confidence
predictions %>%
    ggplot(aes(x=n, y=constrain(.prediction.confidence),
               group=normality, fill=normality)) +
    ylab('Causal Rating') + coord_cartesian(ylim=c(0, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             normalize='panels', geom='slab') +
    stat_pointinterval(aes(y=.prediction.confidence), position=position_dodge(0.75)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ .,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("predictions-confidence.png", width=6, height=4)

## Plot model predictions for causal ratings and confidence
predictions %>%
    filter(.draw == sample(.draw, length(.draw) * 0.5)) %>%
    ggplot(aes(x=constrain(.prediction.confidence),
               y=constrain(.prediction.rating),
               color=normality, fill=normality)) +
    geom_density_bands(aes(alpha=stat(ndensity)), bins=10,
                       h=c(0.45, 0.45),
                       show.legend=FALSE) +
    facet_nested(structure+normality ~ n,
                 labeller=labeller(n=c('1'='1 Cause', '2'='2 Causes',
                                       '3'='3 Causes', '4'='4 Causes'),
                                   normality=c('N'='Normal', 'A'='Abnormal'),
                                   structure=c('C'='Conjunctive',
                                               'D'='Disjunctive'))) +
    xlab('Predicted Confidence Rating') + ylab('Predicted Causal Rating') +
    scale_alpha_continuous(range = c(0, 1)) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75),
                       labels=c('0', '.25', '.5', '.75'),
                       expand=expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75),
                       labels=c('0', '.25', '.5', '.75'),
                       expand=expansion(mult = c(0, 0))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          panel.spacing=unit(0, 'lines'),
          legend.position='none')
ggsave('predictions-bivariate.png', width=5, height=5)


## Extract posterior samples
ratings <- judgments %>% data_grid(n, structure, normality) %>%
    add_fitted_draws(mMulti, re_formula=NA, resp='rating',
                     dpar='sigma', scale='linear')
confidence <- judgments %>% data_grid(n, structure, normality) %>%
    add_fitted_draws(mMulti, re_formula=NA, resp='confidence',
                     dpar='sigma', scale='linear')
rescor <- mMulti %>%
    spread_draws(rescor__rating__confidence) %>%
    rename(rescor=rescor__rating__confidence) %>%
    mutate(model='full')

## Define ROPE ranges as 0.1 * sd of the DV
## NOTE: ROPEs for SD are estimated using posteriors over the data
ROPE = sd(judgments$rating) * 0.1
sdROPE <- sd(add_fitted_draws(judgments, mMulti, resp='rating',
                              dpar='sigma', scale='linear')$sigma) * 0.1
ROPE.conf <- sd(judgments$confidence) * 0.1
sdROPE.conf <- sd(add_fitted_draws(judgments, mMulti, resp='confidence',
                                   dpar='sigma', scale='linear')$confidence) * 0.1
ROPE.rescor <- sd(rescor$rescor) * 0.1



#########################################################################################
##
## Residual Correlation:
##    Test/plot residual correlation between causal judgments and confidence
##
describe_posterior(rescor$rescor, ci=0.95, rope_ci=0.95,
                   rope_range=c(-ROPE.rescor, ROPE.rescor)) %>%
    select(-c(CI, ROPE_CI))
rescor %>%
    mutate(Parameter='Residual Correlation') %>%
    ggplot(aes(y=rescor, x=model, 
               fill=stat(abs(y) > ROPE.rescor))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE.rescor, ROPE.rescor), linetype='dashed') +
    scale_fill_manual(values=c(wes_palette("Darjeeling1", n=5)[5], 'gray80')) +
    ylab('Residual Correlation (Causal Rating and Confidence)') +
    theme_classic() +
    theme(legend.position='none',
          panel.border=element_rect(color='black', fill=NA),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
ggsave('contrast-rescor.png')






fitRatings <- ratings %>%
    ggplot(aes(x=n, y=.value,
               group=normality, fill=normality)) +
    ylab('Estimated Mean Causal Rating') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure, labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

fitSD <- ratings %>%
    ggplot(aes(x=n, y=exp(sigma),
               group=normality, fill=normality)) +
    ylab('Estimated σ(Causal Rating)') + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    coord_cartesian(ylim=c(0.0, 0.5)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))

(fitRatings / fitSD) + plot_layout(guides='collect') +
    plot_annotation(tag_levels='A') &
    theme(legend.position='bottom')
ggsave('figure5.png', width=10, height=7.5)


fitConf <- confidence %>%
    ggplot(aes(x=n, y=.value,
               group=normality, fill=normality)) +
    ylab('Estimated Mean Confidence') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure, labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())

fitSDConf <- confidence %>%
    ggplot(aes(x=n, y=exp(sigma),
               group=normality, fill=normality)) +
    ylab('Estimated σ(Confidence)') + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    coord_cartesian(ylim=c(0.0, 0.5)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))

(fitConf / fitSDConf) + plot_layout(guides='collect') +
    plot_annotation(tag_levels='A') &
    theme(legend.position='bottom')
ggsave('figure6.png', width=10, height=7.5)


#########################################################################################
##
## Mean Causal Judgment:
##    Plot model fits
##
ratings %>%
    ggplot(aes(x=n, y=.value,
               group=normality, fill=normality)) +
    ylab('Estimated Mean Causal Rating') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("fit-ratings.png", width=6, height=4)

## Contrasts for mean rating (Abnormal - Normal)
ratings %>%
    compare_levels(.value, by=normality) %>% ungroup %>%
    mutate(normality=paste0('(', normality, ')')) %>%
    pivot_wider(names_from=c(structure, n, normality),
                names_sep=' ',
                values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-ROPE, ROPE)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
ratings %>%
    compare_levels(.value, by=normality) %>%
    ggplot(aes(x=n, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('Mean Causal Rating Contrasts: Abnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-ratings-normality.png', width=6, height=4)

## Plot model contrasts for mean rating (by N)
ratings %>%
    compare_levels(.value, by=n, comparison='ordered') %>%
    ungroup %>%
    pivot_wider(names_from=c(structure, normality, n), values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-ROPE, ROPE)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
ratings %>%
    compare_levels(.value, by=n, comparison='ordered') %>%
    ggplot(aes(x=n, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(breaks=c('1', '2'), values=c('gray80', PALETTE),
                      name='Normality', labels=c('Normal', 'Abnormal')) +
    ylab('Mean Causal Rating Contrasts: Number of Causes') +
    xlab('Number of Causes') +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-ratings-n.png', width=6, height=4)

#########################################################################################
##
## SD(Causal Judgment):
##    Plot model fits
##
ratings %>%
    ggplot(aes(x=n, y=exp(sigma),
               group=normality, fill=normality)) +
    ylab('Estimated σ(Causal Rating)') + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.75),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.75)) +
    coord_cartesian(ylim=c(0.0, 0.5)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ .,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("fit-sd-ratings-linear.png", width=6, height=4)

## Contrasts for SD ratings (Abnormal - Normal)
ratings %>%
    compare_levels(sigma, by=normality) %>% ungroup %>%
    mutate(normality=paste0('(', normality, ')')) %>%
    pivot_wider(names_from=c(structure, n, normality),
                names_sep=' ',
                values_from=sigma) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE, sdROPE)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
ratings %>%
    compare_levels(sigma, by=normality) %>%
    ggplot(aes(x=n, y=sigma, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('σ(Causal Rating) Contrasts: Abnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-sd-ratings-normality-linear.png', width=6, height=4)

## Contrasts for SD ratings (by N)
ratings %>%
    compare_levels(sigma, by=n, comparison='ordered') %>%
    ungroup %>%
    pivot_wider(names_from=c(structure, normality, n), values_from=sigma) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE, sdROPE)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
ratings %>%
    compare_levels(sigma, by=n, comparison='ordered') %>%
    ggplot(aes(x=n, y=sigma, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(breaks=c('1', '2'), values=c('gray80', PALETTE),
                      name='Normality', labels=c('Normal', 'Abnormal')) +
    ylab('σ(Causal Rating) Contrasts: Number of Causes') +
    xlab('Number of Causes') +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-sd-ratings-n-linear.png', width=6, height=4)





#########################################################################################
##
## Mean Confidence Rating:
##    Plot model fits
##
confidence %>%
    ggplot(aes(x=n, y=.value,
               group=normality, fill=normality)) +
    ylab('Estimated Mean Confidence') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.85),
             n=10000, geom='slab') +
    stat_pointinterval(position=position_dodge(0.85)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("fit-confidence.png", width=6, height=4)

## Contrasts for mean confidence (Abnormal - Normal)
confidence %>%
    compare_levels(.value, by=normality) %>% ungroup  %>%
    mutate(normality=paste0('(', normality, ')')) %>%
    pivot_wider(names_from=c(structure, n, normality),
                names_sep=' ',
                values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-ROPE.conf, ROPE.conf)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
confidence %>%
    compare_levels(.value, by=normality) %>%
    ggplot(aes(x=n, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('Mean Confidence Contrasts: Abnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-confidence-normality.png', width=6, height=4)

## Contrasts for mean confidence (by N)
confidence %>%
    compare_levels(.value, by=n, comparison='control') %>%
    ungroup %>%
    pivot_wider(names_from=n, values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-ROPE.conf, ROPE.conf)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
confidence %>%
    compare_levels(.value, by=n, comparison='ordered') %>%
    ggplot(aes(x=n, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(breaks=c('1', '2'), values=c('gray80', PALETTE),
                      name='Normality', labels=c('Normal', 'Abnormal')) +
    ylab('Mean Confidence Contrasts: Number of Causes') +
    xlab('Number of Causes') +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-confidence-n.png', width=6, height=4)



#########################################################################################
##
## SD(Confidence Ratings):
##    Plot model fits
##
confidence %>%
    ggplot(aes(x=n, y=exp(sigma),
               group=normality, fill=normality)) +
    ylab('Estimated σ(Confidence)') + xlab('Number of Causes') +
    stat_eye(mapping=aes(alpha=n), position=position_dodge(0.85),
             n=10000, geom='slab') + coord_cartesian(ylim=c(0.0, 0.5)) +
    stat_pointinterval(position=position_dodge(0.85)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave("fit-sd-confidence-linear.png", width=6, height=4)

## Contrasts for SD confidence (Abnormal - Normal)
confidence %>%
    compare_levels(sigma, by=normality) %>% ungroup %>%
    mutate(normality=paste0('(', normality, ')')) %>%
    pivot_wider(names_from=c(structure, n, normality),
                names_sep=' ',
                values_from=sigma) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE.conf, sdROPE.conf)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
confidence %>%
    compare_levels(sigma, by=normality) %>%
    ggplot(aes(x=n, y=sigma, group=normality,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('σ(Confidence) Contrasts: Abnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-sd-confidence-normality-linear.png', width=6, height=4)

## Contrasts for SD confidence (by N)
confidence %>%
    compare_levels(sigma, by=n, comparison='ordered') %>%
    ungroup %>%
    pivot_wider(names_from=c(structure, normality, n), values_from=sigma) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE.conf, sdROPE.conf)) %>%
    select(-c(CI, ROPE_CI, ROPE_low, ROPE_high))
confidence %>%
    compare_levels(sigma, by=n, comparison='ordered') %>%
    ggplot(aes(x=n, y=sigma, group=normality,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(breaks=c('1', '2'), values=c('gray80', PALETTE),
                      name='Normality', labels=c('Normal', 'Abnormal')) +
    ylab('σ(Confidence) Contrasts: Number of Causes') +
    xlab('Number of Causes') +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('contrast-sd-confidence-n-linear.png', width=6, height=4)


#########################################################################################
##
## Model Coefficients:
##    Plot model coefficients for mean causal rating, sd(causal rating),
##    mean confidence rating, and sd(confidence rating),
##
mMulti %>% gather_draws(b_rating_Intercept,
                        b_rating_structureD, b_rating_normalityA, `b_rating_structureD:normalityA`,
                        b_rating_n2, `b_rating_n2:structureD`, `b_rating_n2:normalityA`,
                        `b_rating_n2:structureD:normalityA`,
                        b_rating_n3, `b_rating_n3:structureD`, `b_rating_n3:normalityA`,
                        `b_rating_n3:structureD:normalityA`,
                        b_rating_n4, `b_rating_n4:structureD`, `b_rating_n4:normalityA`,
                        `b_rating_n4:structureD:normalityA`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_rating_'='', ':'=' : ',
                                        'structureD'='Causal Structure [disjunctive - conjunctive]',
                                        'normalityA'='Normality [abnormal - normal]',
                                        'n2'='Number of Causes [2 - 1]',
                                        'n3'='Number of Causes [3 - 1]',
                                        'n4'='Number of Causes [4 - 1]')),
                      levels=c('Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [4 - 1]','Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [3 - 1]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [2 - 1]',
                               'Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Normality [abnormal - normal]',
                               'Causal Structure [disjunctive - conjunctive]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > ROPE))) +
    xlab('Estimate') + ylab('') + stat_halfeyeh() +
    geom_vline(xintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                    subtitle='Mean Causal Judgment')
ggsave('coefficients_rating.png', width=10, height=8)

mMulti %>% gather_draws(b_sigma_rating_Intercept,
                        b_sigma_rating_structureD, b_sigma_rating_normalityA,
                        `b_sigma_rating_structureD:normalityA`,
                        b_sigma_rating_n2, `b_sigma_rating_n2:structureD`, `b_sigma_rating_n2:normalityA`,
                        `b_sigma_rating_n2:structureD:normalityA`,
                        b_sigma_rating_n3, `b_sigma_rating_n3:structureD`, `b_sigma_rating_n3:normalityA`,
                        `b_sigma_rating_n3:structureD:normalityA`,
                        b_sigma_rating_n4, `b_sigma_rating_n4:structureD`, `b_sigma_rating_n4:normalityA`,
                        `b_sigma_rating_n4:structureD:normalityA`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_sigma_rating_'='', ':'=' : ',
                                        'structureD'='Causal Structure [disjunctive - conjunctive]',
                                        'normalityA'='Normality [abnormal - normal]',
                                        'n2'='Number of Causes [2 - 1]',
                                        'n3'='Number of Causes [3 - 1]',
                                        'n4'='Number of Causes [4 - 1]')),
                      levels=c('Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [4 - 1]','Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [3 - 1]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [2 - 1]',
                               'Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Normality [abnormal - normal]',
                               'Causal Structure [disjunctive - conjunctive]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > sdROPE))) +
    xlab('Estimate') + ylab('') + stat_halfeyeh() +
    geom_vline(xintercept=c(-sdROPE, sdROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                    subtitle='Standard Deviation of Causal Judgments')
ggsave('coefficients_sd_rating.png', width=10, height=8)

mMulti %>% gather_draws(b_confidence_Intercept,
                        b_confidence_structureD, b_confidence_normalityA,
                        `b_confidence_structureD:normalityA`,
                        b_confidence_n2, `b_confidence_n2:structureD`, `b_confidence_n2:normalityA`,
                        `b_confidence_n2:structureD:normalityA`,
                        b_confidence_n3, `b_confidence_n3:structureD`, `b_confidence_n3:normalityA`,
                        `b_confidence_n3:structureD:normalityA`,
                        b_confidence_n4, `b_confidence_n4:structureD`, `b_confidence_n4:normalityA`,
                        `b_confidence_n4:structureD:normalityA`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_confidence_'='', ':'=' : ',
                                        'structureD'='Causal Structure [disjunctive - conjunctive]',
                                        'normalityA'='Normality [abnormal - normal]',
                                        'n2'='Number of Causes [2 - 1]',
                                        'n3'='Number of Causes [3 - 1]',
                                        'n4'='Number of Causes [4 - 1]')),
                      levels=c('Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [4 - 1]','Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [3 - 1]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [2 - 1]',
                               'Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Normality [abnormal - normal]',
                               'Causal Structure [disjunctive - conjunctive]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > ROPE.conf))) +
    xlab('Estimate') + ylab('') + stat_halfeyeh() +
    geom_vline(xintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Mean Confidence')
ggsave('coefficients_confidence.png', width=10, height=8)

mMulti %>% gather_draws(b_sigma_confidence_Intercept,
                        b_sigma_confidence_structureD, b_sigma_confidence_normalityA,
                        `b_sigma_confidence_structureD:normalityA`,
                        b_sigma_confidence_n2, `b_sigma_confidence_n2:structureD`,
                        `b_sigma_confidence_n2:normalityA`,
                        `b_sigma_confidence_n2:structureD:normalityA`,
                        b_sigma_confidence_n3, `b_sigma_confidence_n3:structureD`,
                        `b_sigma_confidence_n3:normalityA`,
                        `b_sigma_confidence_n3:structureD:normalityA`,
                        b_sigma_confidence_n4, `b_sigma_confidence_n4:structureD`,
                        `b_sigma_confidence_n4:normalityA`,
                        `b_sigma_confidence_n4:structureD:normalityA`) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('b_sigma_confidence_'='', ':'=' : ',
                                        'structureD'='Causal Structure [disjunctive - conjunctive]',
                                        'normalityA'='Normality [abnormal - normal]',
                                        'n2'='Number of Causes [2 - 1]',
                                        'n3'='Number of Causes [3 - 1]',
                                        'n4'='Number of Causes [4 - 1]')),
                      levels=c('Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [4 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [4 - 1]','Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [3 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [3 - 1]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Normality [abnormal - normal]',
                               'Number of Causes [2 - 1] : Causal Structure [disjunctive - conjunctive]',
                               'Number of Causes [2 - 1]',
                               'Causal Structure [disjunctive - conjunctive] : Normality [abnormal - normal]',
                               'Normality [abnormal - normal]',
                               'Causal Structure [disjunctive - conjunctive]',
                               'Intercept'))) %>%
    ggplot(aes(y=.variable, x=.value,
               fill=stat(abs(x) > sdROPE.conf))) +
    xlab('Estimate') + ylab('') + stat_halfeyeh() +
    geom_vline(xintercept=c(-sdROPE.conf, sdROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Standard Deviation of Confidence')
ggsave('coefficients_sd_confidence.png', width=10, height=8)


#########################################################################################
##
## Group-level effects:
##    Plot "random" effects
##
group_effects <- mMulti %>%
    gather_draws(sd_vignette__rating_Intercept,
                 sd_vignette__sigma_rating_Intercept,
                 sd_vignette__confidence_Intercept,
                 sd_vignette__sigma_confidence_Intercept,
                 cor_vignette__rating_Intercept__sigma_rating_Intercept,
                 cor_vignette__rating_Intercept__confidence_Intercept,
                 cor_vignette__sigma_rating_Intercept__confidence_Intercept,
                 cor_vignette__rating_Intercept__sigma_confidence_Intercept,
                 cor_vignette__sigma_rating_Intercept__sigma_confidence_Intercept,
                 cor_vignette__confidence_Intercept__sigma_confidence_Intercept) %>%
    ungroup() %>%
    mutate(.variable=
               factor(str_replace_all(.variable,
                                      c('sd_vignette__rating_Intercept'='σ(Intercept_rating)',
                                        'sd_vignette__sigma_rating_Intercept'='σ(Intercept_σ(rating))',
                                        'sd_vignette__confidence_Intercept'='σ(Intercept_confidence)',
                                        'sd_vignette__sigma_confidence_Intercept'=
                                            'σ(Intercept_σ(confidence))',
                                        'cor_vignette__rating_Intercept__sigma_rating_Intercept'=
                                            'correlation(Intercept_rating, Intercept_σ(rating))',
                                        'cor_vignette__confidence_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_confidence, Intercept_σ(confidence))',
                                        'cor_vignette__rating_Intercept__confidence_Intercept'=
                                            'correlation(Intercept_rating, Intercept_confidence)',
                                        'cor_vignette__sigma_rating_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_σ(rating), Intercept_σ(confidence))',
                                        'cor_vignette__sigma_rating_Intercept__confidence_Intercept'=
                                            'correlation(Intercept_confidence, Intercept_σ(rating))',
                                        'cor_vignette__rating_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_rating, Intercept_σ(confidence))')),
                      levels=c('correlation(Intercept_rating, Intercept_σ(confidence))',
                               'correlation(Intercept_confidence, Intercept_σ(rating))',
                               'correlation(Intercept_σ(rating), Intercept_σ(confidence))',
                               'correlation(Intercept_rating, Intercept_confidence)',
                               'correlation(Intercept_confidence, Intercept_σ(confidence))',
                               'correlation(Intercept_rating, Intercept_σ(rating))',
                               'σ(Intercept_σ(confidence))','σ(Intercept_confidence)',
                               'σ(Intercept_σ(rating))', 'σ(Intercept_rating)')))

group_effects %>%
    filter(!str_detect(.variable, 'correlation')) %>%
    pivot_wider(names_from=.variable, values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(., ci=0.95, test=c()) %>%
    select(-c(CI))

group_effects %>%
    filter(str_detect(.variable, 'correlation')) %>%
    pivot_wider(names_from=.variable, values_from=.value) %>%
    select(-.chain, -.iteration, -.draw) %>%
    describe_posterior(., ci=0.95, rope_ci=0.95,
                       rope_range=c(-sd(unlist(.))*0.1, sd(unlist(.))*0.1)) %>%
    select(-c(CI))

ge1 <- group_effects %>%
    filter(!str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeyeh(normalize='xy', fill=wes_palette("Darjeeling1", n=5)[5]) +
    theme_classic() + theme(legend.position='none') +
    coord_cartesian(xlim=c(0, 0.4))

ge2 <- group_effects %>%
    filter(str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeyeh(fill=wes_palette("Darjeeling1", n=5)[5]) +
    theme_classic() + theme(legend.position='none')


ge1 / ge2 + plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Group-Level Effects',
                            tag_levels = 'A')
ggsave('coefficients_re.png', width=8, height=6)
