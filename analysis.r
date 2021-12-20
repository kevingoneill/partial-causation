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
library(ggisoband)   ## devtools::install_github("clauswilke/ggisoband")
library(patchwork)
library(tibble)
library(ggh4x)      ## devtools::install_github("teunbrand/ggh4x")

## Set a color palette for plotting
PALETTE <- rev(wes_palette("Darjeeling1", n=2))

## Read in and normalize data
judgments <- read.csv('data/processed_data.csv', header=TRUE) %>%
    mutate(rating=rating / 100,
           confidence=confidence / 100,
           structure=factor(structure),
           n=as.factor(n),
           normality=factor(normality, levels=c('N', 'A')),
           norm_type=ifelse(vignette %in% c('St', 'C', 'G'),
                            'descriptive', 'prescriptive'))

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
              prior=c(set_prior('normal(0, 1.0)', resp='rating'),
                      set_prior('normal(0, 5.0)', resp='rating', dpar='sigma'),
                      set_prior('normal(0, 1.0)', resp='confidence'),
                      set_prior('normal(0, 5.0)', resp='confidence', dpar='sigma')),
              data=judgments, file='mMulti', iter=5000,
              inits="0", sample_prior="yes", save_pars=save_pars(all=TRUE),
              cores=4, control=list(adapt_delta=0.99))
summary(mMulti, priors=TRUE)




######################################################################################
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

## null model to compute full correlation b/t ratings and confidence
mNull <- brm(bf(rating ~ 1 + (1 |v| vignette),
                sigma ~ 1 + (1 |v| vignette)) +
             bf(confidence ~ 1 + (1 |v| vignette),
                sigma ~ 1 + (1 |v| vignette)),
             data=judgments, file='mNull', iter=5000,
             inits="0", sample_prior="yes", save_all_pars=TRUE,
             cores=4, control=list(adapt_delta=0.99))

LOO(mMulti, mReduced, mZOIB)
model_weights(mMulti, mReduced, mZOIB)


## Plot raw causal ratings
judgments %>%
    ggplot(aes(x=n, y=rating, group=normality, fill=normality, side=ifelse(normality=='N', 'left', 'right'))) +
    ylab('Causal Rating') + ylim(0, 1) + xlab('Number of Causes') +
    stat_halfeye(mapping=aes(alpha=n), position=position_dodge(0.5), adjust=0.75,
                 n=nrow(judgments), normalize='panels', geom='slab') +
    stat_pointinterval(position=position_dodge(0.5)) +
    scale_alpha_discrete(range = c(0.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/data-ratings.png', width=6, height=4)

## Plot raw confidence ratings
judgments %>%
    ggplot(aes(x=n, y=confidence, group=normality, fill=normality, side=ifelse(normality=='N', 'left', 'right'))) +
    ylab('Confidence') + ylim(0, 1) + xlab('Number of Causes') +
    stat_halfeye(mapping=aes(alpha=n), position=position_dodge(0.5), adjust=0.75,
             n=nrow(judgments), normalize='panels', geom='slab') +
    stat_pointinterval(position=position_dodge(0.5)) +
    scale_alpha_discrete(range = c(0.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/data-confidence.png', width=6, height=4)

## Plot raw causal and confidence ratings
judgments.mean <- judgments %>%
    group_by(normality, structure, n) %>%
    summarize(.rating=mean(rating), .rating.lower=mean_cl_normal(rating)$ymin,
              .rating.upper=mean_cl_normal(rating)$ymax, .confidence=mean(confidence),
              .confidence.lower=mean_cl_normal(confidence)$ymin, .confidence.upper=mean_cl_normal(confidence)$ymax)

plot.bivariate <- judgments %>%
    ggplot(aes(x=confidence, y=rating,
               color=normality, fill=normality)) +
    geom_density_bands(aes(alpha=stat(ndensity)), bins=10,
                       h=c(0.4, 0.4), show.legend=c(alpha=FALSE)) +
    geom_point(aes(x=.confidence, y=.rating), color='black', data=judgments.mean) +
    geom_errorbar(aes(x=.confidence, y=.rating, ymin=.rating.lower, ymax=.rating.upper), color='black', data=judgments.mean) +
    geom_errorbarh(aes(x=.confidence, xmin=.confidence.lower, xmax=.confidence.upper, y=.rating), color='black', data=judgments.mean) +
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
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1),
                       labels=c('', '.25', '', '.75', ''),
                       expand=expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
                       labels=c('', '.25', '', '.75', ''),
                       expand=expansion(mult = c(0, 0))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          panel.spacing=unit(0, 'lines'),
          legend.position='none',
          axis.text=element_text(size=12),
          axis.title=element_text(size=24),
          strip.text=element_text(size=12))
plot.bivariate

ggsave('plots/data-bivariate.png', plot.bivariate, width=5, height=5)


## extract model predictions
constrain <- function (x) {pmax(0, pmin(1, x))}
predictions <- full_join(predicted_draws(mMulti, judgments, resp='rating',
                                         value='.prediction.rating'),
                         predicted_draws(mMulti, judgments, resp='confidence',
                                         value='.prediction.confidence'))
draws <- full_join(epred_draws(mMulti, data_grid(judgments, normality, structure, n),
                                resp='rating', value='.rating', re_formula=NA),
                   epred_draws(mMulti, data_grid(judgments, normality, structure, n),
                                resp='confidence', value='.confidence', re_formula=NA)) %>%
    median_hdi


## Plot model predictions for causal ratings
predictions %>%
    ggplot(aes(x=n, y=constrain(.prediction.rating),
               group=normality, fill=normality, side=ifelse(normality=='N', 'left', 'right'))) +
    ylab('Causal Rating') + coord_cartesian(ylim=c(0, 1)) +
    xlab('Number of Causes') +
    stat_halfeye(mapping=aes(alpha=n), position=position_dodge(0.5),
                 normalize='panels', geom='slab') +
    stat_pointinterval(aes(y=.prediction.rating), position=position_dodge(0.5)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ .,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))

ggsave('plots/predictions-ratings.png', width=6, height=4)

## Plot model predictions for confidence
predictions %>%
    ggplot(aes(x=n, y=constrain(.prediction.confidence),
               group=normality, fill=normality, side=ifelse(normality=='N', 'left', 'right'))) +
    ylab('Causal Rating') + coord_cartesian(ylim=c(0, 1)) +
    xlab('Number of Causes') +
    stat_halfeye(mapping=aes(alpha=n), position=position_dodge(0.5),
             normalize='panels', geom='slab') +
    stat_pointinterval(aes(y=.prediction.confidence), position=position_dodge(0.5)) +
    scale_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid(structure ~ .,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() + theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/predictions-confidence.png', width=6, height=4)

## Plot model predictions for causal ratings and confidence
pred.bivariate <- predictions %>%
    ##filter(.draw == sample(.draw, length(.draw) * 0.75)) %>%
    ggplot(aes(x=constrain(.prediction.confidence),
               y=constrain(.prediction.rating),
               color=normality, fill=normality)) +
    geom_density_bands(aes(alpha=stat(ndensity)), bins=10,
                       h=c(0.25, 0.25), show.legend=FALSE) +
    geom_point(aes(x=.confidence, y=.rating), color='black', data=draws) +
    geom_errorbar(aes(x=.confidence, y=.rating, ymin=.rating.lower, ymax=.rating.upper), color='black', data=draws) +
    geom_errorbarh(aes(x=.confidence, xmin=.confidence.lower, xmax=.confidence.upper, y=.rating), color='black', data=draws) +
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
    scale_x_continuous(breaks=c(0, .25, .5, .75, 1),
                       labels=c('', '.25', '', '.75', ''),
                       expand=expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks=c(0, .25, .5, .75, 1),
                       labels=c('', '.25', '', '.75', ''),
                       expand=expansion(mult = c(0, 0))) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA),
          panel.spacing=unit(0, 'lines'),
          legend.position='none',
          axis.text=element_text(size=12),
          axis.title=element_text(size=24),
          strip.text=element_text(size=12))
pred.bivariate
ggsave('plots/predictions-bivariate.png', pred.bivariate, width=5, height=5)


plot.bivariate + pred.bivariate +
    plot_annotation(tag_levels='A') +
    plot_layout(guides='collect') &
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=24),
          strip.text=element_text(size=12))
ggsave('plots/figure3.png', width=10, height=5)




## Define ROPE ranges as 0.1 * sd of the DV
## NOTE: ROPEs for SD are estimated using posteriors over the data
ROPE = sd(judgments$rating) * 0.1
sdROPE <- sd(add_linpred_draws(judgments, mMulti, resp='rating', dpar='sigma')$sigma) * 0.1
ROPE.conf <- sd(judgments$confidence) * 0.1
sdROPE.conf <- sd(add_linpred_draws(judgments, mMulti, resp='confidence', dpar='sigma')$sigma) * 0.1

## Gather marginal means of all parameters
em.cause <- emmeans(mMulti, ~ normality*n*structure, resp='rating')
em.cause.sd <- emmeans(mMulti, ~ normality*n*structure, resp='rating', dpar='sigma')
em.conf <- emmeans(mMulti, ~ normality*n*structure, resp='confidence')
em.conf.sd <- emmeans(mMulti, ~ normality*n*structure, resp='confidence', dpar='sigma')







fitRatings <- gather_emmeans_draws(em.cause) %>%
    ggplot(aes(x=n, y=.value, side=ifelse(normality=='N', 'left', 'right'),
               group=normality, fill=normality)) +
    ylab('Estimated Mean\nCausal Judgment') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_halfeye(mapping=aes(slab_alpha=n), position=position_dodge(0.5), n=10000) +
    scale_slab_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure, labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16),
          axis.text.y=element_text(size=12),
          strip.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12))

fitSD <- gather_emmeans_draws(em.cause.sd, 'sigma') %>%
    ggplot(aes(x=n, y=exp(sigma), side=ifelse(normality=='N', 'left', 'right'),
               group=normality, fill=normality)) +
    ylab('Estimated\nStandard Deviation\nof Causal Judgments') +
    xlab('Number of Causes') +
    stat_halfeye(mapping=aes(slab_alpha=n), position=position_dodge(0.5), n=10000) +
    coord_cartesian(ylim=c(0.0, 0.5)) +
    scale_slab_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=16),
          axis.text=element_text(size=12),
          strip.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12))

(fitRatings / fitSD) + plot_layout(guides='collect') +
    plot_annotation(tag_levels='A') &
    theme(panel.border=element_rect(color='black', fill=NA),
          legend.position='bottom')
ggsave('plots/fit-ratings.png', width=9, height=6)


fitConf <- gather_emmeans_draws(em.conf) %>%
    ggplot(aes(x=n, y=.value, side=ifelse(normality=='N', 'left', 'right'),
               group=normality, fill=normality)) +
    ylab('Estimated\nMean Confidence') + coord_cartesian(ylim=c(.5, 1)) +
    xlab('Number of Causes') +
    stat_halfeye(mapping=aes(slab_alpha=n), position=position_dodge(0.5), n=10000) +
    scale_slab_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality', labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure, labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16),
          axis.text.y=element_text(size=12),
          strip.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12))

fitSDConf <- gather_emmeans_draws(em.conf.sd, 'sigma') %>%
    ggplot(aes(x=n, y=exp(sigma), side=ifelse(normality=='N', 'left', 'right'),
               group=normality, fill=normality)) +
    ylab('Estimated\nStandard Deviation\nof Confidence') + xlab('Number of Causes') +
    stat_halfeye(mapping=aes(slab_alpha=n), position=position_dodge(0.5), n=10000) +
    coord_cartesian(ylim=c(0.0, 0.5)) +
    scale_slab_alpha_discrete(range = c(.5, 1), guide='none') +
    scale_fill_manual(values=PALETTE, name='Normality',
                      labels=c('Normal', 'Abnormal')) +
    facet_grid( ~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    theme_classic() +
    theme(axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=16),
          axis.text=element_text(size=12),
          strip.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12))

(fitConf / fitSDConf) + plot_layout(guides='collect') +
    plot_annotation(tag_levels='A') &
    theme(panel.border=element_rect(color='black', fill=NA),
          legend.position='bottom')
ggsave('plots/fit-confidence.png', width=9, height=6)


######################################################################################
##
## Mean Causal Judgment:
##
## Contrasts for mean rating (Abnormal - Normal)
emnorm.cause <- contrast(em.cause, 'trt.vs.ctrl', simple='normality')
cnorm.cause <- describe_posterior(emnorm.cause, ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(p=pd_to_p(pd))
cnorm.cause

emnorm.cause %>%
    gather_emmeans_draws() %>%
    ggplot(aes(x=n, y=.value, group=contrast,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('Mean Causal Judgment Contrasts:\nAbnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-ratings-normality.png', width=6, height=4)



## Plot model contrasts for mean rating (by 4 causes - 1 cause)
emn.cause <- emmeans(mMulti, ~ normality*n*structure,
                     resp='rating', at=list(n=c('1', '4'))) %>%
    contrast('trt.vs.ctrl', interaction=TRUE,
             simple=list('n', c('n', 'normality')), combine=TRUE)
cn.cause <- emn.cause %>%
    describe_posterior(ci=0.95, rope_ci=0.95, rope_range=c(-ROPE, ROPE)) %>%
    mutate(Parameter=c('N, C, 4 - 1', 'A, C, 4 - 1',
                       'N, D, 4 - 1', 'A, D, 4 - 1',
                       'A - N, C, 4 - 1', 'A - N, D, 4 - 1'),
           p=pd_to_p(pd))
cn.cause

emn.cause %>%
    gather_emmeans_draws() %>%
    unite(normality, normality, normality_trt.vs.ctrl, sep='') %>%
    mutate(normality=factor(str_remove_all(normality, '\\.'),
                            levels=c('N', 'A', 'A - N'))) %>%
    rename(n=n_trt.vs.ctrl) %>%
    ggplot(aes(x=normality, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), show.legend=FALSE) +
    facet_grid(~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE, wes_palette("Darjeeling1", n=5)[5])) +
    ylab('Mean Causal Judgment Contrasts:\n4 - 1 Causes') +
    scale_x_discrete(name='Normality', labels=c('Normal', 'Abnormal', 'Normal -\nAbnormal')) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-ratings-n.png', width=6, height=3)






######################################################################################
##
## SD(Causal Judgment):
##
## Contrasts for SD ratings (Abnormal - Normal)
emnorm.cause.sd <- contrast(em.cause.sd, 'trt.vs.ctrl', simple='normality')
cnorm.cause.sd <- describe_posterior(emnorm.cause.sd, ci=0.95, rope_ci=0.95,
                                     rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(p=pd_to_p(pd))
cnorm.cause.sd

emnorm.cause.sd %>%
    gather_emmeans_draws(value='sigma') %>%
    ggplot(aes(x=n, y=sigma, group=contrast,
               fill=stat(ifelse(abs(y) < sdROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-sdROPE, sdROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('SD of Causal Judgments Contrasts:\nAbnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-sd-ratings-normality.png', width=6, height=4)



## Contrasts for SD ratings (by N)
emn.cause.sd <- emmeans(mMulti, ~ normality*n*structure,
                        resp='rating', dpar='sigma',
                        at=list(n=c('1', '4'))) %>%
    contrast('trt.vs.ctrl', interaction=TRUE,
             simple=list('n', c('n', 'normality')), combine=TRUE)
cn.cause.sd <- describe_posterior(emn.cause.sd, ci=0.95, rope_ci=0.95,
                                  rope_range=c(-sdROPE, sdROPE)) %>%
    mutate(Parameter=c('N, C, 4 - 1', 'A, C, 4 - 1',
                       'N, D, 4 - 1', 'A, D, 4 - 1',
                       'A - N, C, 4 - 1', 'A - N, D, 4 - 1'),
           p=pd_to_p(pd))
cn.cause.sd


emn.cause.sd %>%
    gather_emmeans_draws() %>%
    unite(normality, normality, normality_trt.vs.ctrl, sep='') %>%
    mutate(normality=factor(str_remove_all(normality, '\\.'),
                            levels=c('N', 'A', 'A - N'))) %>%
    rename(n=n_trt.vs.ctrl) %>%
    ggplot(aes(x=normality, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE, '0', group)))) +
    stat_eye(position=position_dodge(width=1), show.legend=FALSE) +
    facet_grid(~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE, wes_palette("Darjeeling1", n=5)[5])) +
    ylab('SD of Causal Judgments Contrasts:\n4 - 1 Causes') +
    scale_x_discrete(name='Normality', labels=c('Normal', 'Abnormal', 'Normal -\nAbnormal')) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-sd-ratings-n.png', width=6, height=3)




######################################################################################
##
## Mean Confidence Rating:
##
## Contrasts for mean confidence (Abnormal - Normal)
emnorm.conf <- contrast(em.conf, 'trt.vs.ctrl', simple='normality')
cnorm.conf <- describe_posterior(emnorm.conf, ci=0.95, rope_ci=0.95,
                                 rope_range=c(-ROPE.conf, ROPE.conf)) %>%
    mutate(p=pd_to_p(pd))
cnorm.conf

emnorm.conf %>% gather_emmeans_draws %>%
    ggplot(aes(x=n, y=.value, group=contrast,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('Mean Confidence Contrasts:\nAbnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-confidence-normality.png', width=6, height=4)


## Contrasts for mean confidence (by N)
emn.conf <- emmeans(mMulti, ~ normality*n*structure,
                     resp='confidence', at=list(n=c('1', '4'))) %>%
    contrast('trt.vs.ctrl', interaction=TRUE,
             simple=list('n', c('n', 'normality')), combine=TRUE)
cn.conf <- emn.conf %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-ROPE.conf, ROPE.conf)) %>%
    mutate(Parameter=c('N, C, 4 - 1', 'A, C, 4 - 1',
                       'N, D, 4 - 1', 'A, D, 4 - 1',
                       'A - N, C, 4 - 1', 'A - N, D, 4 - 1'),
           p=pd_to_p(pd))
cn.conf


emn.conf %>%
    gather_emmeans_draws() %>%
    unite(normality, normality, normality_trt.vs.ctrl, sep='') %>%
    mutate(normality=factor(str_remove_all(normality, '\\.'),
                            levels=c('N', 'A', 'A - N'))) %>%
    rename(n=n_trt.vs.ctrl) %>%
    ggplot(aes(x=normality, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), show.legend=FALSE) +
    facet_grid(~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE, wes_palette("Darjeeling1", n=5)[5])) +
    ylab('Mean Confidence Contrasts:\n4 - 1 Causes') +
    scale_x_discrete(name='Normality', labels=c('Normal', 'Abnormal', 'Normal -\nAbnormal')) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-confidence-n.png', width=6, height=3)



######################################################################################
##
## SD(Confidence Ratings):
##
## Contrasts for SD confidence (Abnormal - Normal)
emnorm.conf.sd <- contrast(em.conf.sd, 'trt.vs.ctrl', simple='normality')
cnorm.conf.sd <- describe_posterior(emnorm.conf.sd, ci=0.95, rope_ci=0.95,
                                    rope_range=c(-sdROPE.conf, sdROPE.conf)) %>%
    mutate(p=pd_to_p(pd))
cnorm.conf.sd

emnorm.conf.sd %>%
    gather_emmeans_draws('sigma') %>%
    ggplot(aes(x=n, y=sigma, group=contrast,
               fill=stat(ifelse(abs(y) < ROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), n=10000) +
    geom_hline(yintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    facet_grid(structure ~ ., labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    ylab('SD of Confidence Contrasts:\nAbnormal - Normal') +
    xlab('Number of Causes') +
    theme_classic() + theme(legend.position='none',
                            panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-sd-confidence-normality.png', width=6, height=4)

## Contrasts for SD confidence (by N)
emn.conf.sd <- emmeans(mMulti, ~ normality*n*structure,
                     resp='confidence', dpar='sigma', at=list(n=c('1', '4'))) %>%
    contrast('trt.vs.ctrl', interaction=TRUE,
             simple=list('n', c('n', 'normality')), combine=TRUE)
cn.conf.sd <- emn.conf.sd %>%
    describe_posterior(ci=0.95, rope_ci=0.95,
                       rope_range=c(-sdROPE.conf, sdROPE.conf)) %>%
    mutate(Parameter=c('N, C, 4 - 1', 'A, C, 4 - 1',
                       'N, D, 4 - 1', 'A, D, 4 - 1',
                       'A - N, C, 4 - 1', 'A - N, D, 4 - 1'),
           p=pd_to_p(pd))
cn.conf.sd

emn.conf.sd %>%
    gather_emmeans_draws() %>%
    unite(normality, normality, normality_trt.vs.ctrl, sep='') %>%
    mutate(normality=factor(str_remove_all(normality, '\\.'),
                            levels=c('N', 'A', 'A - N'))) %>%
    rename(n=n_trt.vs.ctrl) %>%
    ggplot(aes(x=normality, y=.value, group=normality,
               fill=stat(ifelse(abs(y) < sdROPE.conf, '0', group)))) +
    stat_eye(position=position_dodge(width=1), show.legend=FALSE) +
    facet_grid(~ structure,
               labeller=labeller(structure=c(C='Conjunctive', D='Disjunctive'))) +
    geom_hline(yintercept=c(-sdROPE.conf, sdROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', PALETTE, wes_palette("Darjeeling1", n=5)[5])) +
    ylab('SD of Confidence Contrasts:\n4 - 1 Causes') +
    scale_x_discrete(name='Normality', labels=c('Normal', 'Abnormal', 'Normal -\nAbnormal')) +
    theme_classic() +
    theme(panel.border=element_rect(color='black', fill=NA))
ggsave('plots/contrast-sd-confidence-n.png', width=6, height=4)







######################################################################################
##
## Test for direct effect via model comparison
##
mDE <- brm(bf(confidence ~ n*structure*normality + (1 |v| vignette),
              sigma ~ n*structure*normality + (1 |v| vignette)) +
           bf(rating ~ n*structure*normality + confidence + (1 |v| vignette),
              sigma ~ n*structure*normality + confidence + (1 |v| vignette)) +
           set_rescor(FALSE),
           prior=c(set_prior('normal(0, 1.0)', resp='rating'),
                   set_prior('normal(0, 5.0)', resp='rating', dpar='sigma'),
                   set_prior('normal(0, 1.0)', resp='confidence'),
                   set_prior('normal(0, 5.0)', resp='confidence', dpar='sigma')),
           data=judgments, file='mDE', iter=5000,
           inits="0", sample_prior="yes", save_pars=save_pars(all=TRUE),
           cores=4, control=list(adapt_delta=0.99))

mDE.null <- brm(bf(confidence ~ n*structure*normality + (1 |v| vignette),
                   sigma ~ n*structure*normality + (1 |v| vignette)) +
                bf(rating ~ confidence + (1 |v| vignette),
                   sigma ~ confidence + (1 |v| vignette)) +
                set_rescor(FALSE),
                prior=c(set_prior('normal(0, 1.0)', resp='rating'),
                        set_prior('normal(0, 5.0)', resp='rating', dpar='sigma'),
                        set_prior('normal(0, 1.0)', resp='confidence'),
                        set_prior('normal(0, 5.0)', resp='confidence', dpar='sigma')),
                data=judgments, file='mDEnull', iter=5000,
                inits="0", sample_prior="yes", save_pars=save_pars(all=TRUE),
                cores=4, control=list(adapt_delta=0.99))

loo(mDE, mDE.null, moment_match=TRUE)
model_weights(mDE, mDE.null)






######################################################################################
##
## Contrast Tables

rbind(cnorm.cause, cnorm.cause.sd, cnorm.conf, cnorm.conf.sd) %>%
    separate(Parameter, into=c('Normality', 'N', 'Structure'), sep=', ') %>%
    add_column(Parameter=rep(c('Mean', 'SD'), 2, each=8), .before=1) %>%
    add_column(Variable=rep(c('Causal Judgment', 'Confidence'), each=16),
               .before=1) %>%
    as.data.frame %>%
    mutate(Median=round(Median, 2), pd=round(pd, 2),
           ROPE_Percentage=round(ROPE_Percentage*100, 2),
           CI=sprintf('[%.2f, %.2f]', CI_low, CI_high),
           ROPE=sprintf('[%.2f, %.2f]', ROPE_low, ROPE_high)) %>%
    select(Parameter, Variable, Structure, N,
           Median, CI, pd, p, ROPE, ROPE_Percentage) %>%
    arrange(Structure, Variable, Parameter, N) %>%
    write.csv('contrasts_normality.csv', row.names=FALSE)



rbind(cn.cause, cn.cause.sd, cn.conf, cn.conf.sd) %>%
    separate(Parameter, into=c('Normality', 'Structure', 'N'), sep=', ') %>%
    add_column(Parameter=rep(c('Mean', 'SD'), 2, each=6), .before=1) %>%
    add_column(Variable=rep(c('Causal Judgment', 'Confidence'), each=12),
               .before=1) %>%
    as.data.frame %>%
    mutate(Normality=factor(Normality, levels=c('A', 'N', 'A - N')),
           Median=round(Median, 2), pd=round(pd, 2),
           ROPE_Percentage=round(ROPE_Percentage*100, 2),
           CI=sprintf('[%.2f, %.2f]', CI_low, CI_high),
           ROPE=sprintf('[%.2f, %.2f]', ROPE_low, ROPE_high)) %>%
    select(Parameter, Variable, Structure, Normality,
           Median, CI, pd, p, ROPE, ROPE_Percentage) %>%
    arrange(Structure, Variable, Parameter, Normality) %>%
    write.csv('contrasts_n.csv', row.names=FALSE)


######################################################################################
##
## Vignette-level means:
##    Plot means of causal ratings and confidence for each vignette
##
ggplot(judgments) +
    aes(x=n, y=rating, group=normality, color=normality) +
    stat_summary(fun.data=mean_cl_normal, position=position_dodge(0.75)) +
    facet_grid(vignette ~ structure,
               labeller=labeller(vignette=c('C'='Cafe', 'E'='Electricity',
                                            'G'='Gym', 'R'='Rocket',
                                            'St'='Streaming', 'Se'='Sewage'),
                                 structure=c('C'='Conjunctive', 'D'='Disjunctive'))) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    ylab('Mean Causal Judgment') + coord_cartesian(ylim=c(.25, 1)) +
    xlab('Number of Causes') + theme_classic()
ggsave('plots/vignette-rating.png', height=10, width=7.5)

ggplot(judgments) +
    aes(x=n, y=confidence, group=normality, color=normality) +
    stat_summary(fun.data=mean_cl_normal, position=position_dodge(0.75)) +
    facet_grid(vignette ~ structure,
               labeller=labeller(vignette=c('C'='Cafe', 'E'='Electricity',
                                            'G'='Gym', 'R'='Rocket',
                                            'St'='Streaming', 'Se'='Sewage'),
                                 structure=c('C'='Conjunctive', 'D'='Disjunctive'))) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    ylab('Mean Confidence') + coord_cartesian(ylim=c(.25, 1)) +
    xlab('Number of Causes') + theme_classic()
ggsave('plots/vignette-confidence.png', height=10, width=7.5)


##
## Means by norm type (descriptive vs prescriptive)
##
ggplot(judgments) +
    aes(x=n, y=rating, group=normality, color=normality) +
    stat_summary(fun.data=mean_cl_normal, position=position_dodge(0.75)) +
    facet_grid(norm_type ~ structure,
               labeller=labeller(norm_type=c('descriptive'='Descriptive',
                                             'prescriptive'='Prescriptive'),
                                 structure=c('C'='Conjunctive', 'D'='Disjunctive'))) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    ylab('Mean Causal Judgment') + coord_cartesian(ylim=c(.25, 1)) +
    xlab('Number of Causes') + theme_classic()
ggsave('plots/norm-rating.png', height=5, width=7)

ggplot(judgments) +
    aes(x=n, y=confidence, group=normality, color=normality) +
    stat_summary(fun.data=mean_cl_normal, position=position_dodge(0.75)) +
    facet_grid(norm_type ~ structure,
               labeller=labeller(norm_type=c('descriptive'='Descriptive',
                                             'prescriptive'='Prescriptive'),
                                 structure=c('C'='Conjunctive', 'D'='Disjunctive'))) +
    scale_color_manual(values=PALETTE, name='Normality',
                       labels=c('Normal', 'Abnormal')) +
    ylab('Mean Confidence') + coord_cartesian(ylim=c(.25, 1)) +
    xlab('Number of Causes') + theme_classic()
ggsave('plots/norm-confidence.png', height=5, width=7)


######################################################################################
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
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-ROPE, ROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                    subtitle='Mean Causal Judgment')
ggsave('plots/coefficients_rating.png', width=10, height=8)

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
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-sdROPE, sdROPE), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                    subtitle='Standard Deviation of Causal Judgments')
ggsave('plots/coefficients_sd_rating.png', width=10, height=8)

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
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-ROPE.conf, ROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Mean Confidence')
ggsave('plots/coefficients_confidence.png', width=10, height=8)

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
    xlab('Estimate') + ylab('') + stat_halfeye() +
    geom_vline(xintercept=c(-sdROPE.conf, sdROPE.conf), linetype='dashed') +
    scale_fill_manual(values=c('gray80', wes_palette("Darjeeling1", n=5)[5])) +
    theme_classic() + theme(legend.position='none') +
    plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Standard Deviation of Confidence')
ggsave('plots/coefficients_sd_confidence.png', width=10, height=8)


######################################################################################
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
                                      c('sd_vignette__rating_Intercept'='σ(Intercept_cause)',
                                        'sd_vignette__sigma_rating_Intercept'='σ(Intercept_σ(cause))',
                                        'sd_vignette__confidence_Intercept'='σ(Intercept_confidence)',
                                        'sd_vignette__sigma_confidence_Intercept'=
                                            'σ(Intercept_σ(confidence))',
                                        'cor_vignette__rating_Intercept__sigma_rating_Intercept'=
                                            'correlation(Intercept_cause, Intercept_σ(cause))',
                                        'cor_vignette__confidence_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_confidence, Intercept_σ(confidence))',
                                        'cor_vignette__rating_Intercept__confidence_Intercept'=
                                            'correlation(Intercept_cause, Intercept_confidence)',
                                        'cor_vignette__sigma_rating_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_σ(cause), Intercept_σ(confidence))',
                                        'cor_vignette__sigma_rating_Intercept__confidence_Intercept'=
                                            'correlation(Intercept_confidence, Intercept_σ(cause))',
                                        'cor_vignette__rating_Intercept__sigma_confidence_Intercept'=
                                            'correlation(Intercept_cause, Intercept_σ(confidence))')),
                      levels=c('correlation(Intercept_cause, Intercept_σ(confidence))',
                               'correlation(Intercept_confidence, Intercept_σ(cause))',
                               'correlation(Intercept_σ(cause), Intercept_σ(confidence))',
                               'correlation(Intercept_cause, Intercept_confidence)',
                               'correlation(Intercept_confidence, Intercept_σ(confidence))',
                               'correlation(Intercept_cause, Intercept_σ(cause))',
                               'σ(Intercept_σ(confidence))','σ(Intercept_confidence)',
                               'σ(Intercept_σ(cause))', 'σ(Intercept_cause)')))

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
    describe_posterior(., ci=0.95, rope_ci=0.95) %>%
    select(-c(CI)) %>%
    mutate(p=pd_to_p(pd))

ge1 <- group_effects %>%
    filter(!str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeye(normalize='xy', fill=wes_palette("Darjeeling1", n=5)[5]) +
    theme_classic() + theme(legend.position='none') +
    coord_cartesian(xlim=c(0, 0.4))

ge2 <- group_effects %>%
    filter(str_detect(.variable, 'correlation')) %>%
    ggplot(aes(y=.variable, x=.value)) +
    xlab('Estimate') + ylab('') +
    stat_halfeye(fill=wes_palette("Darjeeling1", n=5)[5]) +
    theme_classic() + theme(legend.position='none')


ge1 / ge2 + plot_annotation(title='Coefficient Estimates for Experiment 2:',
                            subtitle='Group-Level Effects',
                            tag_levels = 'A')
ggsave('plots/coefficients_re.png', width=8, height=6)
