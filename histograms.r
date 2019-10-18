#!/usr/bin/Rscript
library(dplyr)
library(gplots)
library(RColorBrewer)

DIR <- 'histograms/'

args <- commandArgs(trailingOnly=T)
if (length(args) != 1) {
    writeLines('Usage: ./analisis.r <file-name>')
    quit()
}



plot_hist <- function(data) {
    rf <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
    COLORS <- rf(32)

    breaks = (0:101 - 0.5) / 100
    r <- hist(data$rating, breaks=breaks, plot=F)$counts
    c <- hist(data$confidence, breaks=breaks, plot=F)$counts
    
    par(mar=c(12,12,0,0), mgp=c(6,2,0))
    layout(matrix(c(2,0,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
    
    hist2d(data$rating, data$confidence, nbins=50, FUN=function(x) log(length(x)),
           xlim=c(0, 1.0), ylim=c(0, 1.0), xlab='Rating', ylab='Confidence',
           col=COLORS, cex.axis=3.0, cex.lab=4.0)
    par(mar=c(0,10,2,0))
    barplot(r, ylim=c(0, max(r)), xlim=c(0, 100), axes=F, space=0, col='red')
    par(mar=c(10,0,0,2))
    barplot(c, xlim=c(0, max(c)), ylim=c(0, 100), axes=F, space=0, col='red', horiz=T)
}


judgments <- read.csv(args[1], header=TRUE)
judgments$rating <- judgments$rating / 100
judgments$confidence <- judgments$confidence / 100
judgments$condition <- as.factor(judgments$condition)

time <- aggregate(duration ~ id, judgments, function(x) mean(x) / 60.0)
mean(time$duration)
sd(time$duration) / sqrt(nrow(time))

png(paste0(DIR, 'hist_all.png'), height=1000, width=1000)
plot_hist(judgments)
dev.off()

for (nc in unique(judgments$n)) {
    for (c in levels(judgments$condition)) {
        png(paste0(DIR, sprintf('hist_%d_%s.png', nc, c)), height=1000, width=1000)
        plot_hist(subset(judgments, n==nc & condition==c))
        dev.off()
    }
}

for (v in levels(judgments$vignette)) {
    png(paste0(DIR, v, '/hist_all.png'), height=1000, width=1000)
    plot_hist(subset(judgments, vignette==v))
    dev.off()
        
    for (nc in unique(judgments$n)) {
        for (c in levels(judgments$condition)) {
            png(paste0(DIR, v, sprintf('/hist_%d_%s.png', nc, c)), height=1000, width=1000)
            plot_hist(subset(judgments, n==nc & condition==c & vignette==v))
            dev.off()
        }
    }
}
