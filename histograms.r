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

    r <- hist(data$rating, breaks=101, plot=F)$counts
    c <- hist(data$confidence, breaks=101, plot=F)$counts

    par(mar=c(11,11,1,1), mgp=c(8,2.5,0))
    layout(matrix(c(2,0,1,3), 2, 2, byrow=T), c(8,2), c(2,8))

    hist2d(data$rating, data$confidence, nbins=50, FUN=function(x) log(length(x)),
           xlab='Rating', ylab='Confidence', col=COLORS, cex.axis=3.0, cex.lab=4.0)
    par(mar=c(0,8,2,0))
    barplot(r, axes=F, ylim=c(0, max(r)), space=0, col='red')
    par(mar=c(8,0,0,2))
    barplot(c, axes=F, xlim=c(0, max(c)), space=0, col='red', horiz=T)
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

nrow(judgments)
str(judgments)

for (nc in unique(judgments$n)) {
    for (c in levels(judgments$condition)) {
        print(sprintf('Condition: %d %s', nc, c))
        print(nrow(subset(judgments, n==nc & condition==c)))
        png(paste0(DIR, sprintf('hist_%d_%s.png', nc, c)), height=1000, width=1000)
        plot_hist(subset(judgments, n==nc & condition==c))
        dev.off()
    }
}
