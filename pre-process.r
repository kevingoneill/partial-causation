#!/usr/bin/Rscript
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
    writeLines('Usage: ./pre-process.r <input-file-name> <output-file-name>')
    quit()
}
in_file <- args[1]
out_file <- args[2]

data_wide <- read.csv(in_file, header=TRUE, stringsAsFactors=FALSE)

## remove unneeded rows/cols
data_wide <- data_wide[-c(1, 2, 3),]
data_wide <- data_wide[, -c(1:5, 7:20, 54, 56, 58)]

## group the vignettes into two columns
df <- gather(data_wide, question, response, Bridge_1_Resp_1:Moving_Many_Conf_1)
colnames(df) <- c('duration', 'age', 'sex', 'attn_check', 'n', 'nreq', 'id', 'question', 'response')
df <- df[order(df$n, df$nreq, df$id),]
df <- df[-which(df$response == ''),]   # filter out missing responses

## split question identifier into conditions
df <- separate(df, 'question', c('vignette', 'ncond', 'measure'),
               sep='_', extra='drop')
df$measure <- tolower(df$measure)

## convert back to wide format- one row per vignette, columns for resp, conf
df <- df %>% spread(measure, response) %>% rename("rating"=resp, "confidence"=conf)
df$condition <- ifelse(df$nreq == df$n, 'JC', 'OD')
df$condition[1:99] <- 'OD'

df <- df[, c('id', 'attn_check', 'age', 'sex', 'duration',
             'vignette', 'n', 'nreq', 'condition', 'rating', 'confidence')]

## save demographic info
#write.csv(subset(data_wide, select=c('Age', 'Sex', 'AttnCheck')),
#          'demographics.csv', row.names=FALSE)

## filter out by attention check
df <- subset(subset(df, attn_check=='Yes.'), select=-c(attn_check))

write.csv(df, out_file, row.names=FALSE)
