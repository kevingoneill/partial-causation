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
data_wide <- data_wide[, -c(1:5, 7:20)]

## group the vignettes into two columns
df <- pivot_longer(data_wide, E_1_C_N_Rating_1:G_4_D_A_Confidence_1,
                   names_to='question', values_to='response')
colnames(df) <- c('duration', 'age', 'sex', 'attn_check',
                  'id', 'question', 'response')
df <- df[order(df$id),]
df <- df[-which(df$response == ''),]   # filter out missing responses

## split question identifier into conditions
df <- separate(df, 'question', c('vignette', 'n', 'structure',
                                 'normality', 'measure'),
               sep='_', extra='drop')
df$measure <- tolower(df$measure)

## convert back to wide format- one row per vignette, columns for resp, conf
df <- df %>% pivot_wider(names_from=measure, values_from=response)

df <- df[, c('id', 'attn_check', 'age', 'sex', 'duration',
             'vignette', 'n', 'structure', 'normality', 'rating', 'confidence')]

## save demographic info
#write.csv(subset(data_wide, select=c('Age', 'Sex', 'AttnCheck')),
#          'demographics.csv', row.names=FALSE)

## filter out by attention check
df <- subset(subset(df, attn_check=='Yes.'), select=-c(attn_check))
write.csv(df, out_file, row.names=FALSE)
