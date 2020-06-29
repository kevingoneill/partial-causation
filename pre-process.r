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

## remove unneeded cols
data_wide <- data_wide[, -c(1:5, 7:20)]

## group the vignettes into two columns
df <- pivot_longer(data_wide, E_1_C_N_Rating_1:G_4_D_A_Confidence_1,
                   names_to='question', values_to='response')
colnames(df) <- c('duration', 'age', 'sex', 'attn_check',
                  'id', 'question', 'response')
df <- df[order(df$id),]
df <- df[!is.na(df$response),]   # filter out missing responses

## split question identifier into conditions
df <- separate(df, 'question', c('vignette', 'n', 'structure',
                                 'normality', 'measure'),
               sep='_', extra='drop') %>%
    mutate(measure=tolower(measure),
           id=as.numeric(id),
           duration=as.numeric(duration),
           age=as.numeric(age),
           n=as.numeric(n),
           response=as.numeric(response))

## convert back to wide format- one row per vignette, columns for resp, conf
df <- df %>% pivot_wider(names_from=measure, values_from=response)

df <- df[, c('id', 'attn_check', 'age', 'sex', 'duration',
             'vignette', 'n', 'structure', 'normality', 'rating', 'confidence')]

## filter out by attention check
writeLines(sprintf('Excluding %d participants.', sum(df$attn_check != 'Yes.')))
df <- df %>% filter(attn_check=='Yes.') %>%
    select(-attn_check)

## save output file
write.csv(df, out_file, row.names=FALSE)
