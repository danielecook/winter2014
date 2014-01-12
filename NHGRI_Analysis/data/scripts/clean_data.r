#! /usr/bin/env Rscript

DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  setwd("/Users/Dan/Documents/git/winter2014/NHGRI_Analysis/data/")
} else {
  setwd(DIR)
}
df = data.frame(read.table('gwascatalog.txt',header=TRUE,sep="\t",quote="",comment.char="",as.is=TRUE))

# Date Conversions (ISO Standard)
df$Date.Added.to.Catalog <- as.Date(df$Date.Added.to.Catalog,"%m/%d/%Y")
df$Date <- as.Date(df$Date,"%m/%d/%Y")
df$Journal <- factor(df$Journal)
