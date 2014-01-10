#! /usr/bin/env Rscript

DIR = commandArgs(trailingOnly = TRUE)[1]
if (is.na(DIR[1])) {
  setwd("/Users/daniel/Documents/git/winter2014/NHGRI_Analysis/data/scripts/")
} else {
  setwd(commandArgs(trailingOnly = TRUE)[1])
}
gwascatalog = data.frame(read.table('../gwascatalog.txt',header=TRUE,sep="\t",quote="",comment.char="",as.is=TRUE))

