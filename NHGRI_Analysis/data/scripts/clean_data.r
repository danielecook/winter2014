#! /usr/bin/env Rscript

library(stringr)
library(reshape2)
library("annotate")

DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  user <- Sys.info()[["user"]]
  setwd(sprintf("/Users/%s/Documents/git/winter2014/NHGRI_Analysis/data/",user))
} else {
  setwd(DIR)
}
df = data.frame(read.table('gwascatalog.txt',header=TRUE,sep="\t",quote="",comment.char="",stringsAsFactors=TRUE))
GO = data.frame(read.table('GO/GO_reshaped.txt',header=T,sep='\t',stringsAsFactors=T))

# Date Conversions (ISO Standard)
df$Date.Added.to.Catalog <- as.Date(df$Date.Added.to.Catalog,"%m/%d/%Y")
df$Date <- as.Date(df$Date,"%m/%d/%Y")
df$Journal <- factor(df$Journal)

# Ensure that DOID is installed.
source("http://bioconductor.org/biocLite.R")
biocLite("DO.db")
library("DO.db")


pmids <- df$PUBMEDID[1:20]

# Clean up Population Data
extract_study_population <- function(col) {
  x <- df[,col]
  m<- sapply(x,function(i) str_match_all(i,"^([0-9,]{3,7}) (.*?) ancestry (.*?), ([0-9,]{3,7}) (.*?) ancestry (.*)")[[1]][2:7])
  t(m) # AH - you have to transpose these!
}

# Extract case/control information.
df[ c("nCases","case_ancestry","is_case","nControl","control_ancestry","is_control")] <-extract_study_population("Initial.Sample.Size")
