#! /usr/bin/env Rscript

library(stringr)

DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  user <- Sys.info()[["user"]]
  setwd(sprintf("/Users/%s/Documents/git/winter2014/NHGRI_Analysis/data/",user))
} else {
  setwd(DIR)
}
df = data.frame(read.table('gwascatalog.txt',header=TRUE,sep="\t",quote="",comment.char="",stringsAsFactors=TRUE))
GO = data.frame(read.table('GO/GO_reshaped.txt',header=T,sep='\t',stringsAsFactors=T))

# Preliminary GO Annotation Analysis:
df_GO <- merge(df,GO,by.x ="Snp_gene_ids", by.y ="Gene_ID")

# Get Go Usage Counts:
GO_terms <- grep("GO",colnames(df_GO))
f<-apply(df_GO[,GO_terms],2,function(x) sum(x!=""))
names(df_GO[order(f)][1:20])

# Population Extraction
head(df["Initial.Sample.Size"])

x <- df[,"Initial.Sample.Size"][1:200]
m <- gregexpr(c("^([0-9,]{3,7})"),df[,"Initial.Sample.Size"][1:200])
regmatches(x,m=m)

extract_study_population <- function(col) {
	x <- df[,col][1:500]
	m<- str_match_all(x,"^([0-9,]{3,7}) (.*?) ancestry (.*?), ([0-9,]{3,7}) (.*?) ancestry (.*)")
	do.call(rbind.data.frame,m)
}

extract_study_population("Initial.Sample.Size")


q<-do.call(rbind.data.frame,m)

