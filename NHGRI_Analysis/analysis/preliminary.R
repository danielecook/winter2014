#! /usr/bin/env Rscript

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
