#! /usr/bin/env Rscript

library(stringr)
library(reshape2)
library("annotate")
library(RUnit)
library(data.table)
library(ggplot2)
library(source.gist)

source.gist('8809319')


DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  user <- Sys.info()[["user"]]
  setwd(sprintf("/Users/%s/Documents/git/winter2014/NHGRI_Analysis/data/",user))
} else {
  setwd(DIR)
}

# Load NHGRI Catalog
df = data.frame(read.table('gwascatalog.txt',header=TRUE,sep="\t",quote="",comment.char="",stringsAsFactors=F, strip.white=T))
# Create a column for merging the *strongest* SNP allele.
df$rs <- gsub("\\-.*$","",df$Strongest.SNP.Risk.Allele)


GO = data.frame(read.table('GO/GO_reshaped.txt',header=T,sep='\t',stringsAsFactors=T, strip.white=T))


#######################################
# HMAF = HapMap Allele Frequency Data #
#######################################

# Hapmap Allele Frequency Data
HMAF = data.frame(read.table('hapmap/hapmap_allele_freq_reshaped.csv',header=T,sep=",",stringsAsFactors=F, strip.white=T)) 

# Clean up the hapmap allele frequency data.
names(HMAF) <- gsub("\\.\\.","",gsub("\\.\\.\\.\\.","-",gsub("X..","",names(HMAF))))

condense_columns <- function(df,col_name, col_search) {
  df[[col_name]] <- NA
  for (p in colnames(df[,grep(col_search,names(df))])) {
    s <- is.na(df[[col_name]]) & (df[[p]] != "" & !is.na(df[[p]]))
    df[[col_name]][s] <-  as.character(df[[p]][s])
  }
  df <- df[,-grep(sub("-","\\\\-",col_search),colnames(df))]
  df
}

HMAF <- condense_columns(HMAF,"refallele","refallele-")
HMAF <- condense_columns(HMAF,"otherallele","otherallele-")



# Merge in HMAF data
df <- merge(df,HMAF,by=c("rs"), all.x=T,all.y=F)

# Lists SNPs that are in Hapmap but not matching with dataframe
# Now sure how/why these are here.
HMAF$rs[(!(HMAF$rs %in% df$rs))]

length(HMAF$rs[((HMAF$rs %in% df$rs))])



df$risk_allele <- gsub(" ","",do.call(rbind,strsplit(df$Strongest.SNP.Risk.Allele,"-"))[,2])
# Cleanup problems
df$risk_allele[df$risk_allele=="?"] <- NA
# Many difficult to interpret risk alleles remain (e.g. Cxrs12880735)
df$risk_allele[!df$risk_allele %in% c("A","C","T","G")]  <- NA

length(df$risk_allele[!is.na(df$risk_allele)]) # Number SNPs available.

# Number of HMAF SNPs that successfully merged


# Freqs of risk allele bases (tot = 11,168)
#    A    C    G    T 
# 3386 2357 2858 2567

# Number of Forward Strand SNPs
# Forward Strand variants
length(df[complete.cases(df[,c("risk_allele","refallele")]) & ((df$risk_allele == df$refallele) | df$risk_allele == df$otherallele),1])
# Reverse Strand variants
length(df[complete.cases(df[,c("risk_allele","refallele")]) & ((df$risk_allele != df$refallele) & df$risk_allele != df$otherallele),1])

# Create a chrom strand variable
df$chrom_strand <- 0

df <- corder(df,chrom_strand,risk_allele,refallele,otherallele)

# Flip Risk Allele if it appears to be on the reverse (-) strand:
df$chrom_strand[!is.na(df$risk_allele) & (df$risk_allele == df$refallele |  df$risk_allele == df$otherallele)] <- 1
df$chrom_strand[!is.na(df$risk_allele) & (df$risk_allele != df$refallele &  df$risk_allele != df$otherallele)] <- -1 # Reverse

length(df$chrom_strand[df$chrom_strand == -1]) # 282 variants apparently labeled on reverse strand.
comp_base = c(A="T",T="A",C="G",G="C")
# Create variable - risk_allele_corrected, to account for strand flip
# Reverse Strand
df$risk_allele_forward[df$chrom_strand == -1] <- as.character(comp_base[df[df$chrom_strand == -1,c("risk_allele")]])
# Forward Strand
df$risk_allele_forward[df$chrom_strand == 1] <- as.character(df[df$chrom_strand == 1,c("risk_allele")])

df <- corder(df,risk_allele_forward)



# Generate cumulative risk allele freq (for plotting purposes).
df$risk_allele_count <- ifelse(df$risk_allele == df$refallele, rowSums(df[,grep("refallele_count", names(df))],na.rm=T), rowSums(df[,grep("otherallele_count", names(df))],na.rm=T))
df$total_allele_count <- rowSums(df[,grep("totalcount",names(df))],na.rm=T)
df$risk_allele_freq <- df$risk_allele_count / df$total_allele_count

df <- corder(df, risk_allele_count,total_allele_count,risk_allele_freq)
# Generate Risk allele frequency
pops <- gsub('refallele_freq-','',names(df[,grep("refallele_freq",names(df))])) # Populations

for (p in pops) {
  df[[paste0(p,"_risk_allele_freq")]] <- ifelse(df$risk_allele == df$refallele, df[[paste0('refallele_freq-',p)]],df[[paste0('otherallele_freq-',p)]])
}

# Clean up total allele frequency.
df$total_allele_freq[is.nan(df$total_allele_freq)] <- NA


# Plot Risk Allele Frequencies
ggplot(df, aes(x=df$risk_allele_freq, y=df$'JPT_risk_allele_freq')) +
  geom_point(shape=1)      # Use hollow circles
  
  
# Check consistancy of risk and ref/other alleles
sum(df$risk_allele == df$refallele | df$risk_allele == df$otherallele)

# Date Conversions (ISO Standard)
df$Date.Added.to.Catalog <- as.Date(df$Date.Added.to.Catalog,"%m/%d/%Y")
df$Date <- as.Date(df$Date,"%m/%d/%Y")
df$Journal <- factor(df$Journal)


# Clean up Population Data ~ Would love to know if there is a better way to do this...
extract_study_population <- function(col) {
  x <- df[,col]
  m <- t(sapply(x,function(i) str_match_all(i,"^([0-9,]{3,7}) (.*?) ancestry (.*?), ([0-9,]{3,7}) (.*?) ancestry (.*)")[[1]][2:7]))
  m <- as.data.frame(m, stringsAsFactors = F) # AH - you have to transpose these!
  # Add Data Types
  m$V1 <- as.integer(gsub(",","",m$V1))
  m$V2 <- as.factor(m$V2)
  m$V3 <- as.factor(m$V3)
  
  m$V4 <- as.integer(gsub(",","",m$V4))
  m$V5 <- as.factor(m$V5)
  m$V6 <- as.factor(m$V6)
  # Remove GWAS studies that are *not* case/control
  m
}

pop_recodes <- c(
  "African" = "African",
  "African American" = "African",
  "African American, European ancestry, Hispanic/Latin American ancestry,  East Asian ancestry, South Asian" = NA,
  "anti-dsDNA positive European" = "European",
  "Ashkenazi Jewish" = NA,
  "bipolar disorder European" = "European",
  "British cases, 2,365 British controls, 1,145 European" = "European",
  "cases, 1,488 controls, 400 European and South Asian" = NA,
  "Chinese" = "Chinese",
  "Chinese anacestry individuals, 2,431 Malay" = NA,
  "comparatively younger European" = "European",
  "East Asian" = NA,
  "European" = "European",
  "European American" = "European",
  "European ancestry, African American ancestry and Native American" = NA,
  "European and other" = NA,
  "female European" = "European",
  "Han Chinese" = "Chinese",
  "Hispanic" = "Hispanic",
  "Indian" = "Indian",
  "Indo-European" = NA,
  "Indonesian" = NA,
  "Iranian" = NA,
  "irritable mania European" = "European",
  "Japanese" = "Japanese",
  "Jewish" = NA,
  "Korean" = NA,
  "Mexican" = "Mexican",
  "morbidly obese European" = "European",
  "Northern European" = "European",
  "Singaporean Chinese" = "Chinese",
  "Singaporean Chinese cases, 943 Singaporean Chinese controls, 297 Han Chinese cases, 1,044 Han Chinese controls, 573 South Asian" = NA,
  "South Asian" = NA,
  "Thai" = NA,
  "Vietnamese" = NA,
  "African Americans" = "African",
  "anti-dsDNA negative European" = "European",
  "Asian" = NA,
  "Asian Indian" = "Indian",
  "comparatively older European" = "European",
  "elated mania European" = "European",
  "European ancestry, African American" = NA,
  "European and South Asian" = NA,
  "Indian Asian" = "Indian",
  "lean European" = "European",
  "major depressive disorder European" = "European",
  "Malay" = NA,
  "Singaporean Malay" = NA
)

# Test
study_populations <- extract_study_population("Initial.Sample.Size")
study_populations$V2 <- as.factor(pop_recodes[study_populations$V2]) # Recodes Cases
study_populations$V5 <- as.factor(pop_recodes[study_populations$V5]) # Recodes Controls

# Unit Testing
allowed_pops <- c("African","Chinese","European","Hispanic","Indian","Japanese","Mexican")
checkTrue(all(levels(study_populations$V2) %in% allowed_pops), "Error in Population Parsing of Cases")
checkTrue(all(levels(study_populations$V5) %in% allowed_pops), "Error in Population Parsing of Controls")

# Case/Control may still need to be cleaned up.

# Extract case/control information.
df[ c("nCases","case_ancestry","is_case","nControl","control_ancestry","is_control")] <-extract_study_population("Initial.Sample.Size")

# Population Codes
"
ASW (A): African ancestry in Southwest USA
CEU (C): Utah residents with Northern and Western European ancestry from the CEPH collection
CHB (H): Han Chinese in Beijing, China
CHD (D): Chinese in Metropolitan Denver, Colorado
GIH (G): Gujarati Indians in Houston, Texas
JPT (J): Japanese in Tokyo, Japan
LWK (L): Luhya in Webuye, Kenya
MEX (M): Mexican ancestry in Los Angeles, California
MKK (K): Maasai in Kinyawa, Kenya
TSI (T): Toscans in Italy
YRI (Y): Yoruba in Ibadan, Nigeria (West Africa)
"

# Add Hapmap Population Mapping. Not very great; but as close as it will probably get...
hapmap_pop_matches <- c(African = "ASW", Chinese = "CHB", European = "CEU", Hispanic = "MEX", Indian = "GIH", Japanese = "JPT", Mexican = "MEX")

###############################
# Merge in Allele Frequencies #
###############################




#####################
# Merge in EFO data #
#####################
efo = read.csv("EFO/GWAS-EFO-Mappings201302.csv")
efo$efo_terms <- as.factor(ifelse(!grepl("Other *", efo$PARENT), as.character(efo$PARENT), as.character(efo$EFOTRAIT)))

#rename disease trait.
df$DISEASETRAIT <- df$Disease.Trait

# Merge into original data frame
dfm <- merge(df, efo, by = c("PUBMEDID","DISEASETRAIT"))
dfq <- merge(df, efo, by = c("DISEASETRAIT"))

# Generate frequency columns
efo <- ddply(efo, .(efo_terms), mutate, freq.efo = length(efo_terms))

