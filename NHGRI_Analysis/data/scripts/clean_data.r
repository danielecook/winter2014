#! /usr/bin/env Rscript

library(stringr)
library(reshape2)
library(annotate)
library(RUnit)
library(plyr)
library(ggplot2)
library(devtools)
source_url('https://gist.github.com/danielecook/8809319/raw/')


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

# Convert Risk Allele Frequency to Numeric
df$Risk.Allele.Frequency <- as.numeric(df$Risk.Allele.Frequency)

#GO = data.frame(read.table('GO/GO_reshaped.txt',header=T,sep='\t',stringsAsFactors=T, strip.white=T))

#---------------#
# Basic Cleanup #
#---------------#

# Date Conversions (ISO Standard)
df$Date.Added.to.Catalog <- as.Date(df$Date.Added.to.Catalog,"%m/%d/%Y")
df$Date <- as.Date(df$Date,"%m/%d/%Y")
df$Journal <- factor(df$Journal)

#--------#
# Hapmap #
#--------#

#
# BASIC CLEANUP
#

# Hapmap Allele Frequency Data
HMAF = data.frame(read.table('hapmap/hapmap_allele_freq_reshaped.csv',header=T,sep=",",stringsAsFactors=F, strip.white=T)) 

# Clean up the hapmap allele frequency data.
names(HMAF) <- gsub("\\.\\.","",gsub("\\.\\.\\.\\.","-",gsub("X..","",names(HMAF))))

# Special function or cleaning up data
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


#---------------------------------#
# Parse out strongest risk allele #
#---------------------------------#

df$risk_allele <- gsub(" ","",do.call(rbind,strsplit(df$Strongest.SNP.Risk.Allele,"-"))[,2])
# Cleanup problems
df$risk_allele[df$risk_allele=="?"] <- NA

# Many difficult to interpret risk alleles remain (e.g. Cxrs12880735), filter these out.
df$risk_allele[!df$risk_allele %in% c("A","C","T","G")]  <- NA

# Number Risk Alleles from the NHGRI catalog available.
length(df$risk_allele[!is.na(df$risk_allele)]) 
# 11,168 - but we will only be able to match on a few of these.

# Forward Strand
length(subset(df$risk_allele,df$risk_allele == df$refallele))
# Reverse Strand
length(subset(df$risk_allele,df$risk_allele == df$otherallele))
# Neither = 282
length(subset(df$risk_allele,df$risk_allele != df$refallele & df$risk_allele != df$otherallele))

#-----------------------------------#
# Generate Hapmap risk allele freqs #
#-----------------------------------#

# Generate Hapmap risk allele freq (for plotting purposes).
df$hm_risk_allele_count <- ifelse(df$risk_allele == df$refallele, rowSums(df[,grep("refallele_count", names(df))],na.rm=T), rowSums(df[,grep("otherallele_count", names(df))],na.rm=T))
df$hm_total_allele_count <- rowSums(df[,grep("totalcount",names(df))],na.rm=T)
df$hm_risk_allele_freq <- df$hm_risk_allele_count / df$hm_total_allele_count


#---------------------------#
# Load 1kg allele frequency #
#---------------------------#

kg <- read.csv("~/Documents/git/winter2014/NHGRI_Analysis/data/1kg/1kg_formatted.txt")

df <- merge(df,kg,by=c("rs"), all.x=T, all.y=F)

#-------------------------------#
# Deal with NHGRI Strand Issues #
#-------------------------------#

df$refallele <- df$ref_1kg
df$otherallele <- df$oth_1kg

df$chrom_strand <- 0

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

#--------------------------------#
# Plot 1000 g allele frequencies #
#--------------------------------#

# Flip allele frequencies to match risk alleles
df$AF <- ifelse(df$ref_1kg == df$risk_allele_forward, 1-df$AF, df$AF)

# Plot NHGRI vs. 1000 genomes

p <- qplot(df, x=df$AF, y=df$Risk.Allele.Frequency, main="1000 Genomes vs. NHGRI Risk Allele Freq", ylim = c(0,1), xlab="1000 Genomes Allele Frequency", ylab="NHGRI Reported Risk Allele Frequencies") 
p <- p + scale_color_manual(name="Predicted Strand",values=c("#0080ff","#cccccc")) + theme(panel.background = element_rect(fill='white', colour='black'))
p

# Plot Hapmap vs. 1000 genomes
p <- qplot(dfm, x=df$AF, y=df$hm_risk_allele_freq, main="1000 Genomes vs. NHGRI Risk Allele Freq", ylim = c(0,1), xlab="1000 Genomes Allele Frequency", ylab="NHGRI Reported Risk Allele Frequencies") 
p <- p + scale_color_manual(name="Predicted Strand",values=c("#0080ff","#cccccc")) + theme(panel.background = element_rect(fill='white', colour='black'))
p

# Count number of complete cases.
sum(complete.cases(df[,c('AF','Risk.Allele.Frequency')]))

#-------------------------------------------------#
# Generate per population risk allele frequencies #
#-------------------------------------------------#

# Generate Risk allele frequency
pops <- gsub('refallele_freq-','',names(df[,grep("refallele_freq",names(df))])) # Populations

for (p in pops) {
  df[[paste0(p,"_risk_allele_freq")]] <- ifelse(df$risk_allele_forward == df$refallele, df[[paste0('refallele_freq-',p)]],df[[paste0('otherallele_freq-',p)]])
}


#-------------------------#
# Hapmap Allele Frequency #
#-------------------------#

draw_plot <- function(title, var1 = 'hm_risk_allele_freq', var2 = 'Risk.Allele.Frequency', cfactor='chrom_strand') {
  p <- qplot(df, x=df[[var1]], y=df[[var2]], main=title, ylim = c(0,1), xlab="HapMap Computed Risk Allele Frequencies", ylab="NHGRI Reported Risk Allele Frequencies", color=factor(df[[cfactor]], exclude=0)) 
  p <- p + scale_color_manual(name="Predicted Strand",values=c("#0080ff","#cccccc")) + theme(panel.background = element_rect(fill='white', colour='black'))
  p 
}

draw_plot('title',var1='AF')
#--------------------------------------------------------------------------------------------#
# Examine NHGRI reported allele freq. vs. Hapmap allele freq. (Before and after strand flip) #
#--------------------------------------------------------------------------------------------#


# Plot flipped obs.
draw_plot("Before Flip")
ggsave(filename='../analysis/risk_allele_freq/allele_freq_comparison_before_flip.png', plot=last_plot(), width = 10, dpi = 150)


# Plot flipped obs.
draw_plot("After Flip")
ggsave(filename='../analysis/risk_allele_freq/allele_freq_comparison_after_flip.png', plot=last_plot(), width = 10, dpi = 150)

#--------------------------------------#
# Mark studies with large divergences. #
#--------------------------------------#

# Identify hapmap allele freq and reported allele freq discrepancies and flag.qplo
df$risk_alleles_resids <- residuals(lm(df$Risk.Allele.Frequency ~ df$AF, na.action=na.exclude))
df$risk_allele_flag <- 1
df$risk_allele_flag[!is.na(df$risk_alleles_resids) & abs(df$risk_alleles_resids) >0.3] <- 0.5

x <- names(df[,(grep('risk_allele_freq',names(df)))])
x <- x[1:length(x)-1]
df$emptyFreqs <- rowSums(is.na(df[,x]))
df$emptyFreqsFlag <- 0.5
df$emptyFreqsFlag[df$emptyFreqs < 5 & df$emptyFreqs != 11 ] <- 1.0

p <- qplot(df, x=df$hm_risk_allele_freq, y=df$Risk.Allele.Frequency, main="title", ylim = c(0,1), xlab="HapMap Computed Risk Allele Frequencies", ylab="NHGRI Reported Risk Allele Frequencies", color=factor(df$emptyFreqsFlag, exclude=11)) 
p <- p + scale_color_manual(name="Predicted Strand",values=c("#0080ff","#cccccc"))  + opts(panel.background = theme_rect(fill='white', colour='black'))
p 

draw_plot("Residuals > 0.3", cfactor = "risk_allele_flag")
ggsave(filename='../analysis/risk_allele_freq/allele_freq_resids.png', plot=last_plot(), width = 10, dpi = 150)


#--------------------------------#
#  Export Dataframe for Analysis #
#--------------------------------#

save(df ,file="hapmap/hapmap_data.Rda")

save(df ,file="df.Rda")





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



