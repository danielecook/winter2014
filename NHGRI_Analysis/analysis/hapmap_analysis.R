
library(ggplot2)
library(reshape2)


DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  user <- Sys.info()[["user"]]
  setwd(sprintf("/Users/%s/Documents/git/winter2014/NHGRI_Analysis/data/",user))
} else {
  setwd(DIR)
}

# Load data
load(file="hapmap/hapmap_data.Rda")

#-------------------#
# Merge in EFO data #
#-------------------#
efo = read.csv("EFO/GWAS-EFO-Mappings201302.csv")
efo$efo_terms <- as.factor(ifelse(!grepl("Other *", efo$PARENT), as.character(efo$PARENT), as.character(efo$EFOTRAIT)))

#rename disease trait.
df$DISEASETRAIT <- df$Disease.Trait

# Merge into original data frame
dfm <- join(df, efo, by = c("PUBMEDID","DISEASETRAIT"),type = "full")
dfq <- merge(df, efo, by = c("DISEASETRAIT"))

# Generate frequency columns
efo <- ddply(efo, .(efo_terms), mutate, freq.efo = length(efo_terms))

for (i in unique(dfm$PARENT)) {
  dfm_s <- subset(dfm,dfm$PARENT == i)
  p <- qplot(dfm_s, x=dfm_s$hm_risk_allele_freq, y=dfm_s$Risk.Allele.Frequency, main=i, xlim = c(0,1), xlab="HapMap Computed Risk Allele Frequencies", ylab="log(OR)", color=factor(dfm_s$PARENT), size=dfm_s$OR.or.beta )
  p <- p + theme(panel.background = element_rect(fill='white', colour='black'))
  p
  ggsave(filename=sprintf('../analysis/risk_allele_freq/by_group_%s.png',i), plot=last_plot(), width = 10, dpi = 150)
}


download.file("https://raw.github.com/stephenturner/qqman/master/qqman.r", destfile="./qqman.r", method="curl")
source("./qqman.r")

#---------------------------#
# Manhattan Plot by Disease #
#---------------------------#

# Setup Variables
dfm$P <- as.numeric(dfm$p.Value)
dfm$CHR <- dfm$Chr_id
dfm$BP <- dfm$Chr_pos
# Subset
dfmanhattan <-subset(dfm,!is.na(dfm$P))
manhattan(dfmanhattan)

for (i in unique(dfm$PARENT)) {
  png(sprintf('../analysis/risk_allele_freq/manhattan_%s.png',i), width=700)
  manhattan(subset(dfmanhattan,dfmanhattan$PARENT == i),main=i)
  dev.off()
  #ggsave(filename=sprintf('../analysis/risk_allele_freq/manhattan_%s.png',i), plot=p, width = 10, dpi = 150)
}

manhattan(subset(dfmanhattan,dfmanhattan$PARENT == "Cancer"))

#----------------------------------#
# Hypergeometric Tests - Kegg      #
#----------------------------------#

# Import Kegg Data

