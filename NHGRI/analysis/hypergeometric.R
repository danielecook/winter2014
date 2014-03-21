library(stringr)
library(gplots)

# Load Functions
source("scripts/functions.R")

#-------------------#
# Load Initial Data #
#-------------------#

# Load data sets
gwascatalog <- read.delim("data/gwascatalog.txt")
efo <- read.csv("data/EFO/GWAS-EFO-Mappings201306.csv")

# Generate pubmed+disease key in efo and gwascatalog datasets
efo$pubdis <- paste(efo$DISEASETRAIT,efo$PUBMEDID,sep=" %_% ")
gwascatalog$pubdis <- paste(gwascatalog$Disease.Trait,gwascatalog$PUBMEDID,sep=" %_% ")

# Generate split of pubdis by trait
efo_trait <- split(efo$pubdis, efo$EFOTRAIT)
pub_gene <- split(gwascatalog$Snp_gene_ids, gwascatalog$pubdis)

lookup_genes <- function(efo_trait) {
  r <- as.numeric(unique(as.character(unlist(pub_gene[unlist(efo_trait)]))))
  r[!is.na(r)] # Filter out NAs
}

# Generate Trait-Gene List
trait_disease <- Filter(length,lapply(efo_trait, lookup_genes))

# Generate a complete list of genes associated with traits
efo_genes <- as.numeric(unique(unlist(trait_disease)))

#---------------#
# Kegg Analysis #
#---------------#
library(org.Hs.eg.db)
library(KEGG.db)
kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg <- lapply(as.list(kegg[mapped]), as.numeric)
kegg_names <- unlist(as.list(KEGGPATHID2NAME))
# Fix Kegg Names (Use actual name instead of Id)
names(kegg) <- kegg_names[names(kegg)]

# Perform pairwise analysis.
keggr <- pairwise_analysis(kegg,trait_disease)

# Plot heatmaps
pdf(file = "analysis/kegg/kegg_ALL.pdf", width=20, height = 15, pointsize=10)
heatmap.2(keggr,  margins=c(16,16), col=terrain.colors(100), xlab="Kegg Pathway", ylab="EFO Term", trace=c("none"))
dev.off()

# Break down matrix by EFO PARENT
efo_parent <- split(efo$EFOTRAIT,efo$PARENT)

# Filter efo terms by parent, and plot subsetted heatmaps
for(parent in names(efo_parent)) {
  t <- keggr[intersect(as.character(unlist(efo_parent[parent],use.names=F)),rownames(keggr)),]
  pdf(file = sprintf("analysis/kegg/kegg_%s.pdf",parent), width=20, height = 15, pointsize=10)
  heatmap.2(t,  margins=c(16,16), col=terrain.colors(100), title=sprintf("Kegg x EFO Parent (%s)",parent), xlab="Kegg Pathway", ylab="EFO Term", trace=c("none"))
  dev.off()
}

#----------------------#
# Pubmed Results Count #
#----------------------#
library(RISmed)

pub_count <- function(q1,q2) {
  # Counts the number of publications using 2 query words.
  QueryCount(EUtilsSummary(paste(q1,q2,sep=" "), type="esearch", db="pubmed"))
}

pub_results <- sapply(names(comp_set[1]), USE.NAMES=T, function(c1) {
  sapply(names(urn_set), USE.NAMES=T, function(c2) {
         print(c1)
         print(c2)
         pub_count(c1, c2)
         })
})

