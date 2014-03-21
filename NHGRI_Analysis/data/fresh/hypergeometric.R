library(xlsx)
library(stringr)

# Load data sets
gwascatalog <- read.delim("gwascatalog.txt")
efo <- read.xlsx("GWAS-EFO-Mappings201306.xlsx",1)

# Repair names
names(efo) <- str_replace(names(efo),"\\.","")

# Generate pubmed+disease key in efo and gwascatalog datasets
efo$pubdis <- paste(efo$DISEASETRAIT,efo$PUBMEDID,sep=" %_% ")
gwascatalog$pubdis <- paste(gwascatalog$Disease.Trait,gwascatalog$PUBMEDID,sep=" %_% ")

# Generate split of pubdis by trait
efo_trait <- split(efo$pubdis, efo$EFOTRAIT)
pub_gene <- split(gwascatalog$Snp_gene_ids, gwascatalog$pubdis)

# List Trait - Genes
as.numeric(unique(as.character(unlist(pub_gene[unlist(efo_trait[2])]))))

