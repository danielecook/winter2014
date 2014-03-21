library(xlsx)
library(stringr)

#-------------------#
# Load Initial Data #
#-------------------#

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
kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg <- lapply(as.list(kegg[mapped]), as.numeric)



SetMembers <- function(x) {
  # Returns all unique members of a list.
  as.numeric(unique(unlist(x)))
}

FilterEmpty <- function(l) {
  # Removes empty vectors within a list.
  l[lapply(l,length)>0]
} 


pairwise_analysis <- function(comp_set,urn_set) {

  # Generate complete, intersecting set
  i_set <- intersect(SetMembers(comp_set), SetMembers(urn_set))
  
  # Filter each set, retaining only elements in common.
  comp_set <- FilterEmpty(lapply(comp_set, function(c) c[c %in% i_set]))
  urn_set <- FilterEmpty(lapply(urn_set, function(c) c[c %in% i_set]))
  
  # Generate pairwise matrices!
    out <- lapply(comp_set, function(c1) {
      lapply(urn_set, function(c2) intersect(c1, c2))
    })
  out
}

comp_set <- kegg
urn_set <- trait_disease


d <- pairwise_analysis(kegg[1:20],trait_disease)

