# This script imports EFO data and cleans it as is appropriate.

# Has built in EFO services!
library("ontoCAT")
library(plyr)
efo <- getEFO() # Downloads the latest version of EFO Ontology

GWAS_EFO = read.csv("EFO/GWAS-EFO-Mappings201302.csv")

# Pull out EFO Trait
GWAS_EFO$EFOURI <- sub("http://www.ebi.ac.uk/efo/","",GWAS_EFO$EFOURI)

# Unfortunately, they mix a number of ontologies.
GWAS_EFO$direct_parent <- factor(sapply(GWAS_EFO$EFOURI,failwith(NA,function(x) getLabel(getTermParentsById(efo,x)[[1]]))))

GWAS_EFO <- arrange(GWAS_EFO,direct_parent)


corder <- function(df,cols) {
  df[,c(cols,setdiff(names(df),cols))]
}

"
DISEASETRAIT
EFOTRAIT
EFOURI
PARENT
PUBMEDID
AUTHOR
PUBDATE
JOURNAL
"