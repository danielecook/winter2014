#!/bin/bash
# Set directory to current.
cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd

#################
# Download Data #
#################

# Download the GWAS Catalog
#==========================#

# gwascatalog.txt will be timestamped and saved - as the most recent version.
wget --timestamping 'http://www.genome.gov/admin/gwascatalog.txt'

# Download Omim Dataset
#==========================#

wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/mim2gene.txt'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap.key'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap'

# Reformat OMIM dataset

# Download Entrez - Gene Ontology (GO) Mapping (Classifes genes by process, function, etc.)
#==========================#
wget --timestamping --directory-prefix GO 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz'
gunzip 'GO/gene2go.gz'
egrep '(^#|^9606\t)' 'GO/gene2go' | sed 1d > 'GO/gene2go.human.tmp' # Only extracts human GO terms.
echo 'tax_id\tGeneID\tGO_ID\tEvidence\tQualifier\tGO_term\tPubMed Category' | cat - 'GO/gene2go.human.tmp' > 'GO/gene2go.human.txt'
# Remove Temporary Files
rm GO/gene2go.human.tmp
rm GO/gene2go
rm GO/gene2go.gz

# Download KEGG Data (Pathways)
#==============================#
# Download select files from UCSC (hg19)

for var in keggPathway KeggMapDesc knownGene kgXref
do
wget --timestamping --directory-prefix test 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/$var.txt.gz'
gunzip kegg/$var.txt.gz
done

# Join kegg pathway description with ID; keggMapDesc is already sorted; kgXref has 82,960 lines.
sort kegg/keggPathway.txt -k3 | join -1 3 -2 1 -t $'\t' - kegg/keggMapDesc.txt | cut -f 1,2,4 | sort -k 2 > kegg/kegg_tmp.txt # 58,073 lines.
sort kegg/KgXref.txt -k 1 | join -1 1 -2 2 -t $'\t' - kegg/kegg_tmp.txt > kegg/kegg_merged.txt
rm kegg/kegg_tmp.txt

# Because the catalog is frequently updated - create an archive of prior versions, in case
# analysis needs to be verified with older versions.
GWAS_CKSUM=$(md5 -q gwascatalog.txt)
TDATE=$(date +%m%d%Y)
# Check if the file is different and save to the archive
# with md5 hash (to uniquely identify) and todays date.
if [ $(ls catalog_archive | grep "${GWAS_CKSUM:0:5}" | wc -l) == 0 ]; then
	cp gwascatalog.txt "catalog_archive/${TDATE}.${GWAS_CKSUM:0:5}.gwas.catalog.txt"
	echo "New GWAS Catalog saved to archive: ${TDATE}.${GWAS_CKSUM:0:5}.gwas.catalog.txt"
fi

# Run an R Script for further data preperation.
DIR=$(pwd)
Rscript scripts/clean_data.r --args $DIR