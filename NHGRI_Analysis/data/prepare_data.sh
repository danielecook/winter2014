#!/bin/bash
# Set directory to current.
cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd

#################
# Download Data #
#################

# Download the GWAS Catalog
# gwascatalog.txt will be timestamped and saved - as the most recent version.
wget --timestamping 'http://www.genome.gov/admin/gwascatalog.txt'

# Download Omim Dataset
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/mim2gene.txt'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap.key'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap'

# Download Entrez - Gene Ontology (GO) Mapping (Classifes genes by process, function, etc.)
wget --timestamping --directory-prefix GO 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz'
gunzip 'GO/gene2go.gz'
grep '^9606\t' 'GO/gene2go' > 'GO/gene2go.human.txt' # Only extracts human GO terms.
rm gene2go # Remove this file as it is no longer needed.

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