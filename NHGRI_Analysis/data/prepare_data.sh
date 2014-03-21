#!/bin/bash
# Set directory to current.
cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd

# Install Dependencies
brew install ssed
brew install vcftools
pip install csvkit

#---------------#
# Download Data #
#---------------#

#===========================#
# Download the GWAS Catalog #
#===========================#

# gwascatalog.txt will be timestamped and saved - as the most recent version.
wget --timestamping 'http://www.genome.gov/admin/gwascatalog.txt'

# Generate list of *strongest* associated SNPs.
iconv -c -f utf-8 -t ascii  gwascatalog.txt | cut -f 21 | egrep -o "rs[0-9]+" | sort -n | uniq > gwas_catalog_rs_list.txt

#=====================#
# Download MeSH Terms #
#=====================#
wget --timestamping 'ftp://nlmpubs.nlm.nih.gov/online/mesh/.meshtrees/mtrees2014.bin' -O MeSH/mtrees2014.txt

#===========================#
# Download GWAS-EFO mapping #
#===========================#
wget --timestamping --directory-prefix EFO 'http://wwwdev.ebi.ac.uk/fgpt/gwas/ontology/GWAS-EFO-Mappings201302.xlsx'
python in2csv EFO/GWAS-EFO-Mappings201302.xlsx > EFO/GWAS-EFO-Mappings201302.csv

#==========================#
# Download Omim Dataset    #
#==========================#
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/mim2gene.txt'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap.key'
wget --timestamping --directory-prefix omim 'ftp://ftp.omim.org/OMIM/genemap'

# Reformat the genemap file to get the disease name.
cut -f 10,13,14 -d "|" omim/genemap 

# Reformat OMIM dataset, more appropriate format.
# Removes 'disease only entries which have genes without entrez IDs (hence nothing to match on).
egrep "gene" omim/mim2gene.txt | egrep -e "-\t-$" -v | cut -f 1,3,4 > omim/omim.txt


#======================================================#
# Download Entrez - Gene Ontology (GO)                 #
# Mapping (Classifes genes by process, function, etc.) #
#======================================================#

# For a first pass - a wide format of GO term <--> Gene is produced regardless of evidence type.
wget --timestamping --directory-prefix GO 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz'
gunzip GO/*.gz -f
# Filter Human GO Annotation Terms, process file list.
egrep '(^#|^9606\t)' 'GO/gene2go' | sed 1d | cut -f '2,3,6' -d $'\t' - | sort -k 2 | cat <(echo -e 'Gene_ID\tGO_ID\tGO_term') - > 'GO/go_annotations.tmp'
# Reformat with Python Script
./GO/format.py
# Join the files
# Remove Temporary Files
rm -f GO/gene2go
=======
# Remove Temporary Files
rm GO/*.tmp

#===============================#
# Download KEGG Data (Pathways) #
#===============================#
# Download select files from UCSC (hg19)
wget --timestamping --directory-prefix kegg "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownToKeggEntrez.txt.gz"
wget --timestamping --directory-prefix kegg "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/keggMapDesc.txt.gz"
gunzip -f kegg/*.txt.gz

# Using super sed (ssed; install using 'brew install ssed')
# Format and join pathways with their description.
cut -f 3 kegg/knownToKeggEntrez.txt | ssed -e "s/+/\t/g" | sort -k 1 | join -1 1  -2 1  -t $'\t' - 'kegg/keggMapDesc.txt' | uniq | cat <(echo -e "kegg\tGene_ID\tpathway_desc") -  > kegg/entrez_kegg.txt

./kegg/format.py



#----------------------------------------------#
# Hapmap (downloaded from UCSC Genome Browser) #
#----------------------------------------------#
# Download all allele freq. files from hapmap.
wget -nd -r  -A "allele*.gz" -e robots=off --directory-prefix hapmap "http://hapmap.ncbi.nlm.nih.gov/downloads/frequencies/2010-08_phaseII+III/"
gunzip hapmap/*.gz # Unzip all the files
./hapmap/hapmap_create_sqlite.py # This generates an SQL file of the hapmap data.
./hapmap/reshape_hapmap_data.py # Generates a reshaped file with the Hapmap data.


#------------#
# 1k Genomes #
#------------#

wget --timestamping --directory-prefix 1kg 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz'
gunzip '1kg/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz'
vcftools --vcf 1kg/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf --snps gwas_catalog_rs_list.txt --recode-INFO-all  --out 1kg/nhgri_1kg 
python 1kg/vcf_to_tab.py --config 1kg/output.yaml 1kg/nhgri_1kg.recode.vcf 1kg/nhgri_1kg.tab


# Because the catalog is frequently updated - create an archive of prior versions, in case
# analysis needs to be verified with older versions.
GWAS_CKSUM=$(md5 -q gwascatalog.txt)
TDATE=$(date +%m%d%Y)
# Check if the file is different and save to the archive
# with md5 hash (to uni=quely identify) and todays date.
if [ $(ls catalog_archive | grep "${GWAS_CKSUM:0:5}" | wc -l) == 0 ]; then
	cp gwascatalog.txt "catalog_archive/${TDATE}.${GWAS_CKSUM:0:5}.gwas.catalog.txt"
	echo "New GWAS Catalog saved to archive: ${TDATE}.${GWAS_CKSUM:0:5}.gwas.catalog.txt"
fi



# Run an R Script for further data preperation.
DIR=$(pwd)
Rscript scripts/clean_data.r --args $DIR