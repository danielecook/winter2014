# Set to current directory
cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd
pip install csvkit

#===========================#
# Download the GWAS Catalog #
#===========================#

# gwascatalog.txt will be timestamped and saved - as the most recent version.
wget --timestamping 'http://www.genome.gov/admin/gwascatalog.txt'

wget --timestamping --directory-prefix EFO 'http://www.ebi.ac.uk/fgpt/gwas/ontology/GWAS-EFO-Mappings201306.xlsx'
in2csv EFO/GWAS-EFO-Mappings201306.xlsx > EFO/GWAS-EFO-Mappings201306.csv