# NHGRI Analysis
## Structure

* data/ - Scripts for initial data cleaning and the data itself.
* analysis/ - Scripts for various data analysis
* report/ - Integrated reports & discussion
* run.sh - Will run all data cleaning, analysis, and reporting. The goal is for the entire project to be reproducible by running this script.

## Background
__NHGRI Analysis__ is Daniel Cook's winter rotation project at NU. The project will make use of the [NHGRI catalog of published GWAS studies](http://www.genome.gov/26525384) and integrate data from additional data sources such as:

* [HapMap](http://www.hapmap.org)
* [http://www.genome.jp/kegg/pathway.html](KEGG Pathway Database)
* [http://www.cdc.gov/nchs/icd.htm](ICD Disease Classification Database)

## Objective
To examine in a statistically rigorous way how different factors used in GWAS studies contribute to their applicability. The analysis is in some ways broad - and allows for many questions to be answered. For example - how similar/different are diseases in which significant SNPs were found in the same pathway or the same gene? How does population, admixture, and Linkage Disequilibrium contribute to results? 

## Week 1 (1/5/14)
* [ X ] Setup Git Repo
* [ X ] Directory Setup
* [ X ] GWAS Catalog loaded into R

## Week 2 (1/12/14)

Data from external sources is being integrated for analysis.

* [ X ] Kegg (From UCSC)
* [ X ] Gene Ontology (GO)
* [ X ] OMIM (Online Mendelian Inheritance in Man)

# Week 3 (1/19/14)

This week I am continuing to work on data formatting. Hoping to have data ready for preliminary analysis. The work
is mostly focusing on integrating all the data rather than acquisition and base formatting.

* [ X ] Hapmap allele frequency reshaped and ready for merging with GWAS catalog.
* [ X ] GO Annotations reshaped - ready for joining with GWAS catalog. There are 586 GO terms being incorporated.
* [ X ] Kegg Pathways reshaped - ready for joining with GWAS catalog. There are 209 pathways.

# Week 4 (1/26/14)

* [ X ] Learning R.

# Week 5 (2/2/14)

* [ X ] EFO trait data successfully imported.

# Week 6 (2/9/14)

* [ X ] Hapmap strand issues have been resolved.


