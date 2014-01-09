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