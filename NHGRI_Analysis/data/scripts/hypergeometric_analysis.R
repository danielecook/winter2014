library(dplyr)
library(gplots)


DIR = commandArgs(trailingOnly = TRUE)[2]
if (is.na(DIR[1])) {
  user <- Sys.info()[["user"]]
  setwd(sprintf("/Users/%s/Documents/git/winter2014/NHGRI_Analysis/data/",user))
} else {
  setwd(DIR)
}





# Load Kegg pathways
ek <- read.delim("kegg/entrez_kegg.tmp")

# Get Gene and Pathway Counts
ek <- ek %.%
  group_by(Gene_ID) %.%
  mutate(gene_count = n())

ek <- ungroup(ungroup(ek) %.%
                group_by(kegg) %.%
                mutate(pathway_count = n()))

ek <- ek %.%
  group_by(kegg)

#----------------#
# Load EFO data  #
#----------------#

efo = read.csv("EFO/GWAS-EFO-Mappings201302.csv")
efo$efo_term <- as.factor(ifelse(!grepl("Other *", efo$PARENT), as.character(efo$PARENT), as.character(efo$EFOTRAIT)))
efo <- tbl_df(efo)
efo <- efo %.%
  group_by(efo_term) %.%
  mutate(efo_term_count = n()) %.%
  filter(efo_term_count >= 10)

#---------#
# Load df #
#---------#

# Load Gwas Catalog
df <- read.delim("gwascatalog.txt", header=T)
df$Gene_ID <- df$Snp_gene_ids
df$DISEASETRAIT <- df$Disease.Trait 

# Merge into original data frame
df<- left_join(efo, df, by = c("PUBMEDID","DISEASETRAIT"))

# Gene - EFO Term
gt <- df %.%
  group_by(Gene_ID) %.%
  filter(Snp_gene_ids != "") %.%
  mutate(gene_count = n()) %.%
  ungroup() %.%
  select(Snp_gene_ids,gene_count,efo_term) %.%
  arrange(efo_term)

gt <- select(gt, Snp_gene_ids, gene_count, efo_term)

gt <- gt %.%
  group_by(efo_term)

# Remove duplicates
gt <- gt[duplicated(gt),]

library(org.Hs.eg.db)
library(KEGG.db)
xx <- unlist(as.list(KEGGPATHID2NAME))
kegg <- org.Hs.egPATH2EG # Human only
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])


kegg_disease <- data.frame()
kegg_genes <- unique(unlist(kegg2, use.names =F))

for (path in names(kegg2)) {
  for (ef_term in levels(gt$efo_term)) {
    efo_genes <- gt[gt$efo_term == ef_term,c("Snp_gene_ids")]
    intersect_genes <- intersect(efo_genes,as.character(unlist(kegg2[path])))
    q <- length(intersect_genes)
  if (length(efo_genes) > 0) {
  m <- length( kegg2[[path]] )
  n <- length(kegg_genes) - m
  k <- length( intersect(kegg_genes,efo_genes ))
  kegg_disease <- rbind(kegg_disease,data.frame(p = phyper(q, m, n, k, lower.tail = F, log.p=T),
                                                q = q,
                                                k = k, # Balls Drawn
                                                m = m, # White in Urn
                                                N = n+m, # Total number!
                                                efo_count = length(efo_genes),
                                                path = xx[path],
                                                efo_term = ef_term,
                                                kegg_id = path,
                                                intersecting_genes = paste(intersect_genes, collapse=", ")
                                                #kegg_genes = paste(kegg2[[path]], collapse=", "),
                                                #efo_genes = paste(efo_genes, collapse=", "),
                                                ))
    }
  }
}

# Convert INF to next max
kegg_disease$p[kegg_disease$p == -Inf] <- min(kegg_disease$p[kegg_disease$p != -Inf])


library(reshape)

r <- melt(kegg_disease, id=c("path","efo_term"), measure.vars=c("p"))
f<- as.matrix(cast(r, efo_term ~ path, value='p',fill=0))

symm=F

heatmap.2 (f,
           
           # dendrogram control
           Rowv = TRUE,
           Colv=if(symm)"Rowv" else TRUE,
           distfun = dist,
           hclustfun = hclust,
           dendrogram = c("both","row","column","none"),
           symm = FALSE,
           
           # data scaling
           scale = c("none","row", "column"),
           na.rm=TRUE,
           
           # image plot
           revC = identical(Colv, "Rowv"),
           add.expr,
           
           # mapping data to colors
           breaks,
           symbreaks=min(x < 0, na.rm=TRUE) || scale!="none",
           
           # colors
           col="heat.colors",
           
           # block sepration
           colsep,
           rowsep,
           sepcolor="white",
           sepwidth=c(0.05,0.05),
           
           # cell labeling
           cellnote,
           notecex=1.0,
           notecol="cyan",
           na.color=par("bg"),
           
           # level trace
           trace=c("column","row","both","none"),
           tracecol="cyan",
           hline=median(breaks),
           vline=median(breaks),
           linecol=tracecol,
           
           # Row/Column Labeling
           margins = c(5, 5),
           ColSideColors,
           RowSideColors,
           cexRow = 0.2 + 1/log10(nr),
           cexCol = 0.2 + 1/log10(nc),
           labRow = NULL,
           labCol = NULL,
           srtRow = NULL,
           srtCol = NULL,
           adjRow = c(0,NA),
           adjCol = c(NA,0),
           offsetRow = 0.5,
           offsetCol = 0.5,
           
           # color key + density info
           key = TRUE,
           keysize = 1.5,
           density.info=c("histogram","density","none"),
           #denscol=tracecol,
           symkey = min(x < 0, na.rm=TRUE) || symbreaks,
           densadj = 0.25,
           
           # plot labels
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           
           # plot layout
           lmat = NULL,
           lhei = NULL,
           lwid = NULL,
)

# goEnrichment
#--------------------------------------------------------------------------------------------------
getTermToGeneMappings <- function (genes)
{
  goTermMap = mget (genes, env = GOLOCUSID2GO, ifnotfound = NA)
  geneNames = unique (names (goTermMap))
  
  result = new.env ()
  for (gene in geneNames) {
    goTerms = unique (names (goTermMap [gene][[1]]))
    for (goTerm in goTerms) {
      #print (paste (goTerm, gene))
      if (is.na (match (goTerm, ls (result))))
        geneList = c ()
      else
        geneList = get (goTerm, result)
      geneList = c (geneList, gene)
      assign (goTerm, geneList, envir=result)
    } # for goTerm
  } # for gene
  result 
}

# Get list of all genes (Entrez)
analyseGO <- function(genelist,geneUniverse,p.val=0.05,test.dir=c("over","under"), onto=c("BP","MF","CC"), conditional=FALSE){
  require(org.Hs.eg.db)
  require(GOstats)
  selectedGene <- unique(genelist)
  geneUniverse <- unique(geneUniverse)
  hypGobject <- new("GOHyperGParams",geneIds=genelist,universeGeneIds=geneUniverse, annotation="org.Hs.eg.db", ontology=onto,pvalueCutoff=p.val,testDirection=test.dir)
  file <- "results.html"
  if(conditional==FALSE){
    hypGresult <- hyperGTest(hypGobject)
    htmlReport(hypGresult,file=file)
    cat("GOstats test result exported to ",paste0(getwd(),"/",file),"\n")
  } else
  {
    conditional(hypGobject) <- TRUE
    hypGresult <- hyperGTest(hypGobject)
    htmlReport(hypGresult,file=file)
    cat("GOstats test result exported to ",paste0(getwd(),"/",file),"\n")
  }
}

analyseGO( genelist = df$Gene_ID, geneUniverse = unique(entrez_ids) , onto="MF", test.dir="over")



entrez_ids <- mappedkeys(org.Hs.egGO)

df$Gene_ID <- as.numeric(df$Snp_gene_ids) # Drops some.

# Plot global GO term view
plotGOprofiles <- function(genelist, level=2){
  require(goProfiles)
  require(org.Hs.eg.db)
  profileBP <- basicProfile(genelist,onto="BP",level=level, orgPackage="org.Hs.eg.db")
  plotProfiles(profileBP)
  printProfiles(profileBP)
  profileMF <-  basicProfile(genelist,onto="MF",level=level, orgPackage="org.Hs.eg.db")
  plotProfiles(profileMF)
  printProfiles(profileMF)
}

plotGOprofiles(filter(df, Gene_ID != "")$Gene_ID)

#---------------#
# Kegg Analysis #
#---------------#

# Load the Kegg Dataset
kegg <- read.delim("kegg/kegg_reshaped.txt")

df$Gene_ID <- df$Snp_gene_ids

dfk <- join(df,kegg, by = c("Gene_ID"))
