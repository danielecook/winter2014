library(dplyr)
library(gplots)
library(ggplot2)
library(stringr) 
library(reshape2)


#-----------#
# Load Data #
#-----------#
setwd('data')

# Load Kegg pathways
EK <- read.delim("kegg/entrez_kegg.txt")

# GO
GO <- read.delim("~/Documents/git/winter2014/NHGRI_Analysis/data/GO/GO_reshaped.txt")

# GWAS Catalog
df <- read.delim("gwascatalog.txt", header=T)
    # Cleanup Catalog
    df$DISEASETRAIT <- df$Disease.Trait 

# Split obs. with multiple genes.
e <- melt(str_split_fixed(df$Snp_gene_ids,";",n=10))
e<-rename(e[e$value != "",c("X1","value")],c(value="Gene_ID"))
df<- merge(df,e,by.x="row.names",by.y="X1",all.x=T)

# EFO
efo = read.csv("EFO/GWAS-EFO-Mappings201302.csv")

#------------#
# Merge Data #
#------------#

# Merge into original data frame
df <- left_join(efo, df, by = c("PUBMEDID","DISEASETRAIT"))
df$efo_term <- as.factor(df$EFOTRAIT)

# Gene - EFO Term
gt <- df %.%
  filter(Snp_gene_ids != "") %.%
  select(Gene_ID,efo_term,PUBMEDID) %.%
  arrange(efo_term) %.%
  unique() %.%
  group_by(efo_term) %.%
  mutate(efterm_count = n()) %.%
  filter(efterm_count > 10) 

# Filter non-unique efo terms
gt$efo_term <- as.character(gt$efo_term)

# Match with GO terms

GO_disease <- data.frame()
GO_genes <- as.numeric(unique(GO$Gene_ID)

# Reshape Function
rec <- function(df, id.vars = c("path","efo_term"), form = "path ~ efo_term", val = "p",...) {
  m <- melt(df[,c(id.vars,val,...)],id.vars=id.vars)
  as.matrix(cast(m, as.formula(form), value=val, fill=0, fun.aggregate=mean))
}


gt$Gene_ID <- as.numeric(gt$Gene_ID)

#----------------#
# GO Annotations #
#----------------#
GO_NAMES <- names(GO)[2:length(names(GO))]
for (path in GO_NAMES) {
  print(path)
  print(proc.time()-ptm)
  ptm <- proc.time()
  for(ef_term in unique(gt$efo_term))   {
    efo_genes <- as.numeric(gt[gt$efo_term == ef_term,c("Gene_ID")])
    intersect_genes <- intersect(efo_genes,GO[GO[[path]] != "","Gene_ID"])
    q <- length(intersect_genes)
    k <- length( intersect( GO_genes ,efo_genes ))
    # Don't include cases where k = 1
    if (length(efo_genes) > 0 & q > 5 & k != q & k > 1) {
      m <- length(GO[GO[[path]]!="",path])
      n <- length(GO_genes) - m
      pubs <- unique(gt[gt$efo_term == ef_term & gt$Gene_ID %in% intersect_genes,"PUBMEDID"])
      GO_disease <- rbind(GO_disease,data.frame(p = phyper(q, m, n, k, lower.tail = F, log.p=T),
                                                q = q,
                                                k = k, # Balls Drawn
                                                m = m, # White in Urn
                                                N = n+m, # Total number!
                                                efo_count = length(efo_genes),
                                                efo_term = ef_term,
                                                GO_id = path,
                                                path = GO[GO[[path]] != "",path][1],
                                                intersecting_genes = paste(intersect_genes, collapse=", "),
                                                pubs = paste(pubs, collapse=", ")
                                                #kegg_genes = paste(kegg2[[path]], collapse=", "),
                                                #efo_genes = paste(efo_genes, collapse=", "),
      ))
    }
  }
}




save(GO_disease,file="GO.data")


#----------#
# Kegg Set #
#----------#
library(org.Hs.eg.db)
library(KEGG.db)
xx <- unlist(as.list(KEGGPATHID2NAME))
kegg <- org.Hs.egPATH2EG # Human only
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])

kegg_disease <- data.frame()
kegg_genes <- unique(unlist(kegg2, use.names =F))

# Generate p values - hypergeometric test!
for (path in names(kegg2)) {
  for (ef_term in unique(gt$efo_term)) {
    efo_genes <- gt[gt$efo_term == ef_term,c("Snp_gene_ids")]
    intersect_genes <- intersect(efo_genes,as.character(unlist(kegg2[path])))
    q <- length(intersect_genes)
    k <- length( intersect(kegg_genes,efo_genes ))
    # Don't include cases where k = 1
    if (length(efo_genes) > 0 & q > 0 & k != q & k > 1) {
      m <- length( kegg2[[path]] )
      n <- length(kegg_genes) - m
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



#---------------------#
# Plot Results by Pub #
#---------------------#

# Split obs. with multiple genes.
e <- melt(str_split_fixed(GO_disease$pubs,", ",n=40 ))
e<-rename(e[e$value != "",c("X1","value")],c(value="PUB"))
GO_BY_PUB<- merge(GO_disease,e,by.x="row.names",by.y="X1",all.x=T)

GO_PUB_LIST <- GO_BY_PUB[,c("efo_term","PUB","GO_id","path")]
GO_PUB_LIST <- arrange(GO_PUB_LIST, efo_term,PUB)
pub_order <- intersect(unique(GO_PUB_LIST$PUB),dimnames(r)[[2]] )

t<-unique(GO_PUB_LIST[,c("efo_term","PUB")])
t<-group_by(t, efo_term)
breaks<-as.numeric(row.names(t))
r <- rec(GO_BY_PUB, id.vars=c("PUB","path","efo_term"), form = "path ~ PUB", val = "p")

r <- r[,pub_order]

# Reorder matrix appropriately.

png(file = "../analysis/heatmaps/pub_results_sorted_broken.png", width=4000, height=2000,bg = "transparent")
heatmap.2(r[1:50,1:556], col=terrain.colors(100), na.rm = T, Colv=F, dendrogram = c("row"), colsep=breaks[2:length(breaks)], trace=c("none"),sepwidth=c(0.05,0.00), sepcolor="blue", margins=c(16,16), ylab="pathway", xlab="PUB", main="Sorted by EFO Term")
dev.off()


#---------------#
# Gen Heat Maps #
#---------------#

GO_pvals <- as.matrix(cast(melt(GO_disease[,c("path","efo_term","p")],id.vars=c("path","efo_term")), efo_term ~ path, value='p', fill=0))

# Plot heatmap GO pvals
png(file = "../analysis/heatmaps/GO_pvals.png", width=2500, height=1500, pointsize= 18,bg = "transparent")
heatmap.2(GO_pvals, na.rm=T, margins=c(16,16),  ylab="efo_term", xlab="pathway", main="Pvals: GO",trace=c("none"), col=terrain.colors(100))
dev.off()

# Plot heatmap GO pvals
png(file = "../analysis/heatmaps/GO_pvals_scale.png", width=2500, height=1500, pointsize= 18,bg = "transparent")
heatmap.2(GO_pvals, scale=c("row"), na.rm=T, margins=c(16,16),  ylab="efo_term", xlab="pathway", main="Pvals: GO",trace=c("none"), col=terrain.colors(100))
dev.off()

GO_pvals <- as.matrix(cast(melt(GO_disease[,c("path","efo_term","p")],id.vars=c("path","efo_term")), efo_term ~ path, value='p', fill=0))


# Plot like groups
efo$efo_term <- efo$EFOTRAIT

ematch <- efo[,c("efo_term","PARENT","PUBMEDID")]

GO_BY_PUB$PUBMEDID <- GO_BY_PUB$PUB

m <- merge(GO_BY_PUB,ematch, by=c("PUBMEDID","efo_term"))

for ( i in levels(m$PARENT)) { 
  print(i)
GO_pvals <- as.matrix(cast(melt(m[m$PARENT=="Cancer",c("path","efo_term","p")],id.vars=c("path","efo_term")), efo_term ~ path, value='p', fill=0, fun.aggregate=mean))

png(file = sprintf("../analysis/heatmaps/GO_pvals_%s.png", i), width=2500, height=1500, pointsize= 18,bg = "transparent")
heatmap.2(GO_pvals, scale=c("row"), na.rm=F, margins=c(16,16),  ylab="efo_term", xlab="pathway", main=sprintf("Pvals GO %s", i ),trace=c("none"), col=terrain.colors(100))
dev.off()

}


#------------------------------#
#  Break it down by population #
#------------------------------#








# Convert INF to next max
kegg_disease$p[kegg_disease$p == -Inf] <- min(kegg_disease$p[kegg_disease$p != -Inf])

library(reshape)
# Generate a matrix of p-values
kegg_disease_assoc <- as.matrix(cast(
                  melt(kegg_disease, id=c("path","efo_term"), measure.vars=c("p"))
                  , efo_term ~ path, value='p',fill=0))


# Remove general terms
kegg_disease <- filter(kegg_disease, efo_term != "Body measurement" & efo_term != "Biological process")

#-----------------------------------------#
# Kegg ~ efo_term | Pubmed Search Results #
#-----------------------------------------#
# Load
pubmed_r <- melt(read.csv("r_search_results.csv", check.names=F), id.vars=c("efo_term"))
# Rename path variable
colnames(pubmed_r)[2] <- "path"

# Merge pubmed result count with kegg disease
q <- merge(kegg_disease, pubmed_r, by=c("path","efo_term"), type = "left", match="all")


# Plot search results vs. genetic links.
qplot(-q$p, q$value, ylab="Search Results", xlab="-log10(p-val) ~ hypergeometric test; Kegg Pathway", main="(Kegg x Disease) pubmed search results vs. (Kegg x Disease) gene associations", colour=q$efo_term)

# Do some basic filtering
q <- q %.%
  filter(efo_term != "Biological process" & efo_term != "Body measurement" )

# Generate residuals for plotting.
q$resids <- residuals(lm(-q$p ~ q$val))

# Reshape residuals and plot heatmap
pub_heat <- as.matrix(cast(melt(q[,c("path","efo_term","resids")],id.vars=c("path","efo_term")), efo_term ~ path, value='resids', fill=0))
result_c <- as.matrix(cast(melt(q[,c("path","efo_term","value")],id.vars=c("path","efo_term")), efo_term ~ path, value='value', fill=0))
pvals <- as.matrix(cast(melt(q[,c("path","efo_term","p")],id.vars=c("path","efo_term")), efo_term ~ path, value='p', fill=0))


GO_pvals <- as.matrix(cast(melt(GO_disease[,c("GO_name","efo_term","p")],id.vars=c("GO_name","efo_term")), efo_term ~ GO_name, value='p', fill=0))

# Plot heatmap GO pvals
png(file = "../analysis/heatmaps/GO_pvals.png", width=2500, height=1500, pointsize= 18,bg = "transparent")
heatmap.2(GO_pvals, na.rm=T, margins=c(16,16),  ylab="efo_term", xlab="pathway", main="Pvals: GO",trace=c("none"), col=terrain.colors(100))
dev.off()


png(file = "../analysis/heatmaps/pubmed_search_results3.png", width=2500, height=1500, pointsize= 18,bg = "transparent")
heatmap.2(log10(result_c+1), trace=c("none"),  margins=c(16,16),  ylab="efo_term", xlab="pathway", main="Result Count (Kegg & Disease) : Pubmed", col=terrain.colors(100), Rowv=F)
dev.off()

png(file = "../analysis/heatmaps/pub_results_vs_pval3.png", width=2500, height=1000, bg = "transparent")
heatmap.2(pub_heat, scale=c("row"),main="Residual of Corr(enrichment (disease x Gene) ~ pubmed count)",  margins=c(16,16),  ylab="efo_term (scaled on row)", xlab="pathway",title="Publication Results", trace=c("none"), col=terrain.colors(100,alpha=1))
dev.off()

png(file = "../analysis/heatmaps/hyper_pvals2.png", width=2500, height=1000,bg = "transparent")
heatmap.2(pvals, margins=c(16,16), ylab="efo_term", xlab="pathway", main="hypergeometric p-val, scaled by row", trace=c("none"),col=terrain.colors(100))
dev.off()

png(file = "../analysis/heatmaps/pub_results_log.png", width=2500, height=1000,bg = "transparent")
heatmap.2(pvals, na.rm = T, margins=c(16,16), ylab="efo_term", xlab="pathway", main="hypergeometric p-val, scaled by row", trace=c("none"),col=terrain.colors(100))
dev.off()


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
