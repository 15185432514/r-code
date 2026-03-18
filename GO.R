
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("clusterProfiler", quietly = TRUE)) {
 
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}

if (!require("org.Hs.eg.db", quietly = TRUE)) {
 
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}

if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!require("ggsci", quietly = TRUE)) {
  install.packages("ggsci")
}

library(clusterProfiler) 
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(ggsci)



upGene=read.table("corResult.txt",header=T, sep="\t", check.names=F)
upGene=upGene[!duplicated(upGene$mRNA),]
rownames(upGene)=upGene[,1]
upGene=upGene[,-1]
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", 
              toType = c("ENTREZID"), 
              OrgDb = org.Hs.eg.db)

# Enrichment analysis 
GO <- enrichGO(gene = Genes$ENTREZID, 
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "all",     
               pAdjustMethod = "BH", 
               pvalueCutoff = 0.05,   
               qvalueCutoff = 0.05, 
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   
GO_result <- as.data.frame(GO)



barplot(GO)
barplot(GO, drop = TRUE, 
        showCategory =6,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')


