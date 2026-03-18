
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


upGene=read.table("corResult857.txt",header=T, sep="\t", check.names=F)
upGene=upGene[!duplicated(upGene$mRNA),]
rownames(upGene)=upGene[,1]
upGene=upGene[,-1]
Genes <- bitr(rownames(upGene),  
              fromType = "SYMBOL", 
              toType = c("ENTREZID"),
              OrgDb = org.Hs.eg.db) 

# Enrichment analysis 
KEGG <- enrichKEGG(gene = Genes$ENTREZID,
                   organism = "hsa", 
                   keyType = "kegg", 
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

barplot(KEGG,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10,
        title = "KEGG_enrichment")
dotplot(KEGG)