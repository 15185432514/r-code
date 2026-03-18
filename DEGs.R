
####DEGs####
#设置工作目录
setwd("TCGA")
setwd("DEG") 
library(tidyverse)
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
counts<- read.table("lncRNA_counts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","Tumor","Normal"),levels = c("Normal","Tumor"))) %>% 
  column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)


dds <- DESeq(dds)

resultsNames(dds)

res <- results(dds)
save(res,file="DEG_lncRNA.Rda")
DEG <- as.data.frame(res)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

write.table(DEG, file = "exp_lncRNA_diff.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####Volcano plot####
setwd("TCGA")
setwd("DEG")
library(tidyverse)
exp <- read.table("lncRNA_tpms_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#open the DEG_lncRNA.RDA
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

DEG$logP <- -log10(DEG$padj)
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()


ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()

#Add a dividing line
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")#根据需要调整（-1，1）
dev.off()


##Converting differential gene expression data formats
setwd("TCGA")
setwd("DEG")
#open the DEG_lncRNA.RDA

tpms <- read.table("lncRNA_tpms_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

DEG <- as.data.frame(res)
## select Log2FC>2
logFC_cutoff <- 2
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)

library(tidyverse)
a <- filter(DEG,change == 'UP')
up_down <- a
#Obtain the expression matrix of differentially expressed genes
up_down <- rownames_to_column(up_down,var = 'symbol')
tpms <- rownames_to_column(tpms,var = 'symbol')
DEG_UP_down <- inner_join(up_down,tpms,"symbol")

rownames(DEG_UP_down) <- NULL
DEG_UP_down <- DEG_UP_down %>% column_to_rownames("symbol")

DEG_UP_down <-DEG_UP_down [,-c(1:7)]

DEG_UP_down <- t(DEG_UP_down)

DEG_UP_down <-as.data.frame(DEG_UP_down)
DEG_UP_down <- rownames_to_column(DEG_UP_down,var='sample')
x <- ifelse(str_detect(DEG_UP_down$sample,"01A"),"1","0")
DEG_UP_down$group <- x

rownames(DEG_UP_down) <- NULL
DEG_UP_down <- DEG_UP_down %>% column_to_rownames("sample")

write.csv(DEG_UP_down, file = "2UP.csv")
write.table(DEG_UP_down,"2UP.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
