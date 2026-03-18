setwd('TCGA')
setwd("data")
library(tidyverse)
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("remotes")
BiocManager::install("ExperimentHub")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
cancer_type = "TCGA-LUAD"  
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file = "luad.gdc_2022.rda") 

####clinical data####
setwd("TCGA")
setwd("data")
library(tidyverse)
load("luad.gdc_2022.rda")
clinical <- as.data.frame(expquery2@colData) %>%   
  .[!duplicated(.$sample),]

clinical <- filter(clinical,prior_treatment == 'No')
clinical <- filter(clinical,prior_malignancy== 'no')

clinical <-clinical[,c("gender","age_at_index","years_smoked","ajcc_pathologic_stage",
                       "ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m")]

class(clinical$gender)
class(clinical$age_at_index)
class(clinical$years_smoked)
class(clinical$ajcc_pathologic_stage)
class(clinical$ajcc_pathologic_t)
class(clinical$ajcc_pathologic_n)
class(clinical$ajcc_pathologic_m)

table(clinical$gender)
table(clinical$age_at_index)
table(clinical$years_smoked)
table(clinical$ajcc_pathologic_stage)
table(clinical$ajcc_pathologic_t)
table(clinical$ajcc_pathologic_n)
table(clinical$ajcc_pathologic_m)

clinical$ajcc_pathologic_stage <- gsub("A","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_stage <- gsub("B","",clinical$ajcc_pathologic_stage)
clinical$ajcc_pathologic_t <- gsub("a","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_t <- gsub("b","",clinical$ajcc_pathologic_t)
clinical$ajcc_pathologic_m <- gsub("a","",clinical$ajcc_pathologic_m)
clinical$ajcc_pathologic_m <- gsub("b","",clinical$ajcc_pathologic_m)

write.table(clinical,"clinical.all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####counts datas####
setwd('TCGA')
setwd("data")
library(tidyverse)
load("luad.gdc_2022.rda")
load("gene_annotation_2022.rda")
table(gene_annotation_2022$type)

counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES

counts <- counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]


rownames(counts) <- NULL
counts <- counts %>% column_to_rownames("symbol") 

# Retain lncRNA
table(counts$type)
counts <- counts[counts$type == "lncRNA",]
counts <- counts[,-c(1,ncol(counts))]

ncol(counts)
nrow(counts)
#Replace the sample name in clinical with the sample name in counts
clinical <- read.table("clinical.all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
counts <-counts[,rownames(clinical)]

colnames(counts) <- substring(colnames(counts),1,16)
counts <- counts[,!duplicated(colnames(counts))]
table(substring(colnames(counts),14,16))
##Data filtering
counts <- counts[apply(counts, 1, function(x) sum(x >1) > 0.5*ncol(counts)), ]

counts01A <- counts[,substring(colnames(counts),14,16) == c("01A")]

counts11A <- counts[,substring(colnames(counts),14,16) == c("11A")]
table(substring(colnames(counts01A),14,16))
table(substring(colnames(counts11A),14,16))

####tpms datas####
tpms <- expquery2@assays@data@listData[["tpm_unstrand"]]
colnames(tpms) <- expquery2@colData@rownames
rownames(tpms) <- expquery2@rowRanges@ranges@NAMES
tpms <- tpms %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation_2022,"ENSEMBL") %>% 
  .[!duplicated(.$symbol),]
rownames(tpms) <- NULL
tpms <- tpms %>% column_to_rownames("symbol") 
# Retain lncRNA
tpms <- tpms[tpms$type == "lncRNA",]
tpms <- tpms[,-c(1,ncol(tpms))]

clinical <- read.table("clinical.all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tpms <- tpms[,rownames(clinical)]

colnames(tpms) <- substring(colnames(tpms),1,16)
tpms <- tpms[,!duplicated(colnames(tpms))]
##Data filtering
tpms =tpms[apply(tpms, 1, function(x)sum(x > 0) >0.5*ncol(tpms)), ] # tpm用这个。

tpms01A <- tpms[,substring(colnames(tpms),14,16) == c("01A")]

tpms11A <- tpms[,substring(colnames(tpms),14,16) == c("11A")]


write.table(counts01A,"lncRNA_counts01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(counts11A,"lncRNA_counts11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A,"lncRNA_tpms01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A,"lncRNA_tpms11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

counts <- cbind(counts01A,counts11A)
tpms <- cbind(tpms01A,tpms11A)
write.table(counts,"lncRNA_counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms,"lncRNA_tpms.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####tpms_log2####
range(tpms)
range(tpms01A)
range(tpms11A)
tpms_log2 <- log2(tpms+1)
range(tpms_log2)
tpms01A_log2 <- log2(tpms01A+1)
range(tpms01A_log2)
tpms11A_log2 <- log2(tpms11A+1)
range(tpms11A_log2)

write.table(tpms_log2,"lncRNA_tpms_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms01A_log2,"lncRNA_tpms01A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(tpms11A_log2,"lncRNA_tpms11A_log2.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####principal component analysis ####

data <- read.table("lncRNA_tpms_log2.txt", header = TRUE, row.names = 1)

data.pca <- prcomp(t(data), scale. = TRUE)
pc_scores <- data.pca$x
pca_data <- data.frame(PC1 = pc_scores[, 1], PC2 = pc_scores[, 2], Label = ifelse(grepl("01A$", rownames(pc_scores)), "Tumor", "Normal"))

# Plot the distribution of categories
ggplot(pca_data, aes(x = PC1, y = PC2, color = Label)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", title = "PCA Distribution") +
  scale_color_manual(values = c("#339999", "#CC3333")) +
  theme_bw()