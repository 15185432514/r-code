
install.packages("survival")
install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#Download survival information
#xena Website：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Adenocarcinoma%20(LUAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

surv = read.table(file = 'TCGA-LUAD.survival.tsv', sep = '\t', header = TRUE) 

rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
range(surv$OS.time)
#Exclude samples with a follow-up period of less than 30 days
surv <-surv%>%dplyr::select(c('OS','OS.time')) %>% 
  subset(.,OS.time>=30)%>%
  mutate(OS.time=OS.time/30) %>% 
  na.omit()

expr <- read.table("lncRNA_tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
expr <-t(expr)
expr <-as.data.frame(expr)


surv.expr <- cbind(surv,expr)

write.table(surv.expr , file = "surv_expr_01A .txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分

##Survival analysis
setwd("survival")
surv <- read.table("surv_expr_01A .txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

colnames(surv) <- gsub("[-]",".",colnames(surv))
median(surv$KLHDC7B.DT)
surv$group <- ifelse(surv$KLHDC7B.DT > median(surv$KLHDC7B.DT),"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data      = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#Fitting survival curves
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))

#install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE,
           risk.table = F, 
           risk.table.col = "strata",
           palette = "jco", 
           legend.labs = c("Low", "High"), 
           size = 1,
           xlim = c(0,120), 
           break.time.by = 20, 
           legend.title = "",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Months)", 
           
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()
