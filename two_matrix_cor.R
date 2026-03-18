rm(list=ls())
library(limma)
corFilter=0.3           
pvalueFilter=0.05     


lncRNA=read.table("2hub.txt", header=T, sep="\t", check.names=F)
rownames(lncRNA)=lncRNA[,1]
lncRNA=lncRNA[,-1]
lncRNA=as.matrix(lncRNA)

mRNA=read.table("mRNA_tpms_log2.txt", header=T, sep="\t", check.names=F)
mRNA=as.matrix(mRNA)
rownames(mRNA)=mRNA[,1]
exp1=mRNA[,2:ncol(mRNA)]
dimnames1=list(rownames(exp1),colnames(exp1))
mRNA=matrix(as.numeric(as.matrix(exp1)), nrow=nrow(exp1), dimnames=dimnames1)
mRNA=avereps(mRNA)
mRNA=mRNA[rowMeans(mRNA)>0.1,]

#Correlation test
outTab=data.frame()
for(i in row.names(lncRNA)){
  if(sd(lncRNA[i,])>0.1){
    if(1){
      for(j in row.names(mRNA)){
        x=as.numeric(lncRNA[i,])
        y=as.numeric(mRNA[j,])
        corT=cor.test(x,y)
        cor=corT$estimate
        pvalue=corT$p.value
        if((cor>corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(mRNA=j,lncRNA=i,cor,pvalue,Regulation="postive"))
        }
        if((cor< -corFilter) & (pvalue<pvalueFilter)){
          outTab=rbind(outTab,cbind(mRNA=j,lncRNA=i,cor,pvalue,Regulation="negative"))
        }
      }
    }
  }
}

write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)
