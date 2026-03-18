
#install.packages("readr")
library(readr)
library(tidyverse)

mydata=read.table("2hub.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
mydata=t(mydata)
mydata=as.data.frame(mydata)

mydata<- rownames_to_column(mydata,var='sample')
x <- ifelse(str_detect(mydata$sample,"01A"),"1","0")
mydata$class <- x

rownames(mydata) <- NULL
mydata <-mydata%>% column_to_rownames("sample")

mydata<-na.omit(mydata)

View(mydata)

names(mydata)

str(mydata)

mydata$class <- ifelse(mydata$class == "1",1,0)

attach(mydata)

str(mydata)   


colnames(mydata)
formula <- as.formula(class==1~ LINC00857+KLHDC7B.DT )
modelA<-glm(formula,data =mydata,family = binomial(logit),control=list(maxit=100))

mydata$predmodelA<- predict(newdata=mydata,modelA,"response")

View(mydata)
#install.packages("pROC")

library(pROC)

table(mydata$class)

devmodelA <- roc(class~predmodelA, data = mydata,smooth=F)

devmodelA

round(auc(devmodelA),3)
round(ci(auc(devmodelA)),3)

plot(devmodelA, 
     print.auc=F, 
     print.thres=F,
     legacy.axes=TRUE,
     main = "", 
     col= "red",
     print.thres.col="blue",
     identity.col="blue",
     identity.lty=1,
     identity.lwd=1)
roc2 <- roc(mydata$class,mydata$KLHDC7B.DT)
roc3 <- roc(mydata$class,mydata$LINC00857)
plot(roc2, add=TRUE, col="green")
plot(roc3, add=TRUE, col="blue")

auc(roc2)
ci(auc(roc2))
auc(roc3)
ci(auc(roc3))

legend("bottomright",legend=c("7+K(AUC=0.982)",
                              "DT(AUC=0.925)",
                              "857（AUC=0.944)"),
       col=c("red","green","blue"),
       lty=1,
       lwd=3)

