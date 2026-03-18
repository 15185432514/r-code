
if (!require("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}

library(randomForest)
library(tidyverse)

mydata<- read.table("2UP.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

mydata$group<- factor(mydata$group,level=c(0,1),labels=c("Normal","Tumor"))
set.seed(155)
colnames(mydata) <- gsub("[-]",".",colnames(mydata))
rf_model <- randomForest(group~.,data=mydata,ntree = 1000, 
                         importance = TRUE)
rf_model
plot(rf_model)

which.min(rf_model$err.rate[,1])
#Refit the model using the tree with the lowest error
rf_model2 <- randomForest(group~.,data=mydata,ntree=49,importance = TRUE)
rf_model2
#variable importance scores
varImpPlot(rf_model2)
varImpPlot(rf_model2, 
           type=2, 
           main="Variable Importance",
           n.var = 25,
           scale = T,
           cex = 0.8
)

importance <- importance(rf_model2)
importance <- as.data.frame(importance)

top25_MeanDecreaseGini <- head(rownames(importance[order(-importance$`MeanDecreaseGini`),]), 25)
top25_MeanDecreaseGini
write.csv(top25_MeanDecreaseGini,"top25_MeanDecreaseGini_data.csv")
