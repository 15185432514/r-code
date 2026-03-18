library(e1071)
library(kernlab)
library(caret)

set.seed(123)
mydata <- read.csv("2UP.csv",row.names = 1)
X <- mydata[,1:537] 
Y <- as.numeric(as.factor(mydata$group))

#SVM-RFE
ctrl <- rfeControl(functions=caretFuncs, method="cv", number=10)
grid <- expand.grid(.sigma=c(0.1, 0.5, 1), .C=seq(0.5, 1, length=10))

Profile <- rfe(X, Y,
               sizes = c(2,4,6,8, seq(10,40,by=3)),
               rfeControl = ctrl,
               method = "svmRadial",
               tuneGrid = grid)

pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
