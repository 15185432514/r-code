
# Check whether the glmnet package is already installed
if (!require("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}

library(glmnet)

mydata <- read.csv("2UP.csv",row.names = 1,header = TRUE,check.names = F)

colnames(mydata[,1:10])

mydata$group <- ifelse(mydata$group == "1",1,0) 

y <- as.matrix(mydata[,538])  
x <- as.matrix(mydata[, 1:537]) 


# Lasso regression for variable selection 
set.seed(12345)
lasso_model <- glmnet(x, 
                      y, 
                      family = "binomial",
                      alpha = 1) 
print(lasso_model) 

plot(lasso_model,
     xvar = "lambda",
     label = F)

coef_lasso <- coef(lasso_model,
                   s = 0.071930) 
coef_lasso


#Cross-validation for selecting the appropriate Lambda
cv_model <- cv.glmnet(x, y, family = "binomial",alpha = 1,nfolds = 10)

plot(cv_model)

# Based on the cross-validation results, select the lambda value.
lambda_min <- cv_model$lambda.min
lambda_min
lambda_1se <- cv_model$lambda.1se
lambda_1se

coef_cv <- coef(lasso_model, s = lambda_1se)
coef_cv 

# Calculate the odds ratio based on the regression coefficient
exp(coef_cv)

coef_cv <- as.matrix(coef_cv)
coef_cv <- data.frame(coef_cv)

coef_cv$OR <- exp(coef_cv$s1)# Calculate the OR value for each variable
nonzero_vars <- rownames(coef_cv[coef_cv$OR != 1, ])
nonzero_vars <- nonzero_vars [2:26]

lasso_data <- mydata[,nonzero_vars]

write.csv(lasso_data,file = "2UPlasso_data.csv",quote = F)

