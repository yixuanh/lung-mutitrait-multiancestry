library(data.table)
library(tidyverse)
library(bigrquery)
library(glmnet)

df=fread('~/phenotype.csv')
prslist=fread('~/comipledPRS/Asthma/compiled_PRSCSX.csv')
df=merge(df,prslist,by='ID',all.x=T)
df[,c(2:ncol(df))]=scale(df[,c(2:ncol(df))]) #scale all PRS scores

train=fread('/ID_train.csv')
test=fread('/ID_test.csv')

train=na.omit(dfanc[which(dfanc$ID%in%train$x),])
test=na.omit(dfanc[which(dfanc$ID%in%test$x),])

x_train <- as.matrix(train[, grep("^PRS", colnames(train))])
y_train <- train$phenotype

x_test <- as.matrix(train[, grep("^PRS", colnames(x_test))])
y_test <- test$phenotype

num_folds <- 10

#creat cross-validated glmnet model
cv_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0, nfolds = num_folds)

#get the best lambda value selected by cross-validation (minimum)
best_lambda <- cv_model$lambda.min

#fit the final glmnet model using the best lambda
enet_model <- glmnet(x_train, y_train, family = "binomial", alpha = 0, lambda = best_lambda)

#get weights
em=data.frame(as.matrix(enet_model$beta))
write.csv(em,'modelweights.csv')

#predict using EN model
predictions <- predict(enet_model, newx = x_test, type = "response")

dftest=data.frame(cbind(test,predictions))
colnames(dftest)[ncol(dftest)]='prs'

#standardize scores
dftest$prs_std <- dftest$prs
write.csv(dftest,('scores_output.csv',row.names=F)

