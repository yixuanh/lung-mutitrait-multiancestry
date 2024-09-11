library(data.table)
library(tidyverse)
library(bigrquery)
library(glmnet)
library(pROC)
library(gridExtra)
library(ggplot2)


df=fread('path/to/indivdual/copd/data')
prscxs=fread('PRS/COPD/compiled_PRSCSX.csv')
head(prscxs)
df=merge(df,prscxs,by='person_id',all.x=T)

train=fread('/ID_train.csv')
test=fread('/ID_test.csv')

df$age=2023-as.numeric(sapply(dfpop$date_of_birth,function(x) substr(x,1,4)))

#run glmnet
dff=df[-which(df$SERPINA==1),]
pi=1
plot_list <- list()
aucs=c()
ors=c()
sigs=c()
for(smoke in c('all','nonsmoker','smoker')){
  for(anc in c('all','afr','amr','eur')){
    if(anc=='all'){
      dfanc=dff
    }
    if(anc!='all'){
      dfanc=dff[which(dff$ancestry_pred==anc),]
    }
    dfanc=na.omit(dfanc)

    if(smoke=='nonsmoker'){
      dfanc=dfanc[which(dfanc$smokingstatus=='never'),]
    }
    if(smoke=='smoker'){
      dfanc=dfanc[which(dfanc$smokingstatus%in%c('previous','current')),]
    }
    prscol=which(colnames(dfanc)==paste('COPD_',anc,'_prscs',sep=''))
    if(anc=='all'){
      prscol=which(colnames(dfanc)==paste('COPD_eur_prscs',sep='')) 
    }
    
    dfanc[,c(19:ncol(dfanc))]=scale(dfanc[,c(19:ncol(dfanc))]) #scale scores
    
    #make training data
    x_train <- as.matrix(train)
    y_train <- train$phenotype
    
    #make test data
    x_test <- as.matrix(test)
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
    write.csv(em,paste('COPD_analysis/modelweights_',smoke,'_',anc,'.csv',sep=''))
    
    #predict using EN model
    predictions <- predict(enet_model, newx = x_test, type = "response")
    
    dftest=data.frame(cbind(test,predictions))
    colnames(dftest)[ncol(dftest)]='prsxtra'

    #standardize
    dftest$prs_std <- dftest$prs
    dftest$prsxtra_std <- (dftest$prsxtra - mean(dftest$prsxtra, na.rm = TRUE)) / sd(dftest$prsxtra, na.rm = TRUE)
    
    write.csv(dftest,paste('COPD_',smoke,'_',anc,'.csv',sep=''))
    
    roc_obj1 <- roc(dftest$phenotype, dftest$prs_std)
    ci1=ci.auc(roc_obj1)                                   
    print(c('PRS AUC alone:', ci.auc(roc_obj1)))
    aucs=rbind(aucs,c('PRS alone:',auc(roc_obj1), ci1[1], ci1[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
    
    roc_obj2 <- roc(dftest$phenotype, dftest$prsxtra_std)
    ci2=ci.auc(roc_obj2)                                   
    print(c('PRSxtra AUC alone:',ci.auc(roc_obj2)))
    aucs=rbind(aucs,c('PRSxtra alone:',auc(roc_obj2), ci2[1], ci2[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
    
    auc_test_result <- roc.test(roc_obj1, roc_obj2)
    print(auc_test_result)
    sigs=rbind(sigs,c('PRS alone comp:',auc_test_result$p.value,anc,smoke))                                  
    print(dim(dftest))
    print(table(dftest$pheno))                        
    
    if(smoke=='all'){
      
      model <- glm(phenotype ~ age + gender , data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj1 <- auc(roc_obj)
      ci1 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+smokingstatus , data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj2 <- auc(roc_obj)
      ci2 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+smokingstatus+prs_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj3 <- roc(dftest$phenotype, prediction)
      auc_obj3 <- auc(roc_obj3)
      ci3 <- ci.auc(roc_obj3)
      
      model <- glm(phenotype ~ age + gender+smokingstatus+prsxtra_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj4 <- roc(dftest$phenotype, prediction)
      auc_obj4 <- auc(roc_obj4)
      ci4 <- ci.auc(roc_obj4)
      
      auc_test_result <- roc.test(roc_obj3, roc_obj4)
      print(auc_test_result)
      sigs=rbind(sigs,c('PRS alone joint:',auc_test_result$p.value,anc,smoke))                                  
      
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('SmokingS',auc_obj2, ci2[1], ci2[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))                            
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRSxtra',auc_obj4, ci4[1], ci4[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      
    }
    if(smoke!='all'){
      model <- glm(phenotype ~ age + gender , data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj1 <- auc(roc_obj)
      ci1 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+prs_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj3 <- auc(roc_obj)
      ci3 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+prsxtra_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj4 <- auc(roc_obj)
      ci4 <- ci.auc(roc_obj)
      
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRSxtra',auc_obj4, ci4[1], ci4[3],anc,smoke,nrow(dftest),length(which(dftest$pheno==1))))
      
    }
    
    
    
    temp=data.frame(as.matrix(coef(enet_model)))
    temp$v2=row.names(temp)
    print(temp[order(temp[,1],decreasing=T),])                             
  }
}
print(data.frame(aucs))
