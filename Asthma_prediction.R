library(data.table)
library(tidyverse)
library(bigrquery)
library(glmnet)
library(pROC)
library(gridExtra)
library(ggplot2)

df=fread('path/to/indivdual/asthma/data')
prscxs=fread('PRS/Asthma/compiled_PRSCSX.csv')
head(prscxs)
df=merge(df,prscxs,by='person_id',all.x=T)

train=fread('/ID_train.csv')
test=fread('/ID_test.csv')

df$age=2023-as.numeric(sapply(dfpop$date_of_birth,function(x) substr(x,1,4)))

#run glmnet
dff=df
pi=1
plot_list <- list()
aucs=c()
ors=c()
sigs=c()

for(famhx in c('all','yes','no')){
  for(anc in c('all','afr','amr','eur','eas')){
    if(anc=='all'){
      dfanc=dff
    }
    if(anc!='all'){
      dfanc=dff[which(dff$ancestry_pred==anc),]
    }
    dfanc=na.omit(dfanc)
    print(anc)
    print(famhx)
    if(famhx!='all'){
      dfanc=dfanc[which(dfanc$asthmafx==famhx),]
    }
    
    prscol=which(colnames(dfanc)==paste('Asthma',anc,'_prscs',sep=''))
    if(anc=='all'){
      prscol=which(colnames(dfanc)==paste('Asthma_eur_prscs',sep=''))
    }

    dfanc[,c(19:ncol(dfanc))]=scale(dfanc[,c(19:ncol(dfanc))]) #scale scores
    
    train=na.omit(dfanc[which(dfanc$person_id%in%train$x),])
    test=na.omit(dfanc[which(dfanc$person_id%in%test$x),])

    x_train <- as.matrix(train)
    y_train <- train$phenotype
    
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
    write.csv(em,paste('Asthma_analysis/modelweights_',smoke,'_',anc,'.csv',sep=''))
    
    #predict using EN model
    predictions <- predict(enet_model, newx = x_test, type = "response")
    
    dftest=data.frame(cbind(test,predictions))
    colnames(dftest)[ncol(dftest)]='prsxtra'

    #standardize
    dftest$prs_std <- dftest$prs
    dftest$prsxtra_std <- (dftest$prsxtra - mean(dftest$prsxtra, na.rm = TRUE)) / sd(dftest$prsxtra, na.rm = TRUE)

    write.csv(dftest,paste('Asthma_',famhx,'_',anc,'.csv',sep=''))
    
    roc_obj1 <- roc(dftest$phenotype, dftest$prs_std)
    ci1=ci.auc(roc_obj1)                                   
    print(c('PRS AUC alone:', ci.auc(roc_obj1)))
    aucs=rbind(aucs,c('PRS alone:',auc(roc_obj1), ci1[1], ci1[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
    
    roc_obj2 <- roc(dftest$phenotype, dftest$prsxtra_std)
    ci2=ci.auc(roc_obj2)                                   
    print(c('PRSxtra AUC alone:',ci.auc(roc_obj2)))
    aucs=rbind(aucs,c('PRSxtra alone:',auc(roc_obj2), ci2[1], ci2[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
    
    auc_test_result <- roc.test(roc_obj1, roc_obj2)
    print(auc_test_result)
    sigs=rbind(sigs,c('PRS alone:',auc_test_result$p.value,anc,famhx))                                  
    print(dim(dftest))
    print(table(dftest$pheno))        
    
    if(famhx=='all'){
      model <- glm(phenotype ~ age + gender , data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj1 <- auc(roc_obj)
      ci1 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+asthmafx , data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dftest$phenotype, prediction)
      auc_obj2 <- auc(roc_obj)
      ci2 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+asthmafx+prs_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj3 <- roc(dftest$phenotype, prediction)
      auc_obj3 <- auc(roc_obj3)
      ci3 <- ci.auc(roc_obj3)
      
      model <- glm(phenotype ~ age + gender+asthmafx+prsxtra_std, data = dftest, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj4 <- roc(dftest$phenotype, prediction)
      auc_obj4 <- auc(roc_obj4)
      ci4 <- ci.auc(roc_obj4)
      
      auc_test_result <- roc.test(roc_obj3, roc_obj4)
      print(auc_test_result)
      sigs=rbind(sigs,c('PRS joint:',auc_test_result$p.value,anc,famhx))                                  
      
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('Asthma_FX',auc_obj2, ci2[1], ci2[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))                            
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRSxtra',auc_obj4, ci4[1], ci4[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      
    }
    if(famhx!='all'){
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
      
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      aucs=rbind(aucs,c('PRSxtra',auc_obj4, ci4[1], ci4[3],anc,famhx,nrow(dftest),length(which(dftest$pheno==1))))
      
    }
  
    temp=data.frame(as.matrix(coef(enet_model)))
    temp$v2=row.names(temp)
    print(temp[order(temp[,1],decreasing=T),])                             
  }
}
print(aucs)


