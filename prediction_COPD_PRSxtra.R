library(data.table)
library(tidyverse)
library(pROC)

df=fread('~/COPD_PRSxtra_output.csv')
covars=('~/COPD_covars.csv')
df=merge(df,covers,by='ID')

for(smoke in c('all','nonsmoker','smoker')){
  for(anc in c('all','afr','amr','eur')){
    if(anc=='all'){
      dfanc=df
    }
    if(anc!='all'){
      dfanc=df[which(df$ancestry_pred==anc),]
    }
    dfanc=na.omit(dfanc)

    if(smoke=='nonsmoker'){
      dfanc=dfanc[which(dfanc$smokingstatus=='never'),]
    }
    if(smoke=='smoker'){
      dfanc=dfanc[which(dfanc$smokingstatus%in%c('previous','current')),]
    }
    prscol=which(colnames(dfanc)=='prs_std')
       
    roc_obj1 <- roc(dfanc$phenotype, dfanc$prs_std)
    ci1=ci.auc(roc_obj1)                                   
    print(c('PRS AUC alone:', ci.auc(roc_obj1)))
    aucs=rbind(aucs,c('PRS alone:',auc(roc_obj1), ci1[1], ci1[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
    
   
    if(smoke=='all'){
      
      model <- glm(phenotype ~ age + gender , data = dfanc, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dfanc$phenotype, prediction)
      auc_obj1 <- auc(roc_obj)
      ci1 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+smokingstatus , data = dfanc, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dfanc$phenotype, prediction)
      auc_obj2 <- auc(roc_obj)
      ci2 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+smokingstatus+prs_std, data = dfanc, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj3 <- roc(dfanc$phenotype, prediction)
      auc_obj3 <- auc(roc_obj3)
      ci3 <- ci.auc(roc_obj3)
         
      auc_test_result <- roc.test(roc_obj2, roc_obj2)
      print(auc_test_result)
      sigs=rbind(sigs,c('Covariate_vs_PRS:',auc_test_result$p.value,anc,smoke))                                  
      
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
      aucs=rbind(aucs,c('SmokingS',auc_obj2, ci2[1], ci2[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))                            
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
      
    }
    if(smoke!='all'){
      model <- glm(phenotype ~ age + gender , data = dfanc, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dfanc$phenotype, prediction)
      auc_obj1 <- auc(roc_obj)
      ci1 <- ci.auc(roc_obj)
      
      model <- glm(phenotype ~ age + gender+prs_std, data = dfanc, family = "binomial")
      prediction <- predict(model, type = "response")
      roc_obj <- roc(dfanc$phenotype, prediction)
      auc_obj3 <- auc(roc_obj)
      ci3 <- ci.auc(roc_obj)
       
      aucs=rbind(aucs,c('Sex+Age',auc_obj1, ci1[1], ci1[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
      aucs=rbind(aucs,c('PRS',auc_obj3, ci3[1], ci3[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
      aucs=rbind(aucs,c('PRSxtra',auc_obj4, ci4[1], ci4[3],anc,smoke,nrow(dfanc),length(which(dfanc$pheno==1))))
      
    }              
  }
}
write.csv(aucs,'~/AUC_COPD_PRSxtra.csv')
write.csv(sig,'~/compareP_COPD_PRSxtra.csv')
