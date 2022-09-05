library(caret)
#library(MASS)
#library(glmnet)
#library(randomForest)
#library(gbm)

### Load the dataset
BH_gene <- read.table("../Dataset/BH_gene6.txt")
BH_passed_snp <- read.table('../Dataset/BH_passed_snps6_original.txt')
snps <- read.table('../Dataset/snp_inf6.txt')

snp_log1 <- t(snps)

kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold,labels=FALSE))

start <- Sys.time()

cv.error.5_lasso = rep(0,k)
rsq_lasso1 <- rep(0,k)

lasso_mse_results_log1 <- list()
lasso_rsq_results_log1 <- list()

ex_data <- list()
lasso_lambda <- list()
lasso_df <- list()

for(i in dim(BH_gene)[2])
{
  BH_passed_snp_tmp<-snp_log1[,na.omit(match(as.character(BH_passed_snp[,i]),colnames(snp_log1)))]
  ex_data<-cbind(BH_passed_snp_tmp,BH_gene[,i])
  
  for(j in 1:k){
    testIndexes <- which(folds==j,arr.ind=TRUE)
    
    ### Cross-validation and model fit
    cv.out_lasso <- cv.glmnet(as.matrix(ex_data [-testIndexes,-dim(ex_data)[2]]),
                              as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), alpha = 1,
                              nfolds = k)

    lasso.fit <- glmnet(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]),
                        as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), 
                        alpha = 1, lambda = cv.out_lasso$lambda.min)
    
    #plot.glmnet(lasso.fit)
    lasso.pred <- predict(lasso.fit, as.matrix(ex_data[testIndexes,-dim(ex_data)[2]]))
    
    cv.error.5_lasso[j] <- mean((ex_data[testIndexes,dim(ex_data)[2]] - lasso.pred)^2)
    rsq_lasso1[j]<- as.numeric(postResample(lasso.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])

  }
  
  print("Resultados:")
  print(cv.error.5_lasso)
  print(rsq_lasso1)
  print("Parameter")
  print(lasso.fit$lambda)
  lasso_lambda[i] <- lasso.fit$lambda
  print("Number of Predictors")
  print(lasso.fit$df)
  lasso_df[i] <- lasso.fit$df
  
  lasso_mse_results_log1 <- cbind(lasso_mse_results_log1, cv.error.5_lasso)
  lasso_rsq_results_log1 <- cbind(lasso_rsq_results_log1, rsq_lasso1)
}

l_lasso_mse_log1 <-list()
l_lasso_rsq_log1 <-list()

for(i in 1:dim(lasso_mse_results_log1)[2])
{
  l_lasso_mse_log1[i] <- mean(as.numeric(lasso_mse_results_log1[,i]), na.rm = TRUE)
  l_lasso_rsq_log1[i] <- mean(as.numeric(lasso_rsq_results_log1[,i]), na.rm = TRUE)
  print(i)
}

end <- Sys.time()
end-start

### Check Results
mean(as.numeric(l_lasso_mse_log1), na.rm = TRUE) 
mean(as.numeric(l_lasso_rsq_log1), na.rm = TRUE) 
sd(as.numeric(l_lasso_mse_log1), na.rm = TRUE) 
sd(as.numeric(l_lasso_rsq_log1), na.rm = TRUE) 
