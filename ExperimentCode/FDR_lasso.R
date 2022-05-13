#library(caret)
#library(glmnet)

### Load the dataset
BH_gene <- read.table("geneExpPrediction/ResultsData/BH_gene6.txt")
snp = read.table("geneExpPrediction/ResultsData/all_snps.txt")
BH_passed_snp <- read.table("snp_inf22.txt")
snps <- read.table("geneExpPrediction/ResultsData/snp_inf6.txt")
genes68 <- read.table("geneExpPrediction/ResultsData/linear_mse_68.txt")

snp_log1 <- t(snps)
kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold, labels=FALSE))

start <- Sys.time()

cv.error.5_lasso <- rep(0,k)
rsq_lasso1 <- rep(0,k)

lasso_mse_results_log1 <- list()
lasso_rsq_results_log1<- list()

ex_data <- list()
lasso_lambda <- list()
lasso_df <- list()

for(i in 1:dim(BH_gene)[2])#1:dim(BH_gene)[2]
{
  BH_passed_snp_tmp <- snp_log1
  ex_data_refresh <- cbind(BH_passed_snp_tmp,BH_gene[,i])
  
  for(j in 1:k){
    ex_data <- ex_data_refresh
    
    testIndexes <- which(folds==j,arr.ind=TRUE)
    testData <- as.data.frame(ex_data[testIndexes, ])
    trainData <- as.data.frame(ex_data[-testIndexes, ])
    
    ### Make the FDR correction here?
    l_p <- list()
    end <- dim(ex_data)[2]-1
    
    print("START PVALUE")
    print(i)
    
    for(p in 1:end)
    {
      l_p[p] <- anova(lm(ex_data[-testIndexes,end+1] ~ ex_data[-testIndexes,p]))[5][1,]
    }
    
    cis_snps_pvalue <- as.numeric(l_p)
    auxVec <- p.adjust(cis_snps_pvalue, method = "BH")           # Getting the SNPs adjusted pvalue
    snp_names_aux <- c(as.character(colnames(ex_data[,1:end])))  # Getting the SNPs names
    df_aux <- data.frame(snp_names_aux, auxVec)
    names(df_aux) <- c("SNP", "AdjPvalue")
    BH_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.05)     # Applying the p-value threshold
    
    if(dim(BH_snps_aux)[1]>1){
      ### Update ex_data
      BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_snps_aux[,1]),colnames(snp_log1)))]
      ex_data <- cbind(BH_passed_snp_tmp,BH_gene[,i])
      
      ### Cross-validation and model fit
      cv.out_lasso <- cv.glmnet(as.matrix(ex_data [-testIndexes,-dim(ex_data)[2]]),
                                as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), alpha = 1,
                                nfolds = k)
      
      lasso.fit <- glmnet(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]),
                          as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), 
                          alpha = 1, lambda = cv.out_lasso$lambda.min)
      
      #plot.glmnet(lasso.fit)
      lasso.pred<-predict(lasso.fit, as.matrix(ex_data[testIndexes,-dim(ex_data)[2]]))
      
      cv.error.5_lasso[j] <- mean((testData[,dim(testData)[2]]-lasso.pred)^2)
      rsq_lasso1[j] <- as.numeric(postResample(lasso.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])
    }
    
    else if(dim(BH_snps_aux)[1]<=1){
      print("deu ruim aqui")
      cv.error.5_lasso <- 'NA'
      rsq_lasso1 <- 'NA'
      break
    }
    
  }
  print("Resultados:")
  print(cv.error.5_lasso)
  print(rsq_lasso1)
  print("Parameter")
  if(cv.error.5_lasso != 'NA'){
    print(lasso.fit$lambda)
    lasso_lambda[i] <- lasso.fit$lambda
    print("Number of Predictors")
    print(lasso.fit$df)
    lasso_df[i] <- lasso.fit$df
  }
  
  lasso_mse_results_log1 <- cbind(lasso_mse_results_log1, cv.error.5_lasso)
  lasso_rsq_results_log1 <- cbind(lasso_rsq_results_log1, rsq_lasso1)

}

l_lasso_mse_log1 <- list()
l_lasso_rsq_log1 <- list()

for(i in 1:dim(lasso_mse_results_log1)[2])
{
  l_lasso_mse_log1[i] <- mean(as.numeric(lasso_mse_results_log1[,i]))
  l_lasso_rsq_log1[i] <- mean(as.numeric(lasso_rsq_results_log1[,i]))
  
  print(i)
}

end <- Sys.time()
end-start

### Check Results
mean(as.numeric(l_lasso_mse_log1), na.rm = TRUE)
mean(as.numeric(l_lasso_rsq_log1), na.rm = TRUE) 
sd(as.numeric(l_lasso_mse_log1), na.rm = TRUE)
sd(as.numeric(l_lasso_rsq_log1), na.rm = TRUE) 
