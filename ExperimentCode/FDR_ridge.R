install.packages("glmnet")
library(glmnet)
library(caret)
#library(gbm)

BH_gene <- read.table("/home/lffraga/Table Data/BH_gene21.txt")
BH_passed_snp <- read.table("/home/lffraga/Table Data/BH_passed_snps21b.txt")
snps <- read.table("/home/lffraga/Table Data/snp_inf21.txt")

snp_log1 <- t(snps)

kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
#folds <- scan("/home/lffraga/Table Data/folds.txt")
#folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold,labels=FALSE))

start <- Sys.time()

cv.error.5_ridge <- rep(0,k)
rsq_ridge1 <- rep(0,k)

ridge_mse_results_log1 <- list()
ridge_rsq_results_log1 <- list()

ex_data <- list()
ridge_lambda <- list()
ridge_df <- list()

for(i in 244)
{
  BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_passed_snp[,i]),colnames(snp_log1)))]
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
    #dim(df_aux)
    BH_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.05)     # Applying the p-value threshold
    dim(BH_snps_aux)
    
    if(dim(BH_snps_aux)[1]>1){
      ### Update ex_data
      BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_snps_aux[,1]),colnames(snp_log1)))]
      ex_data <- cbind(BH_passed_snp_tmp,BH_gene[,i])
      
      ### Cross-validation and model fit
      cv.out_ridge <- cv.glmnet(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]),
                                as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), 
                                alpha = 0, nfolds = k)
      
      ridge.fit <- glmnet(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]), 
                          as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]),
                          alpha = 0, lambda = cv.out_ridge$lambda.min)
      
      ridge.pred <- predict(ridge.fit,  as.matrix(ex_data[testIndexes,-dim(ex_data)[2]]))
      
      cv.error.5_ridge[j] <- mean((testData[,dim(testData)[2]]-ridge.pred)^2)
      rsq_ridge1[j] <- as.numeric(postResample(ridge.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])
    }
    else if(dim(BH_snps_aux)[1]<=1){
      print("deu ruim aqui")
      cv.error.5_ridge <- 'NA'
      rsq_ridge1 <- 'NA'
      break
    }
  }
  print("Resultados:")
  print(cv.error.5_ridge)
  print(rsq_ridge1)
  print("Parameter")
  print(ridge.fit$lambda)
  ridge_lambda[i] <- ridge.fit$lambda
  print("Number of Predictors")
  print(ridge.fit$df)
  ridge_df[i] <- ridge.fit$df
  
  ridge_mse_results_log1 <- cbind(ridge_mse_results_log1, cv.error.5_ridge)
  ridge_rsq_results_log1 <- cbind(ridge_rsq_results_log1, rsq_ridge1)
}

l_ridge_mse_log1 <- list()
l_ridge_rsq_log1 <- list()

for(i in 1:dim(ridge_mse_results_log1)[2])
{
  l_ridge_mse_log1[i] <- mean(as.numeric(ridge_mse_results_log1[,i]))
  l_ridge_rsq_log1[i] <- mean(as.numeric(ridge_rsq_results_log1[,i]))
  print(i)
}

end <- Sys.time()
end-start

### Check Results
mean(as.numeric(l_ridge_mse_log1), na.rm = TRUE) 
mean(as.numeric(l_ridge_rsq_log1), na.rm = TRUE) 
sd(as.numeric(l_ridge_mse_log1), na.rm = TRUE) 
sd(as.numeric(l_ridge_rsq_log1), na.rm = TRUE) 
