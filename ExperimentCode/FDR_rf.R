library(randomForest)

### Load the dataset
BH_gene <- read.table("/home/lffraga/Table Data/BH_gene6.txt")
BH_passed_snp <- read.table("/home/lffraga/Table Data/all_snps.txt")
snps <- read.table("/home/lffraga/Table Data/all_inf.txt")

snp_log1 <- t(snps)
kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
#folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold, labels=FALSE))

start <- Sys.time()

cv.error.5_rf = rep(0,k)
rsq_rf1<- rep(0,k)
rf_mse_results_log1 <- list()
rf_rsq_results_log1<- list()
ex_data <- list()

for(i in genes68)
{
  BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_passed_snp),colnames(snp_log1)))]
  ex_data_refresh <- cbind(BH_passed_snp_tmp,BH_gene[,i])
  for(j in 1:k){
    ex_data <- ex_data_refresh
    testIndexes <- which(folds==j,arr.ind=TRUE)

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
    
    if(dim(BH_snps_aux)[1]>1){### Update ex_data
      BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_snps_aux[,1]),colnames(snp_log1)))]
      ex_data <- cbind(BH_passed_snp_tmp,BH_gene[,i])
      
      ### Cross-validation and model fit
      rt <- dim(ex_data)[2]/3 
      rf.fit <- randomForest(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]),as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), maxnodes = 5, mtry=rt, ntree = 200)
      rf.pred <- predict(rf.fit, as.matrix(ex_data[testIndexes,-dim(ex_data)[2]]))
      
      cv.error.5_rf[j] <- mean((rf.pred - ex_data[testIndexes,dim(ex_data)[2]])^2)
      rsq_rf1[j] <- as.numeric(postResample(rf.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])
    }
    
    else if(dim(BH_snps_aux)[1]<=1){
      print("deu ruim aqui")
      cv.error.5_rf <- 'NA'
      rsq_rf1 <- 'NA'
      break
    }
    
  }
  
  print("Resultados:")
  print(cv.error.5_rf)
  print(rsq_rf1)
  
  rf_mse_results_log1 <- cbind(rf_mse_results_log1, cv.error.5_rf)
  rf_rsq_results_log1 <- cbind(rf_rsq_results_log1, rsq_rf1)
}

l_rf_mse_log1 <-list()
l_rf_rsq_log1 <-list()

for(i in 1:dim(rf_mse_results_log1)[2])
{
  l_rf_mse_log1[i] <-mean(as.numeric(rf_mse_results_log1[,i]))
  l_rf_rsq_log1[i] <-mean(as.numeric(rf_rsq_results_log1[,i]))
  print(i)
}

end <- Sys.time()
end - start

### Check Results
mean(as.numeric(l_rf_mse_log1), na.rm = TRUE)
mean(as.numeric(l_rf_rsq_log1), na.rm = TRUE)
sd(as.numeric(l_rf_mse_log1), na.rm = TRUE)
sd(as.numeric(l_rf_rsq_log1), na.rm = TRUE)


