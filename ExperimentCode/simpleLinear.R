### Load the dataset
BH_gene <- read.table("../Dataset/BH_gene6.txt")
BH_passed_snp <- read.table("../Dataset/BH_passed_snps6.txt")
snps <- read.table("../Dataset/snp_inf6.txt")

snp_log1 <- t(snps)

kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
#folds <- scan("/home/lffraga/Table Data/folds.txt")
#folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold,labels=FALSE))

start <- Sys.time()

cv.error.5_sl <- rep(0,k)
rsq_sl1 <- rep(0,k)

sl_mse_results_log1 <- list()
sl_rsq_results_log1 <- list()

lm_coef <- list()
ex_data <- list()

for(i in 1:length(BH_gene[,3])){
  print(i)
}

BH_passed_snp <- BH_passed_snp[,3]

for(i in genes68){}

for(i in 1:length(BH_gene[,3])){
  BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_passed_snp[i]),colnames(snp_log1)))]
  ex_data_refresh <- cbind(BH_passed_snp_tmp,BH_gene[i,3])
  

  #BH_passed_snp_tmp <- snp_log1[,na.omit(match(as.character(BH_passed_snp[,i]),colnames(snp_log1)))]
  #ex_data_refresh <- cbind(BH_passed_snp_tmp,BH_gene[,i])
  
  for(j in 1:k){
    ex_data <- ex_data_refresh
    
    testIndexes <- which(folds==j,arr.ind=TRUE)
    testData <- as.data.frame(ex_data[testIndexes, ])
    trainData <- as.data.frame(ex_data[-testIndexes, ])
    
    #sl_mse<-list()
    #sl_rsq<-list()
    
    snp_pvalue<-rep(0,dim(ex_data)[2]-1)
    
    for(sim in 1:(dim(ex_data)[2]-1))
    {
      gd<-c(names(trainData)[sim],names(trainData)[dim(trainData)[2]])
      tm<-as.data.frame(cbind(trainData[,sim], trainData[,dim(trainData)[2]]))
      colnames(tm) <- gd
      
      lm.fit <-lm(BH_gene[, i]~ BH_gene[,266], data = tm) 
      snp_pvalue[sim]<-anova(lm.fit)[5][1,1]
    }
    
    gd <- c(names(trainData)[match(min(as.numeric(snp_pvalue)),snp_pvalue)],names(trainData)[dim(trainData)[2]])
    tm <- as.data.frame(cbind(trainData[,match(min(as.numeric(snp_pvalue)),snp_pvalue)], trainData[,dim(trainData)[2]]))
    colnames(tm) <- gd
    
    lm.fit <- lm(BH_gene[, i]~ BH_gene[,266], data = tm)
    lm.pred <- predict(lm.fit,testData)
    
    cv.error.5_sl[j] <- mean((testData[,dim(testData)[2]]-predict(lm.fit,testData))^2)
    rsq_sl1[j] <- as.numeric(postResample(lm.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])
  }
  
  print("Resultados:")
  print(cv.error.5_sl)
  print(rsq_sl1)
  
  print("Parameter")
  print(lm.fit$coefficients)
  lm_coef[i] <- lm.fit$coefficients
  
  sl_mse_results_log1 <- cbind(sl_mse_results_log1, cv.error.5_sl)
  sl_rsq_results_log1 <- cbind(sl_rsq_results_log1,rsq_sl1)
  print(i)
}

l_sl_mse_log1 <- list()
l_sl_rsq_log1 <- list()

for(i in 1:dim(sl_mse_results_log1)[2])
{
  l_sl_mse_log1[i]<-mean(as.numeric(sl_mse_results_log1[,i]))
  l_sl_rsq_log1[i]<-mean(as.numeric(sl_rsq_results_log1[,i]))
  print(i)
}

end <- Sys.time()
end-start


### Check Results
mean(as.numeric(l_sl_mse_log1))
mean(as.numeric(l_sl_rsq_log1), na.rm = TRUE)
sd(as.numeric(l_sl_mse_log1))
sd(as.numeric(l_sl_rsq_log1))

write.table(l_sl_mse_log1, "/linear_mse_68.txt")
write.table(l_sl_rsq_log1, "/home/lffraga/Results_Data/linear_rsq_68.txt")
