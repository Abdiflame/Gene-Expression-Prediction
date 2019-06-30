library(randomForest)
#library(gbm)

### Load the dataset
BH_gene <- read.table("/home/lffraga/Table Data/BH_gene21.txt")
BH_passed_snp <- read.table("/home/lffraga/Table Data/all_snps.txt")
snps <- read.table("/home/lffraga//Table Data/all_inf.txt")

snp_log1 <- t(snps)

kfold <- 5 # Dataset K-fold CrossValidation
k <- 5 # Iteration K-repetitions
#folds <- sample(cut(seq(1,nrow(BH_gene)),breaks=kfold,labels=FALSE))

start <- Sys.time()

cv.error.5_rf = rep(0,k)
rsq_rf1<- rep(0,k)

rf_mse_results_log1 <- list()
rf_rsq_results_log1<- list()
ex_data <- list()


for(i in genes68)
{
  BH_passed_snp_tmp<-snp_log1[,na.omit(match(as.character(BH_passed_snp[,i]),colnames(snp_log1)))]
  ex_data<-cbind(BH_passed_snp_tmp,BH_gene[,i])
  
  for(j in 1:k){
    testIndexes <- which(folds==j,arr.ind=TRUE)
    testData <- as.data.frame(ex_data[testIndexes, ])
    trainData <- as.data.frame(ex_data[-testIndexes, ])
    
    ### Cross-validation and model fit
    rt<-dim(ex_data)[2]/3 
    rf.fit<-randomForest(as.matrix(ex_data[-testIndexes,-dim(ex_data)[2]]),as.matrix(ex_data[-testIndexes,dim(ex_data)[2]]), mtry=rt, ntree = 500)
    rf.pred<-predict(rf.fit, as.matrix(ex_data[testIndexes,-dim(ex_data)[2]]))
    
    cv.error.5_rf[j] <- mean((rf.pred-testData[,dim(testData)[2]])^2)
    rsq_rf1[j]<- as.numeric(postResample(rf.pred, ex_data[testIndexes,dim(ex_data)[2]])[2])
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

