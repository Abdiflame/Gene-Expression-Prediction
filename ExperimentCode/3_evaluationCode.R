install.packages("qqman")
library(qqman)

### MSE
sim_1_mse <- read.table("geneExpPrediction/ResultsData/linear_mse_68.txt")
sim_2a_mse <- read.table("geneExpPrediction/ResultsData/scenario2a_total_mse.txt")
sim_2b_mse <- read.table("geneExpPrediction/ResultsData/scenario2b_titak_mse.txt")
# sim_3a_mse <- read.table("/home/dmk/geneExpPrediction/ResultsData/scenario3a_total_mse.txt")
sim_3b_mse <- read.table("geneExpPrediction/ResultsData/scenario3b_total_mse.txt")

colnames(sim_2a_mse) <- c("Ridge", "Lasso", "ENET", "RF")
colnames(sim_2b_mse) <- c("Ridge", "Lasso", "ENET", "RF")
# colnames(sim_3a_mse) <- c("Ridge", "Lasso", "ENET", "RF")
colnames(sim_3b_mse) <- c("Ridge", "Lasso", "ENET", "RF")

### RSQ
sim_1_rsq <- read.table("geneExpPrediction/ResultsData/linear_rsq_68.txt")
sim_2a_rsq <- read.table("geneExpPrediction/ResultsData/scenario2a_total_rsq.txt")
sim_2b_rsq <- read.table("geneExpPrediction/ResultsData/scenario2b_total_rsq.txt")
sim_3a_rsq <- read.table("geneExpPrediction/ResultsData/scenario3a_total_rsq.txt")
sim_3b_rsq <- read.table("geneExpPrediction/ResultsData/scenario3b_total_rsq.txt")

colnames(sim_2a_rsq) <- c("Ridge", "Lasso", "ENET", "RF")
colnames(sim_2b_rsq) <- c("Ridge", "Lasso", "ENET", "RF")
colnames(sim_3a_rsq) <- c("Ridge", "Lasso", "ENET")
colnames(sim_3b_rsq) <- c("Ridge", "Lasso", "ENET", "RF")

sim_1_mse <- t(sim_1_mse)
sim_1_rsq <- t(sim_1_rsq)


sd(as.numeric(l_sl_mse_log1))
#mean(as.numeric(l_ols_mse_log1))
sd(as.numeric(l_ridge_mse_log1))
sd(as.numeric(l_lasso_mse_log1))
sd(as.numeric(l_elasticnet_mse_log1))
sd(as.numeric(l_rf_mse_log1))
sd(as.numeric(l_boost_mse_log1))




par(mfrow=c(1,1))


### Plot Comparison LASSO x ENET
t<-c(0,1)

### Scenario 1
plot(sim_2a_rsq[,2],sim_2a_rsq[,3], xlim = t, ylim=t, pch = 21, col="royalblue", cex=2, bg="lightblue", lwd=1, main = "R-Square",
     xlab = "Lasso", ylab = "ENET")
abline(0,1,col="red", lwd = 2)

### Scenario 2
plot(sim_3a_rsq[,2],sim_3a_rsq[,3], xlim = t, ylim=t, pch = 21, col="royalblue", cex=2, bg="lightblue", lwd=1, main = "R-Square",
     xlab = "Lasso", ylab = "ENET")
abline(0,1,col="red", lwd = 2)

### Scenario 3
plot(rsq3[,1],rsq3[,2], xlim = t, ylim=t, pch = 21, col="royalblue", cex=2, bg="lightblue", lwd=1, main = "R-Square",
     xlab = "Lasso", ylab = "ENET")
abline(0,1,col="red", lwd = 2)

### Scenario 4
plot(rsq4[,2],rsq4[,3], xlim = t, ylim=t, pch = 21, col="royalblue", cex=2, bg="lightblue", lwd=1, main = "R-Square",
     xlab = "Lasso", ylab = "ENET")
abline(0,1,col="red", lwd = 2)




### Comparison of prediction perfomance
par(mfrow=c(2,2))
t2 <- c(0,0.35)

### Simulation 2a
plot(sim_2a,sim_1, xlim = t2, ylim= t2, main = "Simulation 2a",
    cex.lab = 1.5, cex.axis=1.5, cex.main = 1, cex = 2, pch = 19, ylab = "Simple Linear Regression", xlab = "R^2 values", 
    col = c("royalblue","royalblue","royalblue","royalblue", "royalblue", "royalblue", "royalblue",
    "orange","orange", "orange", "orange", "orange", "orange", "orange",
    "red", "red", "red", "red", "red", "red", "red",
    "gray", "gray", "gray", "gray" ,"gray", "gray", "gray"))
legend("topleft", c("Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest"),  fill=c("royalblue","orange","red","gray"))
abline(0,1,col="red")

### Simulation 2b
plot(sim_2b_rsq,sim_1, xlim = t2, ylim= t2, main = "Simulation 2b",
     cex.lab = 1.5, cex.axis=1.5, cex.main = 1, cex = 2, pch = 19, ylab = "Simple Linear Regression", xlab = "R^2 values", 
     col = c("royalblue","royalblue","royalblue","royalblue", "royalblue", "royalblue", "royalblue",
             "orange","orange", "orange", "orange", "orange", "orange", "orange",
             "red", "red", "red", "red", "red", "red", "red",
             "gray", "gray", "gray", "gray" ,"gray", "gray", "gray"))
legend("topleft", c("Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest"),  fill=c("royalblue","orange","red","gray"))
abline(0,1,col="red")

### Simulation 3a
plot(sim_3a,sim_1, xlim = t2, ylim= t2, main = "Simulation 3a",
     cex.lab = 1.5, cex.axis=1.5, cex.main = 1, cex = 2, pch = 19, ylab = "Simple Linear Regression", xlab = "R^2 values", 
     col = c("royalblue","royalblue","royalblue","royalblue", "royalblue", "royalblue", "royalblue",
             "orange","orange", "orange", "orange", "orange", "orange", "orange",
             "red", "red", "red", "red", "red", "red", "red",
             "gray", "gray", "gray", "gray" ,"gray", "gray", "gray"))
legend("topleft", c("Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest"),  fill=c("royalblue","orange","red","gray"))
abline(0,1,col="red")

### Simulation 3b
plot(sim_3b,sim_1, xlim = t2, ylim= t2, main = "Simulation 3b",
     cex.lab = 1.5, cex.axis=1.5, cex.main = 1, cex = 2, pch = 19, ylab = "Simple Linear Regression", xlab = "R^2 values", 
     col = c("royalblue","royalblue","royalblue","royalblue", "royalblue", "royalblue", "royalblue",
             "orange","orange", "orange", "orange", "orange", "orange", "orange",
             "red", "red", "red", "red", "red", "red", "red",
             "gray", "gray", "gray", "gray" ,"gray", "gray", "gray"))
legend("topleft", c("Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest"),  fill=c("royalblue","orange","red","gray"))
abline(0,1,col="red")


par(mfrow=c(2,2))

### cis-SNPs only
plot(snps2a[,4], snps2a[,6], pch = 19, ylim = c(-5, 0), col = "royalblue", xlab="Chromosome", ylab="log10(|beta|)", xaxt='n')
axis(1, at = seq(1, 22, 1))
legend("bottomright", c("cis-SNPs","trans-SNPs"), fill = c("royalblue","sienna"))
abline(v=6, col="lightblue", lty = 2, lwd = 2)
text(6, -4, "LINC00044")

plot(snps2b[,4], snps2b[,6], pch = 19, ylim = c(-5, 0), col = "royalblue", xlab="Chromosome", ylab="log10(|beta|)", xaxt='n')
axis(1, at = seq(1, 22, 1))
legend("bottomright", c("cis-SNPs","trans-SNPs"), fill = c("royalblue","sienna"))
abline(v=6, col="lightblue", lty = 2, lwd = 2)
text(6, -4, "LINC00044")

### trans-SNP
plot(snps3b[-highSNPS,4], snps3b[-highSNPS,6], ylim = c(-20, 0), pch = 19, col = "sienna", xlab="Chromosome", ylab="log10(|beta|)", xaxt='n')
axis(1, at = seq(1, 22, 1))
points(snps3b[highSNPS,4], snps3b[highSNPS,6], pch = 19, col = "royalblue")
legend("bottomright", c("cis-SNPs","trans-SNPs"), fill = c("royalblue","sienna"))
abline(v=6, col="lightblue", lty = 2, lwd = )
text(6, -15, "LINC00044")

plot(snps3a[-highSNPS,4], snps3a[-highSNPS,6], ylim = c(-5, 0), pch = 19, col = "sienna", xlab="Chromosome", ylab="log10(|beta|)", xaxt='n')
axis(1, at = seq(1, 22, 1))
points(snps3a[highSNPS,4], snps3a[highSNPS,6], pch = 19, col = "royalblue")
legend("bottomright", c("cis-SNPs","trans-SNPs"), fill = c("royalblue","sienna"))
abline(v=6, col="lightblue", lty = 2, lwd = 2)
text(6, -4, "LINC00044")

highSNPS <- c(22:35)
highSNPS <- c(2:9)
ylim = c(0, 10)
ylim = c(-1, 10)
main = "Simulation I"



plot(A_rsq,B_rsq,xlim = c(0, 0.15), ylim = c(0, 0.15), main="R-Squared", 
     xlab="Method 1", ylab="Method 2", cex.lab = 1.5, cex.axis=1.5, cex.main = 1.5, cex = 3, pch = 19, col = c("red","sienna","darkgreen","royalblue","orange","gray"))
legend("topleft", c("Linear Regression","Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest", "GBM"),  fill=c("red","sienna","darkgreen","royalblue","orange","gray"))
abline(0,1, col="red")


plot(Simple_linear_regression,Lasso_regression, xlim = t2, ylim=t2 , main = "MSE")
abline(lsfit(as.numeric(Simple_linear_regression), as.numeric(Lasso_regression)),col="red")

plot(Simple_linear_regression,Elastic_net_regression, xlim = t2, ylim=t2 , main = "MSE")
abline(lsfit(as.numeric(Simple_linear_regression), as.numeric(Elastic_net_regression)),col="red")

plot(Simple_linear_regression,Random_forest_regression, xlim = t2, ylim=t2 , main = "MSE")
abline(lsfit(as.numeric(Simple_linear_regression), as.numeric(Random_forest_regression)),col="red")



data<-data.frame(Simple_linear_regression=as.list(Simple_linear_regression_rsq),
                 Ridge_regression=as.list(Ridge_regression_rsq),
                 Lasso_regression=as.list(Lasso_regression_rsq),
                 Elastic_net_regression=as.list(Elastic_net_regression_rsq),
                 Random_forest_regression=as.list(Random_forest_regression_rsq))



boxplot(chr21_mse, ylim = c(0, 0.1))
boxplot(chr21b_mse, ylim = c(0, 0.1))

par(mfrow=c(2,2)) #combine multiple plots into one overall graph 

boxplot(chr21_mse, ylim = c(0, 0.1), main ="Method 1", col = 
          c("red","sienna","darkgreen","royalblue","orange","gray"), outline = T,ylab="MSE")
boxplot(chr21b_mse, ylim = c(0, 0.1), main ="Method 2", col = 
          c("red","sienna","darkgreen","royalblue","orange","gray"), outline = T,ylab="MSE")

### R squared
boxplot(sim_1_rsq, ylim = c(0, 0.05), main ="Simulation 1", col = 
          c("royalblue"), outline = T, xlab = "Simple Linear Regression", ylab="R^2")
boxplot(sim_2a_rsq, ylim = c(0.1, 0.35), main ="Simulation 2a", col = 
          c("royalblue","orange","red","gray"), outline = T,ylab="R^2")
boxplot(rsq, ylim= c (0, 0.8), main ="Simulation 2b", col = 
          c("royalblue","orange","red","gray"), outline = T,ylab="R^2")
boxplot(sim_3a_rsq, ylim = c(0.1, 0.35), main ="Simulation 3a", col = 
          c("royalblue","orange","red","gray"), outline = T,ylab="R^2")
boxplot(sim_3b_rsq, ylim = c(0, 0.35), main ="Simulation 3b", col = 
          c("royalblue","orange","red","gray"), outline = T,ylab="R^2")



plot(A_mse,B_mse,xlim = c(0.005, 0.015), ylim = c(0.005, 0.015), main="MSE", 
     xlab="Method 1", ylab="Method 2", cex.lab = 1.5, cex.axis=1.5, cex.main = 1.5, cex = 3, pch = 19, col = c("red","sienna","darkgreen","royalblue","orange","gray"))
legend("topleft", c("Linear Regression","Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest", "GBM"),  fill=c("red","sienna","darkgreen","royalblue","orange","gray"))
abline(0,1, col="red")

plot(A_rsq,B_rsq,xlim = c(0, 0.15), ylim = c(0, 0.15), main="R-Squared", 
     xlab="Method 1", ylab="Method 2", cex.lab = 1.5, cex.axis=1.5, cex.main = 1.5, cex = 3, pch = 19, col = c("red","sienna","darkgreen","royalblue","orange","gray"))
legend("topleft", c("Linear Regression","Ridge Regression", "Lasso Regression", "Elastic Net", "Random Forest", "GBM"),  fill=c("red","sienna","darkgreen","royalblue","orange","gray"))
abline(0,1, col="red")

boxplot(sim_2b_rsq, ylim = c(0, 0.9),col = c("royalblue","orange","red","gray"), outline = T, ylab="R^2")

plot(dif, col="green")
abline(v = y)

barplot(dif2[,2],  xlab="Predicted genes", ylab="Difference between methods", ylim = c(-0.4, 0.2), col="royalblue")
abline(h= -0.1, col="royalblue", lty = 2, lwd = 2 )

### Comparing two methods
A_mse <- c()
A_rsq <- c()
B_mse <- c()
B_rsq <- c()


for(i in 1:6){
  A_mse[i] <- mean(as.numeric(chr21_mse[,i]), na.rm = TRUE)
  A_rsq[i] <- mean(as.numeric(chr21_rsq[,i]), na.rm = TRUE)
  B_mse[i] <- mean(as.numeric(chr21b_mse[,i]), na.rm = TRUE)
  B_rsq[i] <- mean(as.numeric(chr21b_rsq[,i]), na.rm = TRUE)
}
