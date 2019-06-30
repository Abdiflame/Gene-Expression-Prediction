library(FSA)
library(data.table)
library(rowr)

#LOAD DATASET
cis_snps <- read.table("/home/lffraga/Table Data/cis_snps_names_chr21b.txt", header = TRUE)
snps <- read.table("/home/lffraga/Table Data/snp_inf21.txt", header = TRUE)
test <- read.table("/home/lffraga/Table Data/cis_matrix_chr21b.txt", header = TRUE)
gene <- read.table("/home/lffraga/Raw Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",header = TRUE)

### Start Preprocessing
gene_n <- gene[grep("^21$", gene[,3]),-(1:4)] # Select all genes in a determined Chromosome
gene_names <- gene[grep("^21$",gene[,3]),2] # Extract gene names
row.names(gene_n) <- gene_names
gene_n[1:5,1:5] # Genes X Samples (Gene Expression values)

gene_n <- gene_n[,match(colnames(snps),colnames(gene_n))]
dim(gene_n)
n_gene <- t(gene_n)
dim(n_gene)
#n_gene[1:5,1:5] # Samples X Genes (Gene Expression values)

snp <- t(snps)
snp <- snp[,match(rownames(test), colnames(snp))] # Samples x SNPS
#snp[1:5,1:5]
################# FINISH PREPROCESSING #################

### Log Dimension
log_gene<-log(n_gene+(1-min(n_gene)))
log_gene[1:5,1:5]
### End of Log

### Getting P-Value ANOVA
#l <- list()

#for(i in 1:dim(n_gene)[2]) #columns number
#{
#  cis_snp<-match(as.character(unique(cis_snps[,i])),colnames(snp))
#  t<-snp[,na.omit(cis_snp)]
#  l[i]<-dim(t)[2]  
#}

#length <- max(as.numeric(l))
length <- dim(snp_inf)[1]

start<-Sys.time()

c<-list()

l_p <- rep(NA, length)
cis_snp<-match(as.character(unique(snpslist)),colnames(snp))
t <- snp[,na.omit(cis_snp)]
gene_t <- cbind(t, log_gene[,i])
end <- dim(gene_t)[2]-1

for(j in 1:end)
{
  l_p[j] <- anova(lm(gene_t[,end+1] ~ gene_t[,j]))[5][1,]
}

c <- cbind(c,l_p) #column bind, add column


cis_snps_pvalue <- c
colnames(cis_snps) <- gene_names
colnames(cis_snps_pvalue) <- gene_names
dim(cis_snps_pvalue)
dim(cis_snps)

end<-Sys.time()
end-start
##############################################################################################################################
### Benjamini-Hochberg Method
start<-Sys.time()

snp_names_aux <- c()
snp_adjpvalue_aux <- c()
BH_passed_snps_final <- rep(NA, length)

auxVec <- p.adjust(cis_snps_pvalue, method = "BH")
snp_adjpvalue_aux <- c(auxVec)                                      # Getting the adjusted p-value
snp_names_aux <- c(as.character(rownames(snp_inf)))  # Getting the SNPs names
df_aux <- data.frame(snp_names_aux, snp_adjpvalue_aux)
names(df_aux) <- c("SNP", "AdjPvalue")
dim(df_aux)
BH_passed_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.5)       # Applying the p-value threshold
#dim(BH_passed_snps_aux)
#BH_passed_snps_final <- cbind.fill(BH_passed_snps_final, as.character(BH_passed_snps_aux$SNP), fill = NA) # Making the final SNP's matrix


end<-Sys.time()
end-start
### Benjamini-Hochberg Method End

### Matrix Cleaning
BH_passed_snps_final <- BH_passed_snps_final[,-1]
colnames(BH_passed_snps_final) <- colnames(cis_snps) 
BH_passed_snps_final <- BH_passed_snps_final[,colSums(is.na(BH_passed_snps_final)) < (nrow(BH_passed_snps_final)-1)]
dim(BH_passed_snps_final)

### Counting the highest number of SNP in a gene
countSNPS <- 0
gene_size <- c()
tmp = 0

for(i in 1:dim(BH_passed_snps_final_master)[2]) ##BH_passed_snps_final
{
  gene_size[i] = sum(!is.na(BH_passed_snps_final_master[,i])) ##BH_passed_snps_final
  print(gene_size[i])
  countSNPS <- countSNPS + gene_size[i]
}
countSNPS # Check total

max_gene_size <- max(as.numeric(gene_size)) 
max_gene_size # Max gene size

BH_passed_snps_final_master <- BH_passed_snps_final_master[1:max_gene_size-1,] # Remove exceeding rows ## BH_passed_snps_final
dim(BH_passed_snps_final_master)

############################################## TA BOM AQUI TB ##########################

### Saving data
BH_gene <- log_gene[,match(colnames(BH_passed_snps_final_master),colnames(log_gene))]  ### BH_passed_snps_final
write.table(cis_snps_pvalue,"/home/lffraga/R_studio/pValue Data/pvalue_21b.txt")             ### saving pValue 
write.table(BH_gene,"/home/lffraga/R_studio/Table Data/BH_gene21b.txt")                      ### All passed genes after (BH method)
write.table(BH_passed_snps_aux, "/home/lffraga/R_studio/Table Data/BH_passed_snps21c.txt") ### All passed SNPs after (BH method)
