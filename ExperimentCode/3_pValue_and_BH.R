library(FSA)
library(data.table)
library(rowr)

############### Load Dataset ###############

cis_snps <- read.table("../Dataset/cis_snps_names_chr6.txt", header = TRUE)
snp_inf <- read.table("../Dataset/snp_inf6.txt", header = TRUE)
test <- read.table("../Dataset/cis_matrix_chr6.txt", header = TRUE)
gene <- read.table("../Dataset/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",header = TRUE)

############### Start Pre-processing ############### 

gene_n <- gene[grep("^6$", gene[,3]),-(1:4)]          ## Select all genes in a determined Chromosome
gene_names <- gene[grep("^6$",gene[,3]),2]            ## Extract gene names
row.names(gene_n) <- gene_names
gene_n[1:5,1:5] # Genes X Samples (Gene Expression values)

gene_n <- gene_n[,match(colnames(snp_inf),colnames(gene_n))]
dim(gene_n)
n_gene <- t(gene_n)
dim(n_gene)

snp <- t(snp_inf)
snp <- snp[,match(rownames(test), colnames(snp))] # Samples x SNPS

############### Finish Pre-processing ############### 

### Log Dimension
log_gene <- log(n_gene+(1-min(n_gene)))
log_gene[1:5,1:5]
### End of Log

############### This part is used only for cis-SNPs ###############
## Please comment this block if you are using trans-SNPs

l <- list()

for(i in 1:dim(n_gene)[2]) #columns number
{
  cis_snp<-match(as.character(unique(cis_snps[,i])),colnames(snp))
  t<-snp[,na.omit(cis_snp)]
  l[i]<-dim(t)[2]  
}

length <- max(as.numeric(l))
############### This part is used only for cis-SNPs ############### 

#length <- dim(snp_inf)[1]                                    ## for trans-SNPs

############### Prepare data dimensions ############### 

c<-list()

cis_snp <- match(as.character(rownames(snp_inf)),colnames(snp))
t <- snp[,na.omit(cis_snp)]
gene_t <- cbind(t, log_gene[,dim(n_gene)[2]])
colnames(cis_snps) <- gene_names
#end <- dim(gene_t)[2]-1                                      ## for trans-SNPs
#end <- dim(cis_snps)[1]                                      ## for cis-SNPs

############### Benjamini-Hochberg Method - Start ###############
count <- 0
print("START")
BH_passed_snps_final <- rep(NA, length)
snp_names_aux <- c()
snp_adjpvalue_aux <- c()
l_p <- rep(NA, length)

for(i in 1:dim(gene_n)[1])
{
  out <- paste0("Current iteration: ", i)
  print(out)                                                

  BH_passed_snp_tmp <- snp[,na.omit(match(as.character(cis_snps[,i]),colnames(snp)))]
  ex_data <- cbind(BH_passed_snp_tmp, n_gene[,i])
  end <- dim(ex_data)[2]-1
  
  for(p in 1:end)
  {
    l_p[p] <- anova(lm(ex_data[,end+1] ~ ex_data[,p]))[5][1,]
  }
  
  cis_snps_pvalue <- as.numeric(l_p)
  auxVec <- p.adjust(cis_snps_pvalue, method = "BH")                      # Getting the SNPs adjusted pvalue
  #df_aux <- data.frame(auxVec)
  #names(df_aux) <- "AdjPvalue"
  
  #snp_names_aux <- c(as.character(rownames(snp_inf)))                    # Getting the trans-SNPs names 
  snp_names_aux <- c(as.character(cis_snps[,i]))                          # Getting the cis-SNPs names
  
  snp_adjpvalue_aux <- c(auxVec)
  df_aux <- data.frame(snp_names_aux, snp_adjpvalue_aux)
  names(df_aux) <- c("SNP", "AdjPvalue")
  
  BH_passed_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.5)            # Applying the p-value threshold
  BH_passed_snps_final <- cbind.fill(BH_passed_snps_final, as.character(BH_passed_snps_aux$SNP), fill = NA) # Making the final SNP's matrix
  
  out <- paste0("Passed SNPs: ", dim(BH_passed_snps_aux)[1])
  print(out)  
  
  if(dim(BH_passed_snps_aux)[1] > 0){
    count <- count + 1
    out <- paste0("Passed Gene count: ", count)
    print(out) 
  }
}

################################################################

#snp_names_aux <- c()
#snp_adjpvalue_aux <- c()
#BH_passed_snps_final <- rep(NA, length)

#auxVec <- p.adjust(cis_snps_pvalue, method = "BH")
#snp_adjpvalue_aux <- c(auxVec)                                          # Getting the adjusted p-value
#snp_names_aux <- c(as.character(rownames(snp_inf)))                     # Getting the trans-SNPs names 
#snp_names_aux <- c(as.character(rownames(cis_snps)))                    # Getting the cis-SNPs names
#df_aux <- data.frame(snp_names_aux, snp_adjpvalue_aux)
#names(df_aux) <- c("SNP", "AdjPvalue")
#dim(df_aux)
#BH_passed_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.5)            # Applying the p-value threshold
#dim(BH_passed_snps_aux)


############### Benjamini-Hochberg Method - End ###############

### Matrix Cleaning
BH_passed_snps_final <- BH_passed_snps_final[,-1]
colnames(BH_passed_snps_final) <- colnames(cis_snps) 
BH_passed_snps_final <- BH_passed_snps_final[,colSums(is.na(BH_passed_snps_final)) < (nrow(BH_passed_snps_final)-1)]
dim(BH_passed_snps_final)


### Counting the highest number of SNP in a gene
countSNPS <- 0
gene_size <- c()
tmp = 0

for(i in 1:dim(BH_passed_snps_final)[2]) ##BH_passed_snps_final
{
  gene_size[i] = sum(!is.na(BH_passed_snps_final[,i])) ##BH_passed_snps_final
  #print(gene_size[i])
  countSNPS <- countSNPS + gene_size[i]
}
print(countSNPS) # Check total SNPs

max_gene_size <- max(as.numeric(gene_size)) 
max_gene_size # Max gene size

BH_passed_snps_final_master <- BH_passed_snps_final[1:max_gene_size-1,] # Remove exceeding rows ## BH_passed_snps_final
dim(BH_passed_snps_final_master)

############### Saving tables/variables data files ###############

### Saving data
BH_gene <- log_gene[,match(colnames(BH_passed_snps_final_master),colnames(log_gene))]  ### BH_passed_snps_final
write.table(cis_snps_pvalue,'../Dataset/pvalue_6.txt')                                 ### Saving pValue 
write.table(BH_gene, '../Dataset/BH_gene6.txt')                                        ### All passed genes after (BH method)
write.table(BH_passed_snps_final_master, '../Dataset/BH_passed_snps6.txt')             ### All passed SNPs after (BH method)

