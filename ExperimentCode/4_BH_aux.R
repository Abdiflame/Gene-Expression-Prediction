gene <- read.table("../Dataset/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",header = TRUE)
cis_snps2 <- read.table("../Dataset/cis_snps_names_chr21.txt", header = TRUE)
cis_snps_pvalue2 <- read.table("../Dataset/pvalue_21.txt", header = TRUE)

#######################################
### Start Pre-processing
gene_n <- gene[grep("^21$", gene[,3]),-(1:4)] # Select all genes in a determined Chromosome
gene_names <- gene[grep("^21$",gene[,3]),2] # Extract gene names
row.names(gene_n) <- gene_names
gene_n[1:5,1:5] # Genes X Samples (Gene Expression values)

colnames(cis_snps2) <- gene_names

start<-Sys.time()

snp_names_aux <- c()
snp_adjpvalue_aux <- c()
BH_passed_snps_final2 <- rep(NA, dim(cis_snps2)[1])

for(i in 1:dim(cis_snps_pvalue2)[2]) # Benjamini-Hochberg (multiple comparison all snps per gene)
{
  auxVec <- p.adjust(as.numeric(cis_snps_pvalue2[!is.na(cis_snps_pvalue2[,i]),i]), method = "BH")
  snp_adjpvalue_aux <- c(auxVec)                                      # Getting the adjusted p-value
  snp_names_aux <- c(as.character(cis_snps2[!is.na(cis_snps2[,i]),i]))  # Getting the SNPs names
  df_aux <- data.frame(snp_names_aux, snp_adjpvalue_aux)
  names(df_aux) <- c("SNP", "AdjPvalue")
  #dim(df_aux)
  BH_passed_snps_aux <- subset(df_aux, df_aux$AdjPvalue < 0.5)       # Applying the p-value threshold
  #dim(BH_passed_snps_aux)
  BH_passed_snps_final2 <- cbind.fill(BH_passed_snps_final2, as.character(BH_passed_snps_aux$SNP), fill = NA) # Making the final SNP's matrix
  
  print("Iteration")
  print(i)
}

end<-Sys.time()
end-start

BH_passed_snps_final2 <- BH_passed_snps_final2[,-1]
colnames(BH_passed_snps_final2) <- colnames(cis_snps2) 
BH_passed_snps_final2 <- BH_passed_snps_final2[,colSums(is.na(BH_passed_snps_final2)) < (nrow(BH_passed_snps_final2)-1)]

dim(BH_passed_snps_final2)
dim(BH_passed_snps_final)

BH_passed_snps_final_master <- BH_passed_snps_final[,na.omit(match(colnames(BH_passed_snps_final2),colnames(BH_passed_snps_final)))] 

dim(BH_passed_snps_final_master)
