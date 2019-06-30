library(rowr)
install.packages("rowr")

### Load dataset
snp_loc <- read.table("/home/lffraga/Table Data/snp_loc21.txt",header = TRUE)
snp_inf <- read.table("/home/lffraga/Table Data/snp_inf21.txt",header = TRUE)
gene <- read.table("/home/lffraga/Raw Data/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",header = TRUE)

### Gene expression match
start<-Sys.time()

#select chromosome in gene expression dataset
gene_loc <- gene[grep("^21$", gene[,3]),]
geneLocation <- gene_loc[,4]
snps <- t(snp_inf)

#t <- gene_expression_chr[,5:466] # 462 samples
#match_sample_g <- match(rownames(snps),names(gene_expression_chr[,5:466]))
#match_sample_g <- na.omit(match_sample_g)
#length(match_sample_g)

#expr <- data.frame(gene_expression_chr[,1], t[,match_sample_g])
#dim(expr)
#expr[1:5,1:5]

#gene_ex<-expr[-1] #delete first column
#gene_ex[1:5,1:5]
#dim(gene_ex)


#gene_loc <- data.frame(gene_expression_chr[,1],gene_expression_chr[,3:4])
#dim(gene_loc)
#row.names(gene_loc)<-NULL

end<-Sys.time()
end-start

### Cis-SNPs marking
start<-Sys.time()

#geneLocation <- gene_loc[,3]
cis_matrix <- list() 
tmp <- list()

for(i in 1:length(geneLocation))
{
  #tmp <- match(colnames(snp_inf),snp_loc[, 3]) # ALL SNPs
  #tmp <- match(colnames(snps),snp_loc[(geneLocation[i]-10^7 < snp_loc[,2] & snp_loc[,2] < geneLocation[i]+10^7), 3]) # 10Mb
  #tmp <- match(colnames(maf5_snps),maf5_snp_loc[(geneLocation[i]-(10^6)*5 < maf5_snp_loc[,2] & maf5_snp_loc[,2] < geneLocation[i]+(10^6)*5), 3]) # 5Mb
  tmp <- match(colnames(snps),snp_loc[(geneLocation[i]-10^6 < snp_loc[,2] & snp_loc[,2] < geneLocation[i]+10^6), 3]) # 1Mb
  cis_matrix <- cbind(cis_matrix, tmp)
  print(i)
}

dim(cis_matrix)
colnames(cis_matrix) <- gene_loc[,1]
rownames(cis_matrix) <- snp_loc[,3]

end<-Sys.time()
end-start

################################################ FINAL PART START ################################################
#names(cis_matrix[!is.na(cis_matrix[,1]),1])
l <- list()
c <- list()

for(i in 1:dim(cis_matrix)[2])
{
  l[i]<-sum(!is.na(cis_matrix[,i])) #count the number of variants in a gene
}

max_snp_cnt <- max(as.numeric(l)) #get the highest number of variants in a gene

t_size <- rep(NA,max_snp_cnt) #replicates values in NA?
c <- cbind(c,t_size)
dim(c)

for(i in 1:dim(cis_matrix)[2])
{
  c <- cbind.fill(c, names(cis_matrix[!is.na(cis_matrix[,i]),i]), fill = NA) ## errado aqui
}

#c <- c[,-1]
#colnames(c) <- colnames(BH_passed_snp)
#dim(c)
#dim(c[,-1])
#dim(cis_matrix)
#c[1:5,1:5]
#cis_matrix[1:5,1:5]
################################################ FINAL PART END ################################################

### Saving data
write.table(c[,-1],"/home/lffraga/Table Data/cis_snps_names_chr21b.txt")
write.table(cis_matrix, "/home/lffraga/Table Data/cis_matrix_chr21b.txt")
