install.packages("vcfR")
library(vcfR)

#LOAD DATASET
#snp <- read.vcfR("../Dataset/GEUVADIS.chr21.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf")
snp <- read.vcfR("../Dataset/GEUVADIS.chr6.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz")

################## START PREPROCESSING ##################
start <- Sys.time()

tmp <- grep("indel", snp@fix[,3]) # find indels in raw data
length(tmp) # show how many indels in raw data

snpmaf <- INFO2df(snp) # split info into dataframe
snpmaf <- subset(snpmaf, snpmaf$EUR_AF >= 0.05 & snpmaf$AFR_AF >= 0.05) # apply maf filtering
mafIndex <- as.numeric(rownames(snpmaf)) # get the index of filtered snps

end<-Sys.time()
end-start 

start <- Sys.time()

snp_loc <- snp@fix[mafIndex,1:3]  # set location vector
snp_inf <- extract.gt(snp)        # set information vector
snp_inf <- snp_inf[mafIndex,]

tmp <- grep("indel", snp_loc[,3]) # find indels
snp_loc <- snp_loc[-tmp,]         # remove indels from locations
snp_inf <- snp_inf[-tmp,-(1:42)]  # remove indels from genotypes info (-41 first columns) and -1 format
snp_inf <- snp_inf[,-180]         # unmatched column with gene expression
snp_inf <- snp_inf[,-75]          # unmatched column with gene expression

snp_inf[snp_inf=="0|0"] <- 0
snp_inf[snp_inf=="0|1"] <- 1
snp_inf[snp_inf=="1|0"] <- 1
snp_inf[snp_inf =="1|1"] <- 2

snp_inf<-snp_inf[apply(snp_inf,1,function(snp_inf) !all(snp_inf==0)),]
snp_inf<-snp_inf[apply(snp_inf,1,function(snp_inf) !all(snp_inf==1)),] 
snp_inf<-snp_inf[apply(snp_inf,1,function(snp_inf) !all(snp_inf==2)),]

dim(snp_loc) #location dim
dim(snp_inf) #information dim

gMargin <- round(dim(snp_inf)[2]/5) # filtering genotypes <= 20% of the data
snp_inf<-snp_inf[(rowSums(snp_inf=="0") > gMargin),]
snp_inf<-snp_inf[(rowSums(snp_inf=="1") > gMargin),] 
snp_inf<-snp_inf[(rowSums(snp_inf=="2") > gMargin),]

tmp <- match(rownames(snp_inf), snp_loc[,3]) # update snp_loc
snp_loc <- snp_loc[tmp,]

dim(snp_loc) #location dim
dim(snp_inf) #information dim

end<-Sys.time()
end-start 
################## FINISH PREPROCESSING ##################

#write.table(snp_loc, "../Dataset/snp_loc21.txt")
#write.table(snp_inf, "../Dataset/snp_inf21.txt")
write.table(snp_loc, "../Dataset/snp_loc6.txt")
write.table(snp_inf, "../Dataset/snp_inf6.txt")