# Gene Expression Prediction Using a sparse Genome-Wide SNP Dataset
Work on large-scale Machine Learning for gene expression prediction based on genotype only.

# Introduction
In this work, we used millions of SNPs identified by whole-genome sequencing (WGS) to predict gene expression. Due to the sparsity of existing SNP datasets, we applied four regularized regression methods—Ridge, Lasso, Elastic Net (ENET), and Random Forest—combined with or without predictor (SNP) filtering by proximity to and correlation with the expression of the target gene.

# Dataset
Geuvadis consortium data set of 462 unrelated human lymphoblastoid cell line samples from 5 populations from the 1000 Genomes project.
*	CEPH (CEU), Finns (FIN), British (GBR), Toscani (TSI) and Yoruba (YRI)
*	462 individuals with mRNA and 452 individuals with miRNA data
*	421 in the 1000 Genomes Phase 1 dataset + 41 in Phase 2

Dataset format RPKM and VCF file:
*	Gene Expression: RPKM – Reads Per Kilobase Million
*	Genotypes: VCF – Variant Call Format

Downloaded Files:
GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz  (86.6 MB)
https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/

Chromosome 1-22 (44.61 GB)
https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/genotypes/

The dataset is available on the Geuvadis data browser webpage.

Link: https://www.ebi.ac.uk/Tools/geuvadis-das/

# Work Environment
Rstudio: Version 1.1.463

R Language: R version 3.5.2 (2018-12-20)

Link: https://www.rstudio.com/products/rstudio/download/

Used packages:
* Read VCF files: {vcfR}
* Benjamini-Hochberg: p.adjust {stats v. 3.5.2}
* Ridge, Lasso, Enet: glmnet {glmnet v. 2.0-16}
* Random Forest: randomForest {randomForest v. 4.6-14}
