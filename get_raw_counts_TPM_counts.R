library(sleuth)
library(biomaRt)
library("GO.db")
library(org.Hs.eg.db)
library(labdsv)
library(tximport)
library("readr")
library(fastICA)
library(sva)

library(openxlsx)
library(tidyverse)
library(curl)
library(ggplot2)

library(viper)
library(progeny)
library(fgsea)
library(OmnipathR)

library(patchwork)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(limma)
library(WriteXLS)
### rm(list=ls())
#########################
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

t2g <- t2g[!duplicated(t2g$target_id),]
t2g <- t2g[!t2g$ext_gene == "",]
rownames(t2g) <- t2g$target_id
#################################
base_dir <- "/Users/buschlab/ncrna/batch2_b20_2"
samples <- dir(file.path(base_dir))
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- read.table("/Users/buschlab/snsf/snsf0405.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate( s2c, path = kal_dirs)

rownames(s2c) <- s2c$Sample
print(s2c)
#####################################

file <- s2c$path
file <- paste0(file,"/abundance.h5")


txi <- tximport(file, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T, countsFromAbundance = "lengthScaledTPM")
cnts <- txi$abundance
#sel.rn <- rowSums(cnts) > 10
#cnts2 <- cnts[sel.rn,]
cnts2 <- cnts
colnames(cnts2) <- s2c$Sample
cnts2 <- as.data.frame(cnts2)
# cnts2["DNM1",]
# myCPM <- cpm(cnts2)
WriteXLS(c("cnts2"),"TPM_counts_totalRNA2.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

#####################################
#####################################

file <- s2c$path
file <- paste0(file,"/abundance.h5")


txi <- tximport(file, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
cnts <- txi$counts

cnts2 <- cnts
colnames(cnts2) <- s2c$Sample
cnts2 <- as.data.frame(cnts2)

WriteXLS(c("cnts2"),"raw_counts_totalRNA2.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

#####################################
#####################################


cnts2$ext_gene <- rownames(cnts2)
rownames(cnts2) <- NULL
cnts3 <- cnts2
cnts3 <- cnts3[!duplicated(cnts3$ext_gene),]

anno <- read.table("hg38.bf.anno")
colnames(anno) <- c("ensembl", "ext_gene")
anno2 <- anno[!duplicated(anno$ext_gene),]


whole <- merge(cnts2, anno2, by = "ext_gene", all.x = TRUE)

WriteXLS(c("whole"),"raw_counts_totalRNA_add_ensembl.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
















