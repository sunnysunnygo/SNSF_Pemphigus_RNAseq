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

sub1 <- subset(s2c, (Stimulus=="hIgG"|Stimulus=="PX43") & Time=="24h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
#####################3
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T, countsFromAbundance = "lengthScaledTPM")
cnts <- txi$abundance
#sel.rn <- rowSums(cnts) > 10
#cnts2 <- cnts[sel.rn,]
cnts2 <- cnts
colnames(cnts2) <- sub1$Sample
cnts2 <- as.data.frame(cnts2)
cnts2["DNM1",]
myCPM <- cpm(cnts2)
WriteXLS(c("cnts2"),"TPM_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_PX43vshIgG_24h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_PX43vshIgG_24h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="hIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_PX43_vs_hIgG", type="ashr")
# res_px43_24h <- as.data.frame(res[order(res$pvalue),])
res3 <- as.data.frame(res2)
res_px43_24h <- as.data.frame(res2[order(res2$pvalue),])
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
##WriteXLS(c("res_px43_24h"),"DESeq_PX43vshIgG_24h_deg_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

WriteXLS(c("res3"),"DESeq_PX43vshIgG_24h_deg_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#################
pdf("deseq_px43vshIgG_24h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
px43_24_rlog <- rlog(ddsTxi)
px43_24_rlog <- assay(px43_24_rlog)
px43_24_rlog <- as.data.frame(px43_24_rlog)
heat_input2 <- px43_24_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]
pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                #ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'PX43 vs hIgG 24h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


###################
################# 10h PX43 vs hIgG
#####################################
sub1 <- subset(s2c, (Stimulus=="hIgG"|Stimulus=="PX43") & Time=="10h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_PX43vshIgG_10h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_PX43vshIgG_10h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="hIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_PX43_vs_hIgG", type="ashr")
# res_px43_10h <- as.data.frame(res[order(res$pvalue),])
res_px43_10h <- as.data.frame(res2[order(res2$pvalue),])
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
#WriteXLS(c("res_px43_10h"),"DESeq_PX43vshIgG_10h_deg_ncRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#################
pdf("deseq_px43vshIgG_10h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
px43_24_rlog <- rlog(ddsTxi)
px43_24_rlog <- assay(px43_24_rlog)
px43_24_rlog <- as.data.frame(px43_24_rlog)
heat_input2 <- px43_24_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]

pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
######################################3
library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'PX43 vs hIgG 10h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()
############## 5h PX43 vs hIgG 
#######################


sub1 <- subset(s2c, (Stimulus=="hIgG"|Stimulus=="PX43") & Time=="5h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
### txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T,txOut = T)
txi <- tximport(file1, type="kallisto",txOut = T) 
####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_PX43vshIgG_5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_PX43vshIgG_5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="hIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_PX43_vs_hIgG", type="ashr")
# res_px43_5h <- as.data.frame(res[order(res$pvalue),])
res_px43_5h <- as.data.frame(res2[order(res2$pvalue),])
#res3 <- res2[ !is.na(res2$pvalue),]  #  not use na.omit
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
#WriteXLS(c("res_px43_5h"),"DESeq_PX43vshIgG_5h_deg_ncRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

res_px43_5h[c('ENST00000435411.6', "ENST00000669843.1","ENST00000659181.1","ENST00000642559.1","ENST00000656623.1","ENST00000646009.1","ENST00000671084.1","ENST00000645995.1","ENST00000654567.1","ENST00000654342.1","ENST00000424170.5","ENST00000428474.5","ENST00000414750.1","ENST00000424655.1","ENST00000645719.1"),]


#################
pdf("deseq_px43vshIgG_5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
px43_5_rlog <- rlog(ddsTxi)
px43_5_rlog <- assay(px43_5_rlog)
px43_5_rlog <- as.data.frame(px43_5_rlog)
heat_input2 <- px43_5_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]

pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'PX43 vs hIgG 5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


###################
#################### AK23 vs mIgG 5h 
sub1 <- subset(s2c, (Stimulus=="mIgG"|Stimulus=="AK23") & Time=="5h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
txi <- tximport(file1, type="kallisto",txOut = T) 
####
####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_AK23vsmIgG_5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23vsmIgG_5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
#keep <- rowSums(counts(ddsTxi)) >= 3
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="mIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_AK23_vs_mIgG", type="ashr")
# res_AK23_5h <- as.data.frame(res[order(res$pvalue),])
res_AK23_5h <- as.data.frame(res2[order(res2$pvalue),])
#res3 <- res2[ !is.na(res2$pvalue),]  #  not use na.omit
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
#WriteXLS(c("res_AK23_5h"),"DESeq_AK23vsmIgG_5h_deg_ncRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

res_AK23_5h[c('ENST00000435411.6', "ENST00000669843.1","ENST00000659181.1","ENST00000642559.1","ENST00000656623.1","ENST00000646009.1","ENST00000671084.1","ENST00000645995.1","ENST00000654567.1","ENST00000654342.1","ENST00000424170.5","ENST00000428474.5","ENST00000414750.1","ENST00000424655.1","ENST00000645719.1"),]
#################
pdf("deseq_AK23vsmIgG_5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
AK23_5_rlog <- rlog(ddsTxi)
AK23_5_rlog <- assay(AK23_5_rlog)
AK23_5_rlog <- as.data.frame(AK23_5_rlog)
heat_input2 <- AK23_5_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]

pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'AK23 vs mIgG 5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


###################

#################### AK23 vs mIgG 10h 
sub1 <- subset(s2c, (Stimulus=="mIgG"|Stimulus=="AK23") & Time=="10h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_AK23vsmIgG_10h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23vsmIgG_10h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="mIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_AK23_vs_mIgG", type="ashr")
# res_AK23_10h <- as.data.frame(res[order(res$pvalue),])
res_AK23_10h <- as.data.frame(res2[order(res2$pvalue),])
#res3 <- res2[ !is.na(res2$pvalue),]  #  not use na.omit
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
#WriteXLS(c("res_AK23_10h"),"DESeq_AK23vsmIgG_10h_deg_ncRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#################
pdf("deseq_AK23vsmIgG_10h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
AK23_5_rlog <- rlog(ddsTxi)
AK23_5_rlog <- assay(AK23_5_rlog)
AK23_5_rlog <- as.data.frame(AK23_5_rlog)
heat_input2 <- AK23_5_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]
pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'AK23 vs mIgG 10h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


###################
#################### AK23 vs mIgG 24h 
sub1 <- subset(s2c, (Stimulus=="mIgG"|Stimulus=="AK23") & Time=="24h")
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_AK23vsmIgG_24h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23vsmIgG_24h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Stimulus)

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Stimulus <- relevel(ddsTxi$Stimulus, ref="mIgG")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
res2 <- results(ddsTxi) ##### 
# res <- lfcShrink(ddsTxi, coef="Stimulus_AK23_vs_mIgG", type="ashr")
# res_AK23_24h <- as.data.frame(res[order(res$pvalue),])
res_AK23_24h <- as.data.frame(res2[order(res2$pvalue),])
#res3 <- res2[ !is.na(res2$pvalue),]  #  not use na.omit
res3 <- na.omit(res2)
top <- res3[res3$padj < 0.05,]
top2 <- top[abs(top$log2FoldChange) > 1,]
#WriteXLS(c("res_AK23_24h"),"DESeq_AK23vsmIgG_24h_deg_ncRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#################
pdf("deseq_AK23vsmIgG_24h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
AK23_5_rlog <- rlog(ddsTxi)
AK23_5_rlog <- assay(AK23_5_rlog)
AK23_5_rlog <- as.data.frame(AK23_5_rlog)
heat_input2 <- AK23_5_rlog[rownames(top),]
heat_input3 <- heat_input2[,c(1,3,5,2,4,6)]

pheatmap(heat_input3, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
library(EnhancedVolcano)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-4, 4),
                ylim = c(0, 10),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 0.35,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'AK23 vs mIgG 24h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


################### enrichment by gage
###########################################
##############
### GSEA
####
# 
library(gage)
# prep database - Hallmark oncogenic pathways
library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens", , category = "H")
h_data_frame = all_gene_sets %>% dplyr::filter(gs_cat == "H") %>% as.data.frame()
HALLMARK.DB = split(x=h_data_frame$gene_symbol,f=h_data_frame$gs_name)
HALLMARK.DB <- lapply(HALLMARK.DB,function(x) unique(x))


all_gene_sets = msigdbr(species = "Homo sapiens", , category = "C2")
r_data_frame = all_gene_sets %>% dplyr::filter(gs_cat == "C2") %>% as.data.frame()
C2.DB = split(x=r_data_frame$gene_symbol,f=r_data_frame$gs_name)
idx <- grep("REACTOME_|PLURINET",names(C2.DB))
C2.DB <- C2.DB[idx]
C2.DB <- lapply(C2.DB,function(x) unique(x))

all_gene_sets = msigdbr(species = "Homo sapiens", category = "C5",subcategory="GO:BP")
r_data_frame = all_gene_sets %>% dplyr::filter(gs_cat == "C5") %>% as.data.frame()
C5.BP= split(x=r_data_frame$gene_symbol,f=r_data_frame$gs_name)
C5.BP <- lapply(C5.BP,function(x) unique(x))
###
add_genes <- function(gsea_result,pathways,diffgenes) {
  tmp <- gsea_result
  tmp$genes <- NULL
  for(i in rownames(tmp)){
    sym <- sort(pathways[[i]])
    syms <- NULL
    effs <- NULL
    sigs <- NULL
    for(ii in 1:length(sym)){
      myeid <- sym[ii]
      idx <- which(rownames(diffgenes) == myeid)
      if(length(idx)){
        syms[ii] <- myeid
        effs[ii] <- diffgenes$log2FoldChange[idx[1]]
        sigs[ii] <- diffgenes$pvalue[idx[1]]
      }
    }
    idx <- order(effs)
    syms <- syms[idx]
    effs <- effs[idx]
    sigs <- sigs[idx]
    idx <- which(is.na(effs))
    if(length(idx) > 0){
      syms <- syms[-idx]
      effs <- effs[-idx]
      sigs <- sigs[-idx]
    }
    sym <- paste(syms,"(",round(effs,3),",",round(sigs,5),")",sep="",collapse=" ")
    tmp[i,"genes"] <- sym
  }
  return(tmp)
}

##########################
library(gage)
###################PX43vshIgG_5h
#######################################
mylist1 <- list() 
res_px43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_5h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_px43vshIgG$log2FoldChange,row.names=rownames(res_px43vshIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Hallmark_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Hallmark_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Hallmark_Up_hIgG"]] <- add_genes(mylist1[["Hallmark_Up_hIgG"]],HALLMARK.DB,res_px43vshIgG)
mylist1[["Hallmark_Down_hIgG"]] <- add_genes(mylist1[["Hallmark_Down_hIgG"]],HALLMARK.DB,res_px43vshIgG)
rownames(mylist1[["Hallmark_Up_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Up_hIgG"]]))
rownames(mylist1[["Hallmark_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Up_hIgG"]]))
rownames(mylist1[["Hallmark_Down_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Down_hIgG"]]))
rownames(mylist1[["Hallmark_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Down_hIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Reactome_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Reactome_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Reactome_Up_hIgG"]] <- add_genes(mylist1[["Reactome_Up_hIgG"]],C2.DB,res_px43vshIgG)
mylist1[["Reactome_Down_hIgG"]] <- add_genes(mylist1[["Reactome_Down_hIgG"]],C2.DB,res_px43vshIgG)
rownames(mylist1[["Reactome_Up_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Up_hIgG"]]))
rownames(mylist1[["Reactome_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Up_hIgG"]]))
rownames(mylist1[["Reactome_Down_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Down_hIgG"]]))
rownames(mylist1[["Reactome_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Down_hIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["GO_BP_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["GO_BP_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["GO_BP_Up_hIgG"]] <- add_genes(mylist1[["GO_BP_Up_hIgG"]],C5.BP,res_px43vshIgG)
mylist1[["GO_BP_Down_hIgG"]] <- add_genes(mylist1[["GO_BP_Down_hIgG"]],C5.BP,res_px43vshIgG)
rownames(mylist1[["GO_BP_Up_hIgG"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Up_hIgG"]]))
rownames(mylist1[["GO_BP_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Up_hIgG"]]))
rownames(mylist1[["GO_BP_Down_hIgG"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Down_hIgG"]]))
rownames(mylist1[["GO_BP_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Down_hIgG"]]))

WriteXLS(c("mylist1"),"GSEA_PX43vshIgG_5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

library(gage)
###################PX43vshIgG_10h
#######################################
mylist2 <- list() 
res_px43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_10h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_px43vshIgG$log2FoldChange,row.names=rownames(res_px43vshIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist2[["Hallmark_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist2[["Hallmark_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist2[["Hallmark_Up_hIgG"]] <- add_genes(mylist2[["Hallmark_Up_hIgG"]],HALLMARK.DB,res_px43vshIgG)
mylist2[["Hallmark_Down_hIgG"]] <- add_genes(mylist2[["Hallmark_Down_hIgG"]],HALLMARK.DB,res_px43vshIgG)
rownames(mylist2[["Hallmark_Up_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist2[["Hallmark_Up_hIgG"]]))
rownames(mylist2[["Hallmark_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["Hallmark_Up_hIgG"]]))
rownames(mylist2[["Hallmark_Down_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist2[["Hallmark_Down_hIgG"]]))
rownames(mylist2[["Hallmark_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["Hallmark_Down_hIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist2[["Reactome_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist2[["Reactome_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist2[["Reactome_Up_hIgG"]] <- add_genes(mylist2[["Reactome_Up_hIgG"]],C2.DB,res_px43vshIgG)
mylist2[["Reactome_Down_hIgG"]] <- add_genes(mylist2[["Reactome_Down_hIgG"]],C2.DB,res_px43vshIgG)
rownames(mylist2[["Reactome_Up_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist2[["Reactome_Up_hIgG"]]))
rownames(mylist2[["Reactome_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["Reactome_Up_hIgG"]]))
rownames(mylist2[["Reactome_Down_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist2[["Reactome_Down_hIgG"]]))
rownames(mylist2[["Reactome_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["Reactome_Down_hIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist2[["GO_BP_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist2[["GO_BP_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist2[["GO_BP_Up_hIgG"]] <- add_genes(mylist2[["GO_BP_Up_hIgG"]],C5.BP,res_px43vshIgG)
mylist2[["GO_BP_Down_hIgG"]] <- add_genes(mylist2[["GO_BP_Down_hIgG"]],C5.BP,res_px43vshIgG)
rownames(mylist2[["GO_BP_Up_hIgG"]]) <- gsub("GOBP_","",rownames(mylist2[["GO_BP_Up_hIgG"]]))
rownames(mylist2[["GO_BP_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["GO_BP_Up_hIgG"]]))
rownames(mylist2[["GO_BP_Down_hIgG"]]) <- gsub("GOBP_","",rownames(mylist2[["GO_BP_Down_hIgG"]]))
rownames(mylist2[["GO_BP_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist2[["GO_BP_Down_hIgG"]]))

WriteXLS(c("mylist2"),"GSEA_PX43vshIgG_10h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################PX43vshIgG_24h
#######################################
mylist3 <- list() 
res_px43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_24h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_px43vshIgG$log2FoldChange,row.names=rownames(res_px43vshIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["Hallmark_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["Hallmark_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["Hallmark_Up_hIgG"]] <- add_genes(mylist3[["Hallmark_Up_hIgG"]],HALLMARK.DB,res_px43vshIgG)
mylist3[["Hallmark_Down_hIgG"]] <- add_genes(mylist3[["Hallmark_Down_hIgG"]],HALLMARK.DB,res_px43vshIgG)
rownames(mylist3[["Hallmark_Up_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist3[["Hallmark_Up_hIgG"]]))
rownames(mylist3[["Hallmark_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["Hallmark_Up_hIgG"]]))
rownames(mylist3[["Hallmark_Down_hIgG"]]) <- gsub("HALLMARK_","",rownames(mylist3[["Hallmark_Down_hIgG"]]))
rownames(mylist3[["Hallmark_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["Hallmark_Down_hIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["Reactome_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["Reactome_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["Reactome_Up_hIgG"]] <- add_genes(mylist3[["Reactome_Up_hIgG"]],C2.DB,res_px43vshIgG)
mylist3[["Reactome_Down_hIgG"]] <- add_genes(mylist3[["Reactome_Down_hIgG"]],C2.DB,res_px43vshIgG)
rownames(mylist3[["Reactome_Up_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist3[["Reactome_Up_hIgG"]]))
rownames(mylist3[["Reactome_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["Reactome_Up_hIgG"]]))
rownames(mylist3[["Reactome_Down_hIgG"]]) <- gsub("REACTOME_","",rownames(mylist3[["Reactome_Down_hIgG"]]))
rownames(mylist3[["Reactome_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["Reactome_Down_hIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["GO_BP_Up_hIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["GO_BP_Down_hIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["GO_BP_Up_hIgG"]] <- add_genes(mylist3[["GO_BP_Up_hIgG"]],C5.BP,res_px43vshIgG)
mylist3[["GO_BP_Down_hIgG"]] <- add_genes(mylist3[["GO_BP_Down_hIgG"]],C5.BP,res_px43vshIgG)
rownames(mylist3[["GO_BP_Up_hIgG"]]) <- gsub("GOBP_","",rownames(mylist3[["GO_BP_Up_hIgG"]]))
rownames(mylist3[["GO_BP_Up_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["GO_BP_Up_hIgG"]]))
rownames(mylist3[["GO_BP_Down_hIgG"]]) <- gsub("GOBP_","",rownames(mylist3[["GO_BP_Down_hIgG"]]))
rownames(mylist3[["GO_BP_Down_hIgG"]]) <- gsub("_"," ",rownames(mylist3[["GO_BP_Down_hIgG"]]))

WriteXLS(c("mylist3"),"GSEA_PX43vshIgG_24h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################AK23vsmIgG_24h
#######################################
mylist4 <- list() 
res_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_24h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_AK23vsmIgG$log2FoldChange,row.names=rownames(res_AK23vsmIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["Hallmark_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["Hallmark_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["Hallmark_Up_mIgG"]] <- add_genes(mylist4[["Hallmark_Up_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
mylist4[["Hallmark_Down_mIgG"]] <- add_genes(mylist4[["Hallmark_Down_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
rownames(mylist4[["Hallmark_Up_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist4[["Hallmark_Up_mIgG"]]))
rownames(mylist4[["Hallmark_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["Hallmark_Up_mIgG"]]))
rownames(mylist4[["Hallmark_Down_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist4[["Hallmark_Down_mIgG"]]))
rownames(mylist4[["Hallmark_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["Hallmark_Down_mIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["Reactome_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["Reactome_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["Reactome_Up_mIgG"]] <- add_genes(mylist4[["Reactome_Up_mIgG"]],C2.DB,res_AK23vsmIgG)
mylist4[["Reactome_Down_mIgG"]] <- add_genes(mylist4[["Reactome_Down_mIgG"]],C2.DB,res_AK23vsmIgG)
rownames(mylist4[["Reactome_Up_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist4[["Reactome_Up_mIgG"]]))
rownames(mylist4[["Reactome_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["Reactome_Up_mIgG"]]))
rownames(mylist4[["Reactome_Down_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist4[["Reactome_Down_mIgG"]]))
rownames(mylist4[["Reactome_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["Reactome_Down_mIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["GO_BP_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["GO_BP_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["GO_BP_Up_mIgG"]] <- add_genes(mylist4[["GO_BP_Up_mIgG"]],C5.BP,res_AK23vsmIgG)
mylist4[["GO_BP_Down_mIgG"]] <- add_genes(mylist4[["GO_BP_Down_mIgG"]],C5.BP,res_AK23vsmIgG)
rownames(mylist4[["GO_BP_Up_mIgG"]]) <- gsub("GOBP_","",rownames(mylist4[["GO_BP_Up_mIgG"]]))
rownames(mylist4[["GO_BP_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["GO_BP_Up_mIgG"]]))
rownames(mylist4[["GO_BP_Down_mIgG"]]) <- gsub("GOBP_","",rownames(mylist4[["GO_BP_Down_mIgG"]]))
rownames(mylist4[["GO_BP_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist4[["GO_BP_Down_mIgG"]]))

WriteXLS(c("mylist4"),"GSEA_AK23vsmIgG_24h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################
###################AK23vsmIgG_10h
#######################################
mylist5 <- list() 
res_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_10h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_AK23vsmIgG$log2FoldChange,row.names=rownames(res_AK23vsmIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["Hallmark_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["Hallmark_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["Hallmark_Up_mIgG"]] <- add_genes(mylist5[["Hallmark_Up_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
mylist5[["Hallmark_Down_mIgG"]] <- add_genes(mylist5[["Hallmark_Down_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
rownames(mylist5[["Hallmark_Up_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist5[["Hallmark_Up_mIgG"]]))
rownames(mylist5[["Hallmark_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["Hallmark_Up_mIgG"]]))
rownames(mylist5[["Hallmark_Down_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist5[["Hallmark_Down_mIgG"]]))
rownames(mylist5[["Hallmark_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["Hallmark_Down_mIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["Reactome_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["Reactome_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["Reactome_Up_mIgG"]] <- add_genes(mylist5[["Reactome_Up_mIgG"]],C2.DB,res_AK23vsmIgG)
mylist5[["Reactome_Down_mIgG"]] <- add_genes(mylist5[["Reactome_Down_mIgG"]],C2.DB,res_AK23vsmIgG)
rownames(mylist5[["Reactome_Up_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist5[["Reactome_Up_mIgG"]]))
rownames(mylist5[["Reactome_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["Reactome_Up_mIgG"]]))
rownames(mylist5[["Reactome_Down_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist5[["Reactome_Down_mIgG"]]))
rownames(mylist5[["Reactome_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["Reactome_Down_mIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["GO_BP_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["GO_BP_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["GO_BP_Up_mIgG"]] <- add_genes(mylist5[["GO_BP_Up_mIgG"]],C5.BP,res_AK23vsmIgG)
mylist5[["GO_BP_Down_mIgG"]] <- add_genes(mylist5[["GO_BP_Down_mIgG"]],C5.BP,res_AK23vsmIgG)
rownames(mylist5[["GO_BP_Up_mIgG"]]) <- gsub("GOBP_","",rownames(mylist5[["GO_BP_Up_mIgG"]]))
rownames(mylist5[["GO_BP_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["GO_BP_Up_mIgG"]]))
rownames(mylist5[["GO_BP_Down_mIgG"]]) <- gsub("GOBP_","",rownames(mylist5[["GO_BP_Down_mIgG"]]))
rownames(mylist5[["GO_BP_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist5[["GO_BP_Down_mIgG"]]))

WriteXLS(c("mylist5"),"GSEA_AK23vsmIgG_10h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################AK23vsmIgG_5h
#######################################
mylist6 <- list() 
res_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_5h_deg_totalRNA.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_AK23vsmIgG$log2FoldChange,row.names=rownames(res_AK23vsmIgG))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["Hallmark_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["Hallmark_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["Hallmark_Up_mIgG"]] <- add_genes(mylist6[["Hallmark_Up_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
mylist6[["Hallmark_Down_mIgG"]] <- add_genes(mylist6[["Hallmark_Down_mIgG"]],HALLMARK.DB,res_AK23vsmIgG)
rownames(mylist6[["Hallmark_Up_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist6[["Hallmark_Up_mIgG"]]))
rownames(mylist6[["Hallmark_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["Hallmark_Up_mIgG"]]))
rownames(mylist6[["Hallmark_Down_mIgG"]]) <- gsub("HALLMARK_","",rownames(mylist6[["Hallmark_Down_mIgG"]]))
rownames(mylist6[["Hallmark_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["Hallmark_Down_mIgG"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["Reactome_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["Reactome_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["Reactome_Up_mIgG"]] <- add_genes(mylist6[["Reactome_Up_mIgG"]],C2.DB,res_AK23vsmIgG)
mylist6[["Reactome_Down_mIgG"]] <- add_genes(mylist6[["Reactome_Down_mIgG"]],C2.DB,res_AK23vsmIgG)
rownames(mylist6[["Reactome_Up_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist6[["Reactome_Up_mIgG"]]))
rownames(mylist6[["Reactome_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["Reactome_Up_mIgG"]]))
rownames(mylist6[["Reactome_Down_mIgG"]]) <- gsub("REACTOME_","",rownames(mylist6[["Reactome_Down_mIgG"]]))
rownames(mylist6[["Reactome_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["Reactome_Down_mIgG"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["GO_BP_Up_mIgG"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["GO_BP_Down_mIgG"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["GO_BP_Up_mIgG"]] <- add_genes(mylist6[["GO_BP_Up_mIgG"]],C5.BP,res_AK23vsmIgG)
mylist6[["GO_BP_Down_mIgG"]] <- add_genes(mylist6[["GO_BP_Down_mIgG"]],C5.BP,res_AK23vsmIgG)
rownames(mylist6[["GO_BP_Up_mIgG"]]) <- gsub("GOBP_","",rownames(mylist6[["GO_BP_Up_mIgG"]]))
rownames(mylist6[["GO_BP_Up_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["GO_BP_Up_mIgG"]]))
rownames(mylist6[["GO_BP_Down_mIgG"]]) <- gsub("GOBP_","",rownames(mylist6[["GO_BP_Down_mIgG"]]))
rownames(mylist6[["GO_BP_Down_mIgG"]]) <- gsub("_"," ",rownames(mylist6[["GO_BP_Down_mIgG"]]))

WriteXLS(c("mylist6"),"GSEA_AK23vsmIgG_5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################




















