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
#####################
###########################

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

#############################

base_dir <- "/Users/buschlab/ncrna/batch2_b20_2"
samples <- dir(file.path(base_dir))
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- read.table("snsf0405.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate( s2c, path = kal_dirs)

rownames(s2c) <- s2c$Sample
print(s2c)
###########
############################################################
sub1 <- subset(s2c, Stimulus=="PX43") 
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####

###
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Time)

##########3

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Time <- relevel(ddsTxi$Time, ref="5h")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
####################

res_px43_10hvs5h <- results(ddsTxi, name  = "Time_10h_vs_5h")
res_px43_24hvs5h <- results(ddsTxi, name  = "Time_24h_vs_5h")

##############
res_px43_10hvs5h_order <- as.data.frame(res_px43_10hvs5h[order(res_px43_10hvs5h$pvalue),])
res_px43_24hvs5h_order <- as.data.frame(res_px43_24hvs5h[order(res_px43_24hvs5h$pvalue),])

WriteXLS(c("res_px43_10hvs5h_order"),"DESeq_px43_10hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("res_px43_24hvs5h_order"),"DESeq_px43_24hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###########################3
pdf("deseq_px43_10hvs5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
# 
# px43_24_vst <- vst(ddsTxi)
# px43_24_vst <- assay(px43_24_vst)
# px43_24_vst <- as.data.frame(px43_24_vst)
# heat_input <- px43_24_vst[rownames(top2),]
# pheatmap(heat_input, scale = "row", cellwidth = 10, cellheight = 10)
####3

res_rlog <- rlog(ddsTxi)
res_rlog <- assay(res_rlog)
res_rlog <- as.data.frame(res_rlog)
#####
res_px43_10hvs5h <- res_px43_10hvs5h[!is.na(res_px43_10hvs5h$padj),]
res_px43_24hvs5h <- res_px43_24hvs5h[!is.na(res_px43_24hvs5h$padj),]
top1 <- res_px43_10hvs5h[res_px43_10hvs5h$padj < 0.05 & abs(res_px43_10hvs5h$log2FoldChange) > 2 ,]
top2 <- res_px43_24hvs5h[res_px43_24hvs5h$padj < 0.05 & abs(res_px43_24hvs5h$log2FoldChange) > 4.5 ,]


heat_input1 <- res_rlog[rownames(top1),][,c(1,2,3,4,5,6)]
heat_input2 <- res_rlog[rownames(top2),][,c(1,2,3,7,8,9)]
pheatmap(heat_input1, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

###################################
library(EnhancedVolcano)
EnhancedVolcano(res_px43_10hvs5h,
                lab = rownames(res_px43_10hvs5h),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'PX43 10hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()

#############################3

pdf("deseq_px43_24hvs5h_heat_volcano.pdf",width = 9,height = 9)
pheatmap(heat_input2, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

EnhancedVolcano(res_px43_24hvs5h,
                lab = rownames(res_px43_24hvs5h),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 3.75,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'PX43 24hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()
#########################################

#########################################
########## hIgG deg
############
sub1 <- subset(s2c, Stimulus=="hIgG") 
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
###
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
####
#WriteXLS(c("cnts.norm"),"DESeq_hIgG_normalize_counts.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#WriteXLS(c("sub1"),"DESeq_hIgG_sample_info.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

#######################
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Time)

##########3

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Time <- relevel(ddsTxi$Time, ref="5h")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
####################

res_hIgG_10hvs5h <- results(ddsTxi, name  = "Time_10h_vs_5h")
res_hIgG_24hvs5h <- results(ddsTxi, name  = "Time_24h_vs_5h")

##############
res_hIgG_10hvs5h_order <- as.data.frame(res_hIgG_10hvs5h[order(res_hIgG_10hvs5h$pvalue),])
res_hIgG_24hvs5h_order <- as.data.frame(res_hIgG_24hvs5h[order(res_hIgG_24hvs5h$pvalue),])

WriteXLS(c("res_hIgG_10hvs5h_order"),"DESeq_hIgG_10hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("res_hIgG_24hvs5h_order"),"DESeq_hIgG_24hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
######################
pdf("deseq_hIgG_10hvs5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
# 
# hIgG_24_vst <- vst(ddsTxi)
# hIgG_24_vst <- assay(hIgG_24_vst)
# hIgG_24_vst <- as.data.frame(hIgG_24_vst)
# heat_input <- hIgG_24_vst[rownames(top2),]
# pheatmap(heat_input, scale = "row", cellwidth = 10, cellheight = 10)
####3

res_rlog <- rlog(ddsTxi)
res_rlog <- assay(res_rlog)
res_rlog <- as.data.frame(res_rlog)
#####
res_hIgG_10hvs5h <- res_hIgG_10hvs5h[!is.na(res_hIgG_10hvs5h$padj),]
res_hIgG_24hvs5h <- res_hIgG_24hvs5h[!is.na(res_hIgG_24hvs5h$padj),]
top1 <- res_hIgG_10hvs5h[res_hIgG_10hvs5h$padj < 0.05 & abs(res_hIgG_10hvs5h$log2FoldChange) >2 ,]
top2 <- res_hIgG_24hvs5h[res_hIgG_24hvs5h$padj < 0.05 & abs(res_hIgG_24hvs5h$log2FoldChange) >4.5 ,]
heat_input1 <- res_rlog[rownames(top1),][,c(1,2,3,4,5,6)]
heat_input2 <- res_rlog[rownames(top2),][,c(1,2,3,7,8,9)]
pheatmap(heat_input1, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
##pheatmap(heat_input1[, order(colnames(heat_input1))], scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)


##########
library(EnhancedVolcano)
EnhancedVolcano(res_hIgG_10hvs5h,
                lab = rownames(res_hIgG_10hvs5h),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'hIgG 10hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 30,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


pdf("deseq_hIgG_24hvs5h_heat_volcano.pdf",width = 9,height = 9)
pheatmap(heat_input2, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

EnhancedVolcano(res_hIgG_24hvs5h,
                lab = rownames(res_hIgG_24hvs5h),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'hIgG 24hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()
#########################################

#########################################
#######33  AK23 deg
######
sub1 <- subset(s2c, Stimulus=="AK23") 
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
###
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
####
WriteXLS(c("cnts.norm"),"DESeq_AK23_normalize_counts.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23_sample_info.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

#######################
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Time)

##########3

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Time <- relevel(ddsTxi$Time, ref="5h")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
####################

res_AK23_10hvs5h <- results(ddsTxi, name  = "Time_10h_vs_5h")
res_AK23_24hvs5h <- results(ddsTxi, name  = "Time_24h_vs_5h")

##############
res_AK23_10hvs5h_order <- as.data.frame(res_AK23_10hvs5h[order(res_AK23_10hvs5h$pvalue),])
res_AK23_24hvs5h_order <- as.data.frame(res_AK23_24hvs5h[order(res_AK23_24hvs5h$pvalue),])

WriteXLS(c("res_AK23_10hvs5h_order"),"DESeq_AK23_10hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("res_AK23_24hvs5h_order"),"DESeq_AK23_24hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
######################
pdf("deseq_AK23_10hvs5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
# 
# AK23_24_vst <- vst(ddsTxi)
# AK23_24_vst <- assay(AK23_24_vst)
# AK23_24_vst <- as.data.frame(AK23_24_vst)
# heat_input <- AK23_24_vst[rownames(top2),]
# pheatmap(heat_input, scale = "row", cellwidth = 10, cellheight = 10)
####3

res_rlog <- rlog(ddsTxi)
res_rlog <- assay(res_rlog)
res_rlog <- as.data.frame(res_rlog)
#####
res_AK23_10hvs5h <- res_AK23_10hvs5h[!is.na(res_AK23_10hvs5h$padj),]
res_AK23_24hvs5h <- res_AK23_24hvs5h[!is.na(res_AK23_24hvs5h$padj),]
top1 <- res_AK23_10hvs5h[res_AK23_10hvs5h$padj < 0.05 & abs(res_AK23_10hvs5h$log2FoldChange) > 2 ,]
top2 <- res_AK23_24hvs5h[res_AK23_24hvs5h$padj < 0.05 & abs(res_AK23_24hvs5h$log2FoldChange) >4.5 ,]
heat_input1 <- res_rlog[rownames(top1),][,c(1,2,3,4,5,6)]
heat_input2 <- res_rlog[rownames(top2),][,c(1,2,3,7,8,9)]
pheatmap(heat_input1, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
##pheatmap(heat_input1[, order(colnames(heat_input1))], scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)


##########
library(EnhancedVolcano)
EnhancedVolcano(res_AK23_10hvs5h,
                lab = rownames(res_AK23_10hvs5h),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'AK23 10hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


pdf("deseq_AK23_24hvs5h_heat_volcano.pdf",width = 9,height = 9)
pheatmap(heat_input2, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

EnhancedVolcano(res_AK23_24hvs5h,
                lab = rownames(res_AK23_24hvs5h),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'AK23 24hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()
#########################################

#########################################
#######33  mIgG deg
######
sub1 <- subset(s2c, Stimulus=="mIgG") 
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
###
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
####
#WriteXLS(c("cnts.norm"),"DESeq_mIgG_normalize_counts.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#WriteXLS(c("sub1"),"DESeq_mIgG_sample_info.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

#######################
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sub1,
                                   design = ~Batch+Time)

##########3

keep <- rowSums(counts(ddsTxi)) >= 30
ddsTxi <- ddsTxi[keep,]

ddsTxi$Time <- relevel(ddsTxi$Time, ref="5h")
ddsTxi <- DESeq(ddsTxi,test="LRT", reduced=~Batch)
resultsNames(ddsTxi)
####################

res_mIgG_10hvs5h <- results(ddsTxi, name  = "Time_10h_vs_5h")
res_mIgG_24hvs5h <- results(ddsTxi, name  = "Time_24h_vs_5h")

##############
res_mIgG_10hvs5h_order <- as.data.frame(res_mIgG_10hvs5h[order(res_mIgG_10hvs5h$pvalue),])
res_mIgG_24hvs5h_order <- as.data.frame(res_mIgG_24hvs5h[order(res_mIgG_24hvs5h$pvalue),])

WriteXLS(c("res_mIgG_10hvs5h_order"),"DESeq_mIgG_10hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("res_mIgG_24hvs5h_order"),"DESeq_mIgG_24hvs5h_totalRNA_deg.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
######################
pdf("deseq_mIgG_10hvs5h_heat_volcano.pdf",width = 9,height = 9)

library(pheatmap)
# 
# mIgG_24_vst <- vst(ddsTxi)
# mIgG_24_vst <- assay(mIgG_24_vst)
# mIgG_24_vst <- as.data.frame(mIgG_24_vst)
# heat_input <- mIgG_24_vst[rownames(top2),]
# pheatmap(heat_input, scale = "row", cellwidth = 10, cellheight = 10)
####3

res_rlog <- rlog(ddsTxi)
res_rlog <- assay(res_rlog)
res_rlog <- as.data.frame(res_rlog)
#####
res_mIgG_10hvs5h <- res_mIgG_10hvs5h[!is.na(res_mIgG_10hvs5h$padj),]
res_mIgG_24hvs5h <- res_mIgG_24hvs5h[!is.na(res_mIgG_24hvs5h$padj),]
top1 <- res_mIgG_10hvs5h[res_mIgG_10hvs5h$padj < 0.05 & abs(res_mIgG_10hvs5h$log2FoldChange) > 2 ,]
top2 <- res_mIgG_24hvs5h[res_mIgG_24hvs5h$padj < 0.05 & abs(res_mIgG_24hvs5h$log2FoldChange) >4.5 ,]
heat_input1 <- res_rlog[rownames(top1),][,c(1,2,3,4,5,6)]
heat_input2 <- res_rlog[rownames(top2),][,c(1,2,3,7,8,9)]
pheatmap(heat_input1, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)
##pheatmap(heat_input1[, order(colnames(heat_input1))], scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)


##########
library(EnhancedVolcano)
EnhancedVolcano(res_mIgG_10hvs5h,
                lab = rownames(res_mIgG_10hvs5h),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'mIgG 10hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()


pdf("deseq_mIgG_24hvs5h_heat_volcano.pdf",width = 9,height = 9)
pheatmap(heat_input2, scale = "row", cellwidth = 10, cellheight = 10,cluster_cols = F)

EnhancedVolcano(res_mIgG_24hvs5h,
                lab = rownames(res_mIgG_24hvs5h),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10, 10),
                #ylim = c(0, 12),
                legendPosition = 'top', arrowheads = F,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                title = 'mIgG 24hvs5h',
                subtitle = ' ',
                xlab = "log2FC",
                max.overlaps = 20,
                legendLabSize = 9,
                legendIconSize = 4,)

dev.off()
#########################################

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


##########################
library(gage)
###################PX43vshIgG
#######################################

###################px43_10hvs5h
#######################################
mylist1 <- list() 
res_px43_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_px43_10hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_px43_10hvs5h$log2FoldChange,row.names=rownames(res_px43_10hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Hallmark_Up"]] <- add_genes(mylist1[["Hallmark_Up"]],HALLMARK.DB,res_px43_10hvs5h)
mylist1[["Hallmark_Down"]] <- add_genes(mylist1[["Hallmark_Down"]],HALLMARK.DB,res_px43_10hvs5h)
rownames(mylist1[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Up"]]))
rownames(mylist1[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Up"]]))
rownames(mylist1[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Down"]]))
rownames(mylist1[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Reactome_Up"]] <- add_genes(mylist1[["Reactome_Up"]],C2.DB,res_px43_10hvs5h)
mylist1[["Reactome_Down"]] <- add_genes(mylist1[["Reactome_Down"]],C2.DB,res_px43_10hvs5h)
rownames(mylist1[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Up"]]))
rownames(mylist1[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Up"]]))
rownames(mylist1[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Down"]]))
rownames(mylist1[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["GO_BP_Up"]] <- add_genes(mylist1[["GO_BP_Up"]],C5.BP,res_px43_10hvs5h)
mylist1[["GO_BP_Down"]] <- add_genes(mylist1[["GO_BP_Down"]],C5.BP,res_px43_10hvs5h)
rownames(mylist1[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Up"]]))
rownames(mylist1[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Up"]]))
rownames(mylist1[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Down"]]))
rownames(mylist1[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Down"]]))

WriteXLS(c("mylist1"),"GSEA_px43_10hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################px43_24hvs5h
#######################################
mylist1 <- list() 
res_px43_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_px43_24hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_px43_24hvs5h$log2FoldChange,row.names=rownames(res_px43_24hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Hallmark_Up"]] <- add_genes(mylist1[["Hallmark_Up"]],HALLMARK.DB,res_px43_24hvs5h)
mylist1[["Hallmark_Down"]] <- add_genes(mylist1[["Hallmark_Down"]],HALLMARK.DB,res_px43_24hvs5h)
rownames(mylist1[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Up"]]))
rownames(mylist1[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Up"]]))
rownames(mylist1[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist1[["Hallmark_Down"]]))
rownames(mylist1[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist1[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["Reactome_Up"]] <- add_genes(mylist1[["Reactome_Up"]],C2.DB,res_px43_24hvs5h)
mylist1[["Reactome_Down"]] <- add_genes(mylist1[["Reactome_Down"]],C2.DB,res_px43_24hvs5h)
rownames(mylist1[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Up"]]))
rownames(mylist1[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Up"]]))
rownames(mylist1[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist1[["Reactome_Down"]]))
rownames(mylist1[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist1[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist1[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist1[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist1[["GO_BP_Up"]] <- add_genes(mylist1[["GO_BP_Up"]],C5.BP,res_px43_24hvs5h)
mylist1[["GO_BP_Down"]] <- add_genes(mylist1[["GO_BP_Down"]],C5.BP,res_px43_24hvs5h)
rownames(mylist1[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Up"]]))
rownames(mylist1[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Up"]]))
rownames(mylist1[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist1[["GO_BP_Down"]]))
rownames(mylist1[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist1[["GO_BP_Down"]]))

WriteXLS(c("mylist1"),"GSEA_px43_24hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################hIgG_10hvs5h
#######################################
mylist3 <- list() 
res_hIgG_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_10hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_hIgG_10hvs5h$log2FoldChange,row.names=rownames(res_hIgG_10hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["Hallmark_Up"]] <- add_genes(mylist3[["Hallmark_Up"]],HALLMARK.DB,res_hIgG_10hvs5h)
mylist3[["Hallmark_Down"]] <- add_genes(mylist3[["Hallmark_Down"]],HALLMARK.DB,res_hIgG_10hvs5h)
rownames(mylist3[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist3[["Hallmark_Up"]]))
rownames(mylist3[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist3[["Hallmark_Up"]]))
rownames(mylist3[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist3[["Hallmark_Down"]]))
rownames(mylist3[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist3[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["Reactome_Up"]] <- add_genes(mylist3[["Reactome_Up"]],C2.DB,res_hIgG_10hvs5h)
mylist3[["Reactome_Down"]] <- add_genes(mylist3[["Reactome_Down"]],C2.DB,res_hIgG_10hvs5h)
rownames(mylist3[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist3[["Reactome_Up"]]))
rownames(mylist3[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist3[["Reactome_Up"]]))
rownames(mylist3[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist3[["Reactome_Down"]]))
rownames(mylist3[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist3[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist3[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist3[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist3[["GO_BP_Up"]] <- add_genes(mylist3[["GO_BP_Up"]],C5.BP,res_hIgG_10hvs5h)
mylist3[["GO_BP_Down"]] <- add_genes(mylist3[["GO_BP_Down"]],C5.BP,res_hIgG_10hvs5h)
rownames(mylist3[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist3[["GO_BP_Up"]]))
rownames(mylist3[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist3[["GO_BP_Up"]]))
rownames(mylist3[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist3[["GO_BP_Down"]]))
rownames(mylist3[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist3[["GO_BP_Down"]]))

WriteXLS(c("mylist3"),"GSEA_hIgG_10hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################hIgG_24hvs5h
#######################################
mylist4 <- list() 
res_hIgG_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_24hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_hIgG_24hvs5h$log2FoldChange,row.names=rownames(res_hIgG_24hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["Hallmark_Up"]] <- add_genes(mylist4[["Hallmark_Up"]],HALLMARK.DB,res_hIgG_24hvs5h)
mylist4[["Hallmark_Down"]] <- add_genes(mylist4[["Hallmark_Down"]],HALLMARK.DB,res_hIgG_24hvs5h)
rownames(mylist4[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist4[["Hallmark_Up"]]))
rownames(mylist4[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist4[["Hallmark_Up"]]))
rownames(mylist4[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist4[["Hallmark_Down"]]))
rownames(mylist4[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist4[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["Reactome_Up"]] <- add_genes(mylist4[["Reactome_Up"]],C2.DB,res_hIgG_24hvs5h)
mylist4[["Reactome_Down"]] <- add_genes(mylist4[["Reactome_Down"]],C2.DB,res_hIgG_24hvs5h)
rownames(mylist4[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist4[["Reactome_Up"]]))
rownames(mylist4[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist4[["Reactome_Up"]]))
rownames(mylist4[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist4[["Reactome_Down"]]))
rownames(mylist4[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist4[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist4[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist4[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist4[["GO_BP_Up"]] <- add_genes(mylist4[["GO_BP_Up"]],C5.BP,res_hIgG_24hvs5h)
mylist4[["GO_BP_Down"]] <- add_genes(mylist4[["GO_BP_Down"]],C5.BP,res_hIgG_24hvs5h)
rownames(mylist4[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist4[["GO_BP_Up"]]))
rownames(mylist4[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist4[["GO_BP_Up"]]))
rownames(mylist4[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist4[["GO_BP_Down"]]))
rownames(mylist4[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist4[["GO_BP_Down"]]))

WriteXLS(c("mylist4"),"GSEA_hIgG_24hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################AK23_10hvs5h
#######################################
mylist5 <- list() 
res_AK23_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_10hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_AK23_10hvs5h$log2FoldChange,row.names=rownames(res_AK23_10hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["Hallmark_Up"]] <- add_genes(mylist5[["Hallmark_Up"]],HALLMARK.DB,res_AK23_10hvs5h)
mylist5[["Hallmark_Down"]] <- add_genes(mylist5[["Hallmark_Down"]],HALLMARK.DB,res_AK23_10hvs5h)
rownames(mylist5[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist5[["Hallmark_Up"]]))
rownames(mylist5[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist5[["Hallmark_Up"]]))
rownames(mylist5[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist5[["Hallmark_Down"]]))
rownames(mylist5[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist5[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["Reactome_Up"]] <- add_genes(mylist5[["Reactome_Up"]],C2.DB,res_AK23_10hvs5h)
mylist5[["Reactome_Down"]] <- add_genes(mylist5[["Reactome_Down"]],C2.DB,res_AK23_10hvs5h)
rownames(mylist5[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist5[["Reactome_Up"]]))
rownames(mylist5[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist5[["Reactome_Up"]]))
rownames(mylist5[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist5[["Reactome_Down"]]))
rownames(mylist5[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist5[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist5[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist5[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist5[["GO_BP_Up"]] <- add_genes(mylist5[["GO_BP_Up"]],C5.BP,res_AK23_10hvs5h)
mylist5[["GO_BP_Down"]] <- add_genes(mylist5[["GO_BP_Down"]],C5.BP,res_AK23_10hvs5h)
rownames(mylist5[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist5[["GO_BP_Up"]]))
rownames(mylist5[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist5[["GO_BP_Up"]]))
rownames(mylist5[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist5[["GO_BP_Down"]]))
rownames(mylist5[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist5[["GO_BP_Down"]]))

WriteXLS(c("mylist5"),"GSEA_AK23_10hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################AK23_24hvs5h
#######################################
mylist6 <- list() 
res_AK23_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_24hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_AK23_24hvs5h$log2FoldChange,row.names=rownames(res_AK23_24hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["Hallmark_Up"]] <- add_genes(mylist6[["Hallmark_Up"]],HALLMARK.DB,res_AK23_24hvs5h)
mylist6[["Hallmark_Down"]] <- add_genes(mylist6[["Hallmark_Down"]],HALLMARK.DB,res_AK23_24hvs5h)
rownames(mylist6[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist6[["Hallmark_Up"]]))
rownames(mylist6[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist6[["Hallmark_Up"]]))
rownames(mylist6[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist6[["Hallmark_Down"]]))
rownames(mylist6[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist6[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["Reactome_Up"]] <- add_genes(mylist6[["Reactome_Up"]],C2.DB,res_AK23_24hvs5h)
mylist6[["Reactome_Down"]] <- add_genes(mylist6[["Reactome_Down"]],C2.DB,res_AK23_24hvs5h)
rownames(mylist6[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist6[["Reactome_Up"]]))
rownames(mylist6[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist6[["Reactome_Up"]]))
rownames(mylist6[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist6[["Reactome_Down"]]))
rownames(mylist6[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist6[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist6[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist6[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist6[["GO_BP_Up"]] <- add_genes(mylist6[["GO_BP_Up"]],C5.BP,res_AK23_24hvs5h)
mylist6[["GO_BP_Down"]] <- add_genes(mylist6[["GO_BP_Down"]],C5.BP,res_AK23_24hvs5h)
rownames(mylist6[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist6[["GO_BP_Up"]]))
rownames(mylist6[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist6[["GO_BP_Up"]]))
rownames(mylist6[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist6[["GO_BP_Down"]]))
rownames(mylist6[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist6[["GO_BP_Down"]]))

WriteXLS(c("mylist6"),"GSEA_AK23_24hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################
###################mIgG_10hvs5h
#######################################
mylist7 <- list() 
res_mIgG_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_10hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_mIgG_10hvs5h$log2FoldChange,row.names=rownames(res_mIgG_10hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist7[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist7[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist7[["Hallmark_Up"]] <- add_genes(mylist7[["Hallmark_Up"]],HALLMARK.DB,res_mIgG_10hvs5h)
mylist7[["Hallmark_Down"]] <- add_genes(mylist7[["Hallmark_Down"]],HALLMARK.DB,res_mIgG_10hvs5h)
rownames(mylist7[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist7[["Hallmark_Up"]]))
rownames(mylist7[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist7[["Hallmark_Up"]]))
rownames(mylist7[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist7[["Hallmark_Down"]]))
rownames(mylist7[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist7[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist7[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist7[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist7[["Reactome_Up"]] <- add_genes(mylist7[["Reactome_Up"]],C2.DB,res_mIgG_10hvs5h)
mylist7[["Reactome_Down"]] <- add_genes(mylist7[["Reactome_Down"]],C2.DB,res_mIgG_10hvs5h)
rownames(mylist7[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist7[["Reactome_Up"]]))
rownames(mylist7[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist7[["Reactome_Up"]]))
rownames(mylist7[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist7[["Reactome_Down"]]))
rownames(mylist7[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist7[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist7[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist7[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist7[["GO_BP_Up"]] <- add_genes(mylist7[["GO_BP_Up"]],C5.BP,res_mIgG_10hvs5h)
mylist7[["GO_BP_Down"]] <- add_genes(mylist7[["GO_BP_Down"]],C5.BP,res_mIgG_10hvs5h)
rownames(mylist7[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist7[["GO_BP_Up"]]))
rownames(mylist7[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist7[["GO_BP_Up"]]))
rownames(mylist7[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist7[["GO_BP_Down"]]))
rownames(mylist7[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist7[["GO_BP_Down"]]))

WriteXLS(c("mylist7"),"GSEA_mIgG_10hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################

###################mIgG_24hvs5h
#######################################
mylist8 <- list() 
res_mIgG_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_24hvs5h_totalRNA_deg.xlsx",1,rowNames = T) ####### 

tmp <- data.frame(B=res_mIgG_24hvs5h$log2FoldChange,row.names=rownames(res_mIgG_24hvs5h))
# Hallmark
resgo <- gage(tmp,gsets=HALLMARK.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist8[["Hallmark_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist8[["Hallmark_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist8[["Hallmark_Up"]] <- add_genes(mylist8[["Hallmark_Up"]],HALLMARK.DB,res_mIgG_24hvs5h)
mylist8[["Hallmark_Down"]] <- add_genes(mylist8[["Hallmark_Down"]],HALLMARK.DB,res_mIgG_24hvs5h)
rownames(mylist8[["Hallmark_Up"]]) <- gsub("HALLMARK_","",rownames(mylist8[["Hallmark_Up"]]))
rownames(mylist8[["Hallmark_Up"]]) <- gsub("_"," ",rownames(mylist8[["Hallmark_Up"]]))
rownames(mylist8[["Hallmark_Down"]]) <- gsub("HALLMARK_","",rownames(mylist8[["Hallmark_Down"]]))
rownames(mylist8[["Hallmark_Down"]]) <- gsub("_"," ",rownames(mylist8[["Hallmark_Down"]]))

# Reactome
resgo <- gage(tmp,gsets=C2.DB)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist8[["Reactome_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist8[["Reactome_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist8[["Reactome_Up"]] <- add_genes(mylist8[["Reactome_Up"]],C2.DB,res_mIgG_24hvs5h)
mylist8[["Reactome_Down"]] <- add_genes(mylist8[["Reactome_Down"]],C2.DB,res_mIgG_24hvs5h)
rownames(mylist8[["Reactome_Up"]]) <- gsub("REACTOME_","",rownames(mylist8[["Reactome_Up"]]))
rownames(mylist8[["Reactome_Up"]]) <- gsub("_"," ",rownames(mylist8[["Reactome_Up"]]))
rownames(mylist8[["Reactome_Down"]]) <- gsub("REACTOME_","",rownames(mylist8[["Reactome_Down"]]))
rownames(mylist8[["Reactome_Down"]]) <- gsub("_"," ",rownames(mylist8[["Reactome_Down"]]))

# GO
resgo <- gage(tmp,gsets=C5.BP)
sig_go <-sigGeneSet(resgo,cutoff=0.1,qpval = "p.val")
mylist8[["GO_BP_Up"]] <- subset(as.data.frame(sig_go$greater),select=-c(exp1))
mylist8[["GO_BP_Down"]] <- subset(as.data.frame(sig_go$less),select=-c(exp1))
# Add Genes
mylist8[["GO_BP_Up"]] <- add_genes(mylist8[["GO_BP_Up"]],C5.BP,res_mIgG_24hvs5h)
mylist8[["GO_BP_Down"]] <- add_genes(mylist8[["GO_BP_Down"]],C5.BP,res_mIgG_24hvs5h)
rownames(mylist8[["GO_BP_Up"]]) <- gsub("GOBP_","",rownames(mylist8[["GO_BP_Up"]]))
rownames(mylist8[["GO_BP_Up"]]) <- gsub("_"," ",rownames(mylist8[["GO_BP_Up"]]))
rownames(mylist8[["GO_BP_Down"]]) <- gsub("GOBP_","",rownames(mylist8[["GO_BP_Down"]]))
rownames(mylist8[["GO_BP_Down"]]) <- gsub("_"," ",rownames(mylist8[["GO_BP_Down"]]))

WriteXLS(c("mylist8"),"GSEA_mIgG_24hvs5h_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###########################
#######################







