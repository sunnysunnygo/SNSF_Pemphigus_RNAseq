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
sub1 <- subset(s2c, (Stimulus=="PX43" & Time == "5h" ) |  (Stimulus=="PX43" & Time == "10h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## px43 10h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_PX43_10hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_PX43_10hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3
############################################################
sub1 <- subset(s2c, (Stimulus=="PX43" & Time == "5h" ) |  (Stimulus=="PX43" & Time == "24h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## px43 24h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_PX43_24hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_PX43_24hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############3
############################################################
###########
############################################################
sub1 <- subset(s2c, (Stimulus=="hIgG" & Time == "5h" ) |  (Stimulus=="hIgG" & Time == "10h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## hIgG 10h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_hIgG_10hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_hIgG_10hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3
############################################################
sub1 <- subset(s2c, (Stimulus=="hIgG" & Time == "5h" ) |  (Stimulus=="hIgG" & Time == "24h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## hIgG 24h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_hIgG_24hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_hIgG_24hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############3
############################################################
###########
############################################################
sub1 <- subset(s2c, (Stimulus=="AK23" & Time == "5h" ) |  (Stimulus=="AK23" & Time == "10h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## AK23 10h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_AK23_10hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23_10hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3
############################################################
sub1 <- subset(s2c, (Stimulus=="AK23" & Time == "5h" ) |  (Stimulus=="AK23" & Time == "24h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## AK23 24h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_AK23_24hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_AK23_24hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############3
############################################################
###########
############################################################
sub1 <- subset(s2c, (Stimulus=="mIgG" & Time == "5h" ) |  (Stimulus=="mIgG" & Time == "10h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## mIgG 10h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_mIgG_10hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_mIgG_10hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
###############3
############################################################
sub1 <- subset(s2c, (Stimulus=="mIgG" & Time == "5h" ) |  (Stimulus=="mIgG" & Time == "24h" ))
file1 <- sub1$path
file1 <- paste0(file1,"/abundance.h5")
txi <- tximport(file1, type="kallisto", tx2gene=t2g[,c(1,3)],ignoreTxVersion=T)
####
####
######## mIgG 24h vs 5h 
cnts <- as.data.frame(txi$counts)
colnames(cnts) <- sub1$Sample
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
cnts.norm <- round(cnts.norm)
range(cnts.norm)
cnts.norm <- as.data.frame(cnts.norm)
WriteXLS(c("cnts.norm"),"DESeq_mIgG_24hvs5h_normalize_counts_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
WriteXLS(c("sub1"),"DESeq_mIgG_24hvs5h_sample_info_totalRNA.xlsx",row.names=T,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############3
############################################################
##rm(list=ls())
############
#############     px43 10h vs 5h TF
#############################
###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("PX43_10hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_PX43_10hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_PX43_10hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_PX43_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_px43_10hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_PX43_10hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_PX43_10hvs5h_matrix <- ttop_PX43_10hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_PX43_10hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_PX43_10hvs5h_matrix_df <- ttop_PX43_10hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_PX43_10hvs5h_matrix_df[abs(ttop_PX43_10hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_PX43_10hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_PX43_10hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "PX43_10hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("PX43_10hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_PX43_10hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_PX43_10hvs5h[ttop_PX43_10hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_PX43_10hvs5h[ttop_PX43_10hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "PX43_10hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()
################
###############
#####################
####################
############
#############     px43 24h vs 5h TF
#############################

###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("PX43_24hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_PX43_24hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_PX43_24hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_PX43_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_px43_24hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_PX43_24hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_PX43_24hvs5h_matrix <- ttop_PX43_24hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_PX43_24hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_PX43_24hvs5h_matrix_df <- ttop_PX43_24hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_PX43_24hvs5h_matrix_df[abs(ttop_PX43_24hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_PX43_24hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_PX43_24hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "PX43_24hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("PX43_24hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_PX43_24hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_PX43_24hvs5h[ttop_PX43_24hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_PX43_24hvs5h[ttop_PX43_24hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "PX43_24hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
###############3
############################################################
##rm(list=ls())
############
#############     hIgG 10h vs 5h TF
#############################
###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("hIgG_10hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_10hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_10hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_hIgG_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_10hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_hIgG_10hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_hIgG_10hvs5h_matrix <- ttop_hIgG_10hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_hIgG_10hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_hIgG_10hvs5h_matrix_df <- ttop_hIgG_10hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_hIgG_10hvs5h_matrix_df[abs(ttop_hIgG_10hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_hIgG_10hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_hIgG_10hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "hIgG_10hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("hIgG_10hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_hIgG_10hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_hIgG_10hvs5h[ttop_hIgG_10hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_hIgG_10hvs5h[ttop_hIgG_10hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "hIgG_10hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()
################
###############
#####################
####################
############
#############     hIgG 24h vs 5h TF
#############################

###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("hIgG_24hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_24hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_24hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_hIgG_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_hIgG_24hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_hIgG_24hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_hIgG_24hvs5h_matrix <- ttop_hIgG_24hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_hIgG_24hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_hIgG_24hvs5h_matrix_df <- ttop_hIgG_24hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_hIgG_24hvs5h_matrix_df[abs(ttop_hIgG_24hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_hIgG_24hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_hIgG_24hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "hIgG_24hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("hIgG_24hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_hIgG_24hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_hIgG_24hvs5h[ttop_hIgG_24hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_hIgG_24hvs5h[ttop_hIgG_24hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "hIgG_24hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
###############3
############################################################
##rm(list=ls())
############
#############     AK23 10h vs 5h TF
#############################
###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("AK23_10hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_10hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_10hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_AK23_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_10hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_AK23_10hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_AK23_10hvs5h_matrix <- ttop_AK23_10hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_AK23_10hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_AK23_10hvs5h_matrix_df <- ttop_AK23_10hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vsAK23_matrix_df2 <- ttop_AK23_10hvs5h_matrix_df[abs(ttop_AK23_10hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_AK23_10hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vsAK23_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_AK23_10hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "AK23_10hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("AK23_10hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_AK23_10hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_AK23_10hvs5h[ttop_AK23_10hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_AK23_10hvs5h[ttop_AK23_10hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "AK23_10hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()
################
###############
#####################
####################
############
#############     AK23 24h vs 5h TF
#############################

###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("AK23_24hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_24hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_24hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_AK23_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23_24hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_AK23_24hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_AK23_24hvs5h_matrix <- ttop_AK23_24hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_AK23_24hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_AK23_24hvs5h_matrix_df <- ttop_AK23_24hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vsAK23_matrix_df2 <- ttop_AK23_24hvs5h_matrix_df[abs(ttop_AK23_24hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_AK23_24hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vsAK23_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_AK23_24hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "AK23_24hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("AK23_24hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_AK23_24hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_AK23_24hvs5h[ttop_AK23_24hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_AK23_24hvs5h[ttop_AK23_24hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "AK23_24hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################

###############3
############################################################
##rm(list=ls())
############
#############     mIgG 10h vs 5h TF
#############################
###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("mIgG_10hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_mIgGvsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_10hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_10hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_mIgG_10hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_10hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_mIgG_10hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_mIgG_10hvs5h_matrix <- ttop_mIgG_10hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_mIgG_10hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_mIgG_10hvs5h_matrix_df <- ttop_mIgG_10hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vsmIgG_matrix_df2 <- ttop_mIgG_10hvs5h_matrix_df[abs(ttop_mIgG_10hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_mIgG_10hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vsmIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_mIgG_10hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "mIgG_10hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("mIgG_10hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_mIgG_10hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_mIgG_10hvs5h[ttop_mIgG_10hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_mIgG_10hvs5h[ttop_mIgG_10hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "mIgG_10hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()
################
###############
#####################
####################
############
#############     mIgG 24h vs 5h TF
#############################

###03 progeny process
source("/Users/buschlab/Downloads/transcriptutorial/scripts/support_functions.R")
####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(openxlsx)
####
pdf("mIgG_24hvs5h_totalRNA_pathway.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_mIgGvsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_24hvs5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_24hvs5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_mIgG_24hvs5h <- read.xlsx("/Users/buschlab/ncrna/DESeq_mIgG_24hvs5h_totalRNA_deg.xlsx",1,rowNames = F)
names(ttop_mIgG_24hvs5h)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_mIgG_24hvs5h_matrix <- ttop_mIgG_24hvs5h %>% 
  dplyr::select(ID, stat) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  column_to_rownames(var = "ID") %>%
  as.matrix()

#####

PathwayActivity_counts <- progeny(Normalised_counts_matrix, scale=TRUE, 
                                  organism="Human", top = 100)
head(PathwayActivity_counts)

Activity_counts <- as.vector(PathwayActivity_counts)
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))
###progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)],fontsize=14, 
progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_mIgG_24hvs5h_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()

#####
head(PathwayActivity_zscore)

colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))


############
head(PathwayActivity_zscore_df)
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

####################

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

ttop_mIgG_24hvs5h_matrix_df <- ttop_mIgG_24hvs5h_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vsmIgG_matrix_df2 <- ttop_mIgG_24hvs5h_matrix_df[abs(ttop_mIgG_24hvs5h_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_mIgG_24hvs5h_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`JAK-STAT`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vsmIgG_matrix_df, aes(x=GeneID,y=t))

dev.off()
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_mIgG_24hvs5h_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "mIgG_24hvs5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

######################

##### TF by Dorothea 


library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## For the volcano plot (related to support functions)
library(ggrepel)
#######
####
pdf("mIgG_24hvs5h_totalRNA_res_TF.pdf", width = 9, height = 9)
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_mIgG_24hvs5h_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "stat") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

######################
###sub for volcano plot


####
# targets_SPI1 <- regulons$target[regulons$tf == "SPI1"]
# #volcano_nice(as.data.frame(ttop_mIgG_24hvs5h[ttop_mIgG_24hvs5h$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_mIgG_24hvs5h[ttop_mIgG_24hvs5h$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "mIgG_24hvs5h_totalRNA_TFActivity_CARNIVALinput.csv")

#############
tf_activities_counts <- 
  dorothea::run_viper(Normalised_counts_matrix, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

#####
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
tf_activities_counts_filter <- as.data.frame(tf_activities_counts_filter)
#tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################

