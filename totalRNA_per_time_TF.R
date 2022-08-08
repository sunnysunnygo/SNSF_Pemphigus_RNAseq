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
pdf("AK23vsmIgG_24h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_24h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_24h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_24h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_AK23vsmIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_AK23vsmIgG_matrix <- ttop_AK23vsmIgG %>% 
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
progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_AK23vsmIgG_matrix, 
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

ttop_AK23vsmIgG_matrix_df <- ttop_AK23vsmIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_AK23vsmIgG_matrix_df[abs(ttop_AK23vsmIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_AK23vsmIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`TNFa`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_AK23vsmIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "AK23vsmIgG_24h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_AK23vsmIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "AK23vsmIgG_24h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############


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
pdf("AK23vsmIgG_10h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_10h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_10h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_10h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_AK23vsmIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_AK23vsmIgG_matrix <- ttop_AK23vsmIgG %>% 
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

progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_AK23vsmIgG_matrix, 
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

ttop_AK23vsmIgG_matrix_df <- ttop_AK23vsmIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_AK23vsmIgG_matrix_df[abs(ttop_AK23vsmIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_AK23vsmIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`Androgen`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_AK23vsmIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "AK23vsmIgG_10h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_AK23vsmIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "AK23vsmIgG_10h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]

dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############
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
pdf("AK23vsmIgG_5h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_AK23vsmIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_AK23vsmIgG_5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_AK23vsmIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_AK23vsmIgG_5h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_AK23vsmIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_AK23vsmIgG_matrix <- ttop_AK23vsmIgG %>% 
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

progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_AK23vsmIgG_matrix, 
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

ttop_AK23vsmIgG_matrix_df <- ttop_AK23vsmIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_AK23vsmIgG_matrix_df[abs(ttop_AK23vsmIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_AK23vsmIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`TNFa`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_AK23vsmIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "AK23vsmIgG_5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_AK23vsmIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_AK23vsmIgG[ttop_AK23vsmIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "AK23vsmIgG_5h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############
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
pdf("PX43vshIgG_5h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_PX43vshIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_5h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_5h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_PX43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_5h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_PX43vshIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_PX43vshIgG_matrix <- ttop_PX43vshIgG %>% 
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

progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_PX43vshIgG_matrix, 
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

ttop_PX43vshIgG_matrix_df <- ttop_PX43vshIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_PX43vshIgG_matrix_df[abs(ttop_PX43vshIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_PX43vshIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`Hypoxia`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_PX43vshIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "PX43vshIgG_5h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_PX43vshIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "PX43vshIgG_5h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############
#############3
#################

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
pdf("PX43vshIgG_10h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_PX43vshIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_10h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_10h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_PX43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_10h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_PX43vshIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_PX43vshIgG_matrix <- ttop_PX43vshIgG %>% 
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

progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_PX43vshIgG_matrix, 
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

ttop_PX43vshIgG_matrix_df <- ttop_PX43vshIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_PX43vshIgG_matrix_df[abs(ttop_PX43vshIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_PX43vshIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`TNFa`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_PX43vshIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "PX43vshIgG_10h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_PX43vshIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "PX43vshIgG_10h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############
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
pdf("PX43vshIgG_24h_totalRNA_tf.pdf", width = 9, height = 9)
####Normalised_counts <- read_csv("/Users/buschlab/snsf/DESeq_PX43vshIgG_general_normalize_counts.xlsx")
Normalised_counts <- read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_24h_normalize_counts_totalRNA.xlsx",1,rowNames = F)
head(Normalised_counts)
names(Normalised_counts)[1] <- "gene"

Experimental_design <-  read.xlsx("/Users/buschlab/ncrna/DESeq_PX43vshIgG_24h_sample_info_totalRNA.xlsx",1,rowNames = F)

ttop_PX43vshIgG <- read.xlsx("/Users/buschlab/ncrna/per_timepoint/DESeq_PX43vshIgG_24h_deg_totalRNA.xlsx",1,rowNames = F)
names(ttop_PX43vshIgG)[1] <- "ID"
####

Normalised_counts_matrix <- Normalised_counts %>% 
  dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

ttop_PX43vshIgG_matrix <- ttop_PX43vshIgG %>% 
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

progeny_hmap <- pheatmap(t(PathwayActivity_counts)[, c(1,3,5,2,4,6)],fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA,cluster_cols=F )
######

PathwayActivity_zscore <- progeny(ttop_PX43vshIgG_matrix, 
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

ttop_PX43vshIgG_matrix_df <- ttop_PX43vshIgG_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")
#ttop_PX4_3vshIgG_matrix_df2 <- ttop_PX43vshIgG_matrix_df[abs(ttop_PX43vshIgG_matrix_df$stat) < 2,]
scat_plots <- progeny::progenyScatter(df = ttop_PX43vshIgG_matrix_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "t_values", verbose = FALSE)

plot(scat_plots[[1]]$`TNFa`) ### ylim not works, ????????? 
#ggplot(ttop_PX4_3vshIgG_matrix_df, aes(x=GeneID,y=t))
#########################3
PathwayActivity_CARNIVALinput <- progeny(ttop_PX43vshIgG_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

head(PathwayActivity_CARNIVALinput)

write_csv(PathwayActivity_CARNIVALinput, 
          "PX43vshIgG_24h_totalRNA_PathwayActivity_CARNIVALinput.csv")

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

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(ttop_PX43vshIgG_matrix, regulons,
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
# #volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_SPI1,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

# targets_JUN <- regulons$target[regulons$tf == "JUN"]
# 
# volcano_nice(as.data.frame(ttop_PX43vshIgG[ttop_PX43vshIgG$ID %in% targets_JUN,]), 
#              FCIndex = 2, pValIndex = 4, IDIndex = 1,nlabels = 20, label = T, 
#              straight = FALSE) 

##################3
tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 

write_csv(tf_activities_CARNIVALinput, "PX43vshIgG_24h_totalRNA_TFActivity_CARNIVALinput.csv")

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
tf_activities_counts_filter <- tf_activities_counts_filter[, c(1,3,5,2,4,6)]
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA,scale = "row",cellwidth = 10, cellheight = 10,cluster_cols = F)



dev.off()

####################
############





























