#!/data/sb_service_01/guo/soft/conda/envs/R4/bin Rscript
library(sleuth)
library(tidyverse)
library(biomaRt)
library(dplyr)
library(gridExtra)
library(tximport)
library(org.Hs.eg.db)
library(DESeq2)
#################################

base_dir <- "/Users/buschlab/ncrna/batch2_b20_2"
samples <- dir(file.path(base_dir))
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

#s2c <- read.table( "/data/lied_guo/snsf1.txt", header = TRUE, stringsAsFactors=FALSE)

#######################
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#############
t2g <- t2g[!duplicated(t2g$target_id),]
rownames(t2g) <- t2g$target_id

t2g.red <- t2g[, 2:3]
t2g.red <- unique(t2g.red)
#####################################
s2c <- read.table( "/data/lied_guo/snsf3.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate( s2c, path = kal_dirs)


files <- file.path(kal_dirs, "abundance.h5")
all(file.exists(files))
names(files) <- paste0( s2c$sample )
# 1.2 import files
txi <- tximport(files, type = "kallisto",
                tx2gene = t2g, countsFromAbundance = "lengthScaledTPM",
                txOut=TRUE)
txi.sum <- summarizeToGene(txi, t2g, ignoreTxVersion=TRUE)
head(txi$counts); head(txi.sum$counts)
txi <- txi.sum

# 2. Preprocessing: normalize data
cnts <- txi$counts
sel.rn <- rowSums(cnts) != 0
cnts <- cnts[sel.rn,]

write.csv(cnts, "raw_count.csv") ###   save raw count to a table
libsizes <- colSums(cnts)
size.factor <- libsizes/exp(mean(log(libsizes)))
cnts.norm <- t(t(cnts)/size.factor)
range(cnts.norm)
cnts.norm <- log2(cnts.norm + 8) # add 8 to stabalize the variance at low expression
range(cnts.norm)

head(cnts.norm)
write.csv(cnts, "normalize_count.csv") ###   save normalized count to a table


iqr <- apply(cnts.norm, 1, IQR)
cnts.norm <- cnts.norm[order(iqr,decreasing=T),]

PCA <- prcomp(t(cnts.norm[1:1000,]), scale = T)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4])

#######colvec <- c( rep("cornflowerblue", 6), rep("darkgoldenrod2", 6) )
###shapevec <- c( rep(15, 3), rep(16, 3), rep(17, 3), rep(18, 3) )
shapevec <- c(rep(15,12), rep(16,12), rep(17,12))

colvec <- rep(rep(c("red", "green","blue","yellow"),3),3)


pdf(file="pca.snsf.pdf", height=6, width=6, useDingbats = F)
par(las=1)
plot(dataGG$PC1, dataGG$PC2, xlim=c(-30,40), ylim=c(-20,25),
     xlab=paste("PC1 [", percentVar[1], "%]", sep=""),
     ylab=paste("PC2 [", percentVar[2], "%]", sep=""),
      pch= shapevec, col = colvec, cex = 1.5,
     main="PCA (1000 most variable genes)",
     axes=F)
text(dataGG$PC1, dataGG$PC2, labels = s2c$sample, pos = 1, cex = 0.5)
abline(v=0, h=0, lty=2, col="grey50")
axis(1); axis(2)
legend("topright",
		 col=colvec,pch =shapevec,
       legend = s2c$sample, cex=0.5,
       bty = "n")
dev.off()











