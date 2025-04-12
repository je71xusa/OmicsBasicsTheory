## sample script for tutorial using Mov10 dataset
## adapted from HBC tutorial https://hbctraining.github.io/DGE_workshop/schedule/1.5-day.html

library(ggplot2)
library(ggrepel)
library(DESeq2)
library(magrittr)
library(pheatmap)
library(fdrtool)
library(sva)
#### 1. Loading data and  first look ####
data <- read.table(paste0("data/Mov10_full_counts.txt"), header = T, row.names = 1)
meta <- read.table(paste0("meta/Mov10_full_meta.txt"), header = T, row.names = 1)

View(data)
View(meta)

# What is the experimental design?
# Answer:

#### 2. Distribution of counts and the problem of low counts ####
mov10_kd <- rowMeans(data[,c(1,2)]) # mean of expression in Mov10 Knockdown
mov10_oe <- rowMeans(data[,c(3,4,5)]) # mean of expression in Mov10 Overexpression
control <- rowMeans(data[,c(6,7,8)]) # mean of expression in Control

mov10_df <- data.frame(
  gene = rownames(data),
  mov10_kd = mov10_kd,
  mov10_oe = mov10_oe,
  control = control
)

mov10_df$diffOEctrl <- mov10_df$mov10_oe - mov10_df$control # calculating difference in expression between Overexpression and control

head(mov10_df[order(mov10_df$diffOEctrl),]) 
tail(mov10_df[order(mov10_df$diffOEctrl),])

# Why can we not just rank the genes according to difference in expression?

# Answer

## plotting distribution of counts for one sample
hist(data$Mov10_oe_1, breaks = 200)
hist(data$Mov10_oe_1, breaks = 200000, xlim = c(0,100))

## look at low counts here (genes were chose already)
mov10_df[which(mov10_df$gene %in% c("CCL5", "CCDC60")),]

## Why normalization
colSums(mov10_df[,-c(1,5)]) # check number of reads in each condition
colSums(data) # check number of reads in each sample

## We will do normalization and low count filtering etc. all in DESeq2!

#### 3.Deseq2 pipeline basic ####
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

## Low count filtering
keep <- rowSums(counts(dds)) >= 10 # everything above 10 across samples
dds <- dds[keep,]

## Normalization
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

norm_counts <- counts(dds, normalized=TRUE)
norm_counts[1,1]

scPlotdf <- data.frame(sizeFactors(dds),colSums(counts(dds)))
colnames(scPlotdf) <- c("sizeF","seqDepth")
ggplot(scPlotdf,aes(x=sizeF,y=seqDepth)) + geom_point() + 
  geom_text(label = rownames(scPlotdf), nudge_x = 0.15, check_overlap = T)

## data visualization
vs <- vst(dds)
plotDispEsts(dds) 
#### 4. Exploration for batch effects ####
vs_meta <- as.data.frame(colData(vs))
colnames(vs_meta)
## pca
plotPCA(vs, intgroup = "sampletype", ) + geom_label_repel(aes(label = name))
## hierachical clustering
vs_mat <- assay(vs)
vs_cor <- cor(vs_mat)
pheatmap(vs_cor)

## for example of a batch effect we will add some fake ones
vs$batch1 <- "A"
vs$batch1[which(vs$sampletype %in% c("control", "MOV10_knockdown"))] <- "B"

vs$batch2 <- "A"
vs$batch2[c(1,3,6)] <- "B"

# This data set is very clean so in general I would not do any batch correction

#### 5. SVA ####
norm.cts <- counts(dds, normalized=TRUE)
mm <- model.matrix(~ sampletype, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))
norm.cts <- norm.cts[rowSums(norm.cts) > 10,]
fit <- svaseq(norm.cts, mod=mm, mod0=mm0)
stripchart(fit$sv[,1] ~ dds$sampletype,vertical = TRUE,main="SV1")
cleaningY = function(y, mod, svaobj) {
  X = cbind(mod, svaobj$sv)
  Hat = solve(t(X)%*%X)%*%t(X)
  beta = (Hat%*%t(y))
  P = ncol(mod)
  cleany = y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
  return(cleany)
}
cleaned <- cleaningY(norm.cts,mod=mm,fit)
pca = prcomp(t( cleaned ), scale=FALSE)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
ggplot(data=d, aes(x = PC1, y = PC2, label = rownames(d))) + 
  geom_point(aes(colour = rownames(d))) + geom_text(hjust=0,vjust=0.5)
ggplot(data=d, aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = rownames(d))) 

dds$SV1 <- fit$sv[,1]
dds$SV2 <- fit$sv[,2]

#### 5. Run DESEQ2 ####
## define your design
design(dds) <- ~ SV1 + SV2 + sampletype
## run the test using Wald, this will do normalization etc. all by itself
dds <- DESeq(dds, test = "Wald")

## classic running and comparing between including SVs and no SVs
res <- results(dds)
resultsNames(dds)
res_tbl_oneSVs <- as.data.frame(results(dds, 
                          contrast = c("sampletype", "MOV10_overexpression", "control")))
res_tbl_twoSVs_sig <- res_tbl_twoSVs[which(res_tbl_twoSVs$padj < 0.05),]
res_tbl_noSVs_sig <- res_tbl_noSVs[which(res_tbl_noSVs$padj < 0.05),]
res_tbl_oneSVs_sig <- res_tbl_oneSVs[which(res_tbl_oneSVs$padj < 0.05),]

# running for loop for later with differing SV numbers
designs <- c("twoSVs", "oneSVs", "noSVs")
res_list <- list()
for(des in designs){
  if(des == "twoSVs"){
    design(dds) <- ~ SV1 + SV2 + sampletype
  }else if(des == "oneSVs"){
    design(dds) <- ~ SV1 + sampletype
  }else{
    design(dds) <- ~ sampletype
  }
  dds <- DESeq(dds, test = "Wald")
  res_tbl <- as.data.frame(results(dds, 
                                          contrast = c("sampletype", "MOV10_overexpression", "control")))
  res_tbl_sig <- res_tbl[which(res_tbl$padj < 0.05),]
  res_list[["allres"]][[des]] <- res_tbl
  res_list[["sigres"]][[des]] <- res_tbl_sig
}

## explaining pvalue histograms, log2Fold shrinkage
hist(res_tbl_twoSVs$pvalue, breaks=200)
# after deseq2 and shrinkage
plotMA(res)
plotMA(allres2)
## fdrtool

## Interpretation of output and visualization

# volcano plot
allres <- res_list$allres$twoSVs
sigres <- res_list$sigres$twoSVs

allres$sig <- "No"
allres$sig[which(allres$padj < 0.05)] <- "Yes"

allres2 <- DESeq2::lfcShrink(dds, type = "normal", res = results(dds), coef = 3)
allres2 <- as.data.frame(allres2)
allres2$sig <- "No"
allres2$sig[which(allres2$padj < 0.05)] <- "Yes"


ggplot(data = allres2) +
  geom_point(mapping =aes(x = log2FoldChange, y=-log10(padj), color = sig))+
  geom_hline(yintercept = -log10(0.05))


# heatmap

# singular gene plots






