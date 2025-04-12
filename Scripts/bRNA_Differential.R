## sample script for tutorial using Mov10 dataset
## adapted from HBC tutorial https://hbctraining.github.io/DGE_workshop/schedule/1.5-day.html

library(ggplot2)
library(ggrepel)
library(DESeq2)
library(magrittr)
library(pheatmap)
library(fdrtool)
library(dplyr)
library(tidyr)
library(sva)
library(tibble)
#### 1. Loading data and  first look ####
data <- read.table(paste0("./OmicsBasicsTheory/data/Mov10_full_counts.txt"), header = T, row.names = 1)
meta <- read.table(paste0("./OmicsBasicsTheory/data/Mov10_full_meta.txt"), header = T, row.names = 1)

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

## classic running and first look into results
resultsNames(dds)
res <- results(dds, alpha = 0.05, contrast = c("sampletype", "MOV10_overexpression", "control"))
summary(res)
# you can also get the results in a table format
res_tbl <- as.data.frame(results(dds, contrast = c("sampletype", "MOV10_overexpression", "control")))


## explaining pvalue histograms
hist(res$pvalue, breaks=200)

## fdrtool
FDR.ddsRes <- fdrtool(res$stat, statistic= "normal", plot = T)
res_afterFDRtool <- res
res_afterFDRtool$pvalue <- FDR.ddsRes$pval
res_afterFDRtool[,"padj"]  <- p.adjust(res_afterFDRtool$pvalue, method = "BH")

# comparing before and after fdrtool
hist(res$pvalue, breaks=200)
hist(res_afterFDRtool$pvalue, breaks=200)

summary(res)
summary(res_afterFDRtool)

# how many of after fdrtool significant genes are in the significant gene list from befor fdrtool?
table(rownames(res_afterFDRtool)[which(res_afterFDRtool$padj<0.05)] %in% rownames(res)[which(res$padj<0.05)])
# what is the ranking of the still significant genes in the before fdrtool gene list?
plot(which(rownames(res_afterFDRtool)[which(res_afterFDRtool$padj<0.05)] %in% rownames(res)[which(res$padj<0.05)]))

## log2Fold shrinkage
# shrinkage also has the types apeglm and ashr. Check Deseq2 documentation
res_afterFDRtool_lfcShrinked <- lfcShrink(dds, res = res_afterFDRtool, type = "normal", coef = "sampletype_MOV10_overexpression_vs_control")
# after deseq2 and shrinkage
plotMA(res_afterFDRtool)
plotMA(res_afterFDRtool_lfcShrinked)

#### 6. Interpretation of output and visualization ####
res_afterFDRtool_tbl <- res_afterFDRtool %>% data.frame() %>% rownames_to_column(var="gene")
res_afterFDRtool_lfcShrinked_tbl <- res_afterFDRtool_lfcShrinked %>% data.frame() %>% rownames_to_column(var="gene")
# volcano plot
res_afterFDRtool_tbl$sig <- "No"
res_afterFDRtool_tbl$sig[which(res_afterFDRtool_tbl$padj < 0.05)] <- "Yes"

res_afterFDRtool_lfcShrinked_tbl$sig <- "No"
res_afterFDRtool_lfcShrinked_tbl$sig[which(res_afterFDRtool_lfcShrinked_tbl$padj < 0.05)] <- "Yes"

# before lfc shrinkage
ggplot(data = res_afterFDRtool_tbl) +
  geom_point(mapping =aes(x = log2FoldChange, y=-log10(padj), color = sig))+
  geom_hline(yintercept = -log10(0.05))

# after lfc shrinkage
ggplot(data = res_afterFDRtool_lfcShrinked_tbl) +
  geom_point(mapping =aes(x = log2FoldChange, y=-log10(padj), color = sig))+
  geom_hline(yintercept = -log10(0.05))

## heatmap: good for checking if effect is replicated in replicates
normalized_counts <- counts(dds, normalized = TRUE)
sig_norm <- data.frame(normalized_counts) %>% rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% res_afterFDRtool_lfcShrinked_tbl$gene[which(res_afterFDRtool_lfcShrinked_tbl$padj < 0.05)])
# filtering for control and Mov10 overexpressed
sig_norm <- sig_norm[,-c(2,3)]
pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
         cluster_rows = T, 
         show_rownames = F,
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)    

# singular gene plots
# with default plotCounts
plotCounts(dds, gene="MOV10", intgroup="sampletype", returnData = FALSE, normalized = TRUE)

# with more advanced labelling options in ggplot2
d <- plotCounts(dds, gene="MOV10", intgroup="sampletype", returnData = TRUE, normalized = TRUE)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d)), max.overlaps = 15) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))





