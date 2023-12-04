# GO_v1.R 

## inspired by 
# http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
# https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_2021_June/Markdowns/12_Gene_set_testing.pdf
# https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html


# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/GDM-SigPooled")
load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/GDM-SigPooled/gdm.placenta_data-final.RData")


# https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwiCm7CIsuGAAxVNlGoFHcKEApkQFnoECA4QAQ&url=https%3A%2F%2Fedu.sib.swiss%2Fpluginfile.php%2F8237%2Fmod_folder%2Fcontent%2F0%2FRNASEQ20_Day3_HandsOn.pdf%3Fforcedownload%3D1&usg=AOvVaw3hBKrMwrUzvprrh1peZCZu&opi=89978449
library(DESeq2)
library(pheatmap)
# library(org.Mm.eg.db) # BiocManager::install("org.Mm.eg.db")
library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
library(DOSE) # BiocManager::install("DOSE")
library(pathview) # BiocManager::install("pathview")
library(clusterProfiler) # BiocManager::install("clusterProfiler")
library(AnnotationHub) # BiocManager::install("AnnotationHub")
library(ensembldb) # BiocManager::install("ensembldb")
library(tidyverse)

setwd("..")
dir.create("functional-enrichment_v2")
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2")

# deseq_counts <- DESeqDataSetFromMatrix(countData = df5.1, 
 #                                      colData =dd_meta,
 #                                      design = ~Diabetes) 
# dds <-DESeq(deseq_counts,parallel = T)
# res <- results(dds, contrast=c("Diabetes","Yes","No"), alpha=0.05, cooksCutoff=FALSE)

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = df5.1, 
                                      colData =dd_meta,
                                      design = ~DiabetesSubtype)
view(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# Transform normalized counts using the rlog function
# To improve the distances/clustering for the PCA and heirarchical clustering visualization methods, we need to moderate the variance across the mean by applying the rlog transformation to the normalized counts.
rld <- rlog(dds, blind=TRUE) # This may take a long time

#Principal components analysis (PCA)
# DESeq2 has a built-in function for plotting PCA plots, that uses ggplot2 under the hood. This is great because it saves us having to type out lines of code and having to fiddle with the different ggplot2 layers. In addition, it takes the rlog object as an input directly, hence saving us the trouble of extracting the relevant information from it.
# The function plotPCA() requires two arguments as input: an rlog object and the intgroup (the column in our metadata that we are interested in).
plotPCA(rld, intgroup="DiabetesSubtype")

## save it
png(file=paste0("PCAplot_DiabetesSubtype_500.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="DiabetesSubtype")
dev.off()
png(file=paste0("PCAplot_SampleID_500.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="SampleID")
dev.off()
png(file=paste0("PCAplot_DiabetesYesNo_500.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="Diabetes")
dev.off()
# What does this plot tell you about the similarity of samples? Does it fit the expectation from the experimental design? 
  # By default the function uses the top 500 most variable genes. You can change this by adding the ntop argument and specifying how many genes you want to use to draw the plot.

plotPCA(rld, intgroup="DiabetesSubtype", ntop=3000)



png(file=paste0("PCAplot_DiabetesSubtype_3k.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="DiabetesSubtype", ntop=3000)
dev.off()
png(file=paste0("PCAplot_SampleID_3k.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="SampleID", ntop=3000)
dev.off()
png(file=paste0("PCAplot_DiabetesYesNo_3k.png"),res=300, width=2500, height=1500)
plotPCA(rld, intgroup="Diabetes", ntop=3000)
dev.off()


plotPCA(rld, intgroup="DiabetesSubtype", ntop=3000)

pcaData <- plotPCA(rld, intgroup="DiabetesSubtype", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p1<-ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype), ellipse = T) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  stat_ellipse()
p1

png(file=paste0("PCAplot_DiabetesSubtype_3k_ellipse.png"),res=300, width=2500, height=1500)
p1
dev.off()

### Extract the rlog matrix from the object
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)
# Overall, we observe pretty high correlations across the board ( > 0.999) suggesting no outlying sample(s). 
# Also, similar to the PCA plot you see the samples clustering together by sample group. 
# Together, these plots suggest to us that the data are of good quality and we have the green light to proceed to differential expression analysis.

png(file=paste0("Heatmap_SampleID_pairwisecorrelation-default.png"),res=300, width=2500, height=1500)
pheatmap(rld_cor)
dev.off()

pheatmap(mat = rld_cor,
cutree_rows = 6,cutree_cols = 5,border_color = "black")

## Save with the code
png(file=paste0("Heatmap_SampleID_pairwisecorrelation-cutmodules.png"), res=300, width=3000, height=2000)
pheatmap(mat = rld_cor, cutree_rows = 6,cutree_cols = 5,border_color = "black")
dev.off()

## annotate by DiabetesSubtype and sample
# anot_col<-data.frame(df)%>%select(DiabetesSubtype)
# rownames(anot_col)<-df$SampleID

# pheatmap(mat = rld_cor, annotation_col = anot_col, cutree_rows = 6,cutree_cols = 5,border_color = "black")

## Save with the code
#png(file=paste0("Heatmap_DiabetesSubtype-SampleID_pairwisecorrelation-cutmodules.png"),res=300, width=3000, height=2000)
#pheatmap(mat = rld_cor, annotation_col = anot_col,cutree_rows = 3,cutree_cols = 7,border_color = "black")
#dev.off()



# dds <- DESeq(dds)
dds <-DESeq(dds,parallel = T)


plotDispEsts(dds)

##################################################################################################################################################################################
########################################################### Control vs A1 #####################################
##################################################################################################################################################################################

dir.create("A1")
setwd("A1")

# res <- results(dds)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","A1"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05

res <- results(dds, contrast=c("DiabetesSubtype","Control","A1"))
summary(res)

plotMA(res, ylim=c(-2,2)) # The genes that are significantly DE are colored to be easily identified.

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("all")
setwd("all")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
A1.all.sig.genes <- row.names(signif_res)
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
A1.up.sig.genes <- row.names(signif_res)
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
A1.down.sig.genes <- row.names(signif_res)

signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("histone H3-K4 methylation", "methylation",
                    "chromatin remodeling", "epigenetic regulation of gene expression")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[1:3000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

dir.create("up")
setwd("up")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("histone modification", "epigenetic regulation of gene expression",
                    "myeloid cell differentiation", "regulation of mRNA metabolic process")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
signif_res <- signif_res[1:3000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("down")
setwd("down")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
	# 2359 genes
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
	############################################ Manually select categories from a node of interest ############################################
categorys <- c("mRNA processing", "tRNA modification",
                    "methylation", "protein targeting")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
signif_res <- signif_res[1:2000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually select categories from a node of interest ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")
setwd("..")

##################################################################################################################################################################################
########################################################### Control vs A2 #####################################
##################################################################################################################################################################################

dir.create("A2")
setwd("A2")

# res <- results(dds)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","A2"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05

res <- results(dds, contrast=c("DiabetesSubtype","Control","A2"))
summary(res)

plotMA(res, ylim=c(-2,2)) # The genes that are significantly DE are colored to be easily identified.

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("all")
setwd("all")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
A2.all.sig.genes <- row.names(signif_res)
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
A2.up.sig.genes <- row.names(signif_res)
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
A2.down.sig.genes <- row.names(signif_res)

signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("RNA splicing", "aerobic respiration",
                    "protein folding", "establishment of protein localization to organelle")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[1:3000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

dir.create("up")
setwd("up")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
A2.up.sig.genes <- signif_res
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("regulation of actin filiment based process", "cell-substrate adhesion",
                    "antigen processing and presentation of peptide antigen", "ribonucleoprotein complex assembly")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
signif_res <- signif_res[1:1500,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("down")
setwd("down")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
	# 2359 genes
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
	############################################ Manually select categories from a node of interest ############################################
categorys <- c("regulation of protein ubiquitination", "ncRNA processing",
                    "mitochondrial respiratory chain complex assembly", "mitochondrion transport")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
A2.down.sig.genes <- signif_res 
signif_res <- signif_res[1:1500,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually select categories from a node of interest ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")
setwd("..")

##################################################################################################################################################################################
########################################################### Control vs T2DM #####################################
##################################################################################################################################################################################

dir.create("T2DM")
setwd("T2DM")

# res <- results(dds)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","T2DM"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05

res <- results(dds, contrast=c("DiabetesSubtype","Control","T2DM"))
summary(res)

plotMA(res, ylim=c(-2,2)) # The genes that are significantly DE are colored to be easily identified.

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("all")
setwd("all")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
T2DM.all.sig.genes <- row.names(signif_res)
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
T2DM.up.sig.genes <- row.names(signif_res)
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
T2DM.down.sig.genes <- row.names(signif_res)

signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("chromatin remodeling", "methylation",
                    "RNA splicing", "nuclear transport")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[1:2500,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

dir.create("up")
setwd("up")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
		############################################ Manually select categories from a node of interest ############################################

categorys <- c("histone modification", "negative regulation of gene expression, epigenetic",
                    "myeloid cell differentiation", "mRNA processing")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
signif_res <- signif_res[1:3000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually selected categories from a node of interest above ############################################
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")

############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler
dir.create("down")
setwd("down")
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
	# 2359 genes
signif_genes <- as.character(rownames(signif_res))

## I have common gene names not ENSEMBL ID
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)


# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.table(cluster_summary, file="GOEnrichmentAnalysis_results.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)
emapplot(x2, cex_category=1.5)
	## based on this, select node(s) of interest 
	############################################ Manually select categories from a node of interest ############################################
categorys <- c("mRNA processing", "cellular respiration",
                    "ncRNA processing", "mitochondrial translation")


png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 6714 genes, too many to plot


## subset the top  genes for visualization
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange < -0.301 & !is.na(signif_res$log2FoldChange), ]
signif_res <- signif_res[1:3000,]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "SYMBOL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

png(file=paste0("cnetplot-top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego,categorySize="pvalue",showCategory = 5,foldChange= signif_res_lFC, vertex.label.font=6)
dev.off()

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

png(file=paste0("cnetplot-top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 
dev.off()
############################################ Manually select categories from a node of interest ############################################

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories ################## Rename plot if changing node
png(file=paste0("cnetplot-cat_top3k_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-cat_top3k_default.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()

treeplot(x2, hclust_method = "average")

png(file=paste0("treeplot_q.05_avg.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
png(file=paste0("treeplot_default.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()

# while (!is.null(dev.list()))  dev.off()
setwd("..")
setwd("..")

####################################################################################################################################
###################################### make a venn diagram and table of the up and down genes for each group #######################################################
####################################################################################################################################
library(VennDiagram)



venn_data <- list(
  A1.up = A1.up.sig.genes,
    A1.down = A1.down.sig.genes,
    A2.up = A2.up.sig.genes,
    A2.down = A2.down.sig.genes
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1.up", "A1.down", "A2.up", "A2.down"),
  filename = NULL
)
grid.draw(venn.plot)

png(file=paste0("A1.A2_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

venn_data <- list(
    A2.up = A2.up.sig.genes,
    A2.down = A2.down.sig.genes,
    T2DM.up = T2DM.up.sig.genes,
    T2DM.down = T2DM.down.sig.genes
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A2.up", "A2.down", "T2DM.up", "T2DM.down"),
  filename = NULL)

grid.draw(venn.plot)


png(file=paste0("A2.T2DM_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

venn_data <- list(
  A1.up = A1.up.sig.genes,
    A1.down = A1.down.sig.genes,
    T2DM.up = T2DM.up.sig.genes,
    T2DM.down = T2DM.down.sig.genes
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1.up", "A1.down","T2DM.up", "T2DM.down"),
  filename = NULL
)

grid.draw(venn.plot)

png(file=paste0("A1.T2DM_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()


## 
# Create a list of gene lists
gene_lists <- list(
    A1.up = A1.up.sig.genes,
    A1.down = A1.down.sig.genes,
    A2.up = A2.up.sig.genes,
    A2.down = A2.down.sig.genes,
    T2DM.up = T2DM.up.sig.genes,
    T2DM.down = T2DM.down.sig.genes
)

# Perform intersections and set differences
results <- list()
for (i in 1:(length(gene_lists) - 1)) {
  for (j in (i + 1):length(gene_lists)) {
    list1_name <- names(gene_lists)[i]
    list2_name <- names(gene_lists)[j]
    list1 <- gene_lists[[i]]
    list2 <- gene_lists[[j]]
    
    intersect_genes <- intersect(list1, list2)
    setdiff_list1 <- setdiff(list1, intersect_genes)
    setdiff_list2 <- setdiff(list2, intersect_genes)
    
    results[[paste(list1_name, "_and_", list2_name, "_Intersection", sep = "")]] <- intersect_genes
    results[[paste(list1_name, "_minus_", list2_name, "_SetDiff", sep = "")]] <- setdiff_list1
    results[[paste(list2_name, "_minus_", list1_name, "_SetDiff", sep = "")]] <- setdiff_list2
  }
}

# Write the results to individual CSV files
for (result_name in names(results)) {
  result_genes <- results[[result_name]]
  result_df <- data.frame(Genes = result_genes)
  write.csv(result_df, paste(result_name, ".csv", sep = ""), row.names = FALSE)
}

# Perform intersections and set differences and make a combined results table
results <- list()
for (i in 1:(length(gene_lists) - 1)) {
  for (j in (i + 1):length(gene_lists)) {
    list1_name <- names(gene_lists)[i]
    list2_name <- names(gene_lists)[j]
    list1 <- gene_lists[[i]]
    list2 <- gene_lists[[j]]
    
    intersect_genes <- intersect(list1, list2)
    setdiff_list1 <- setdiff(list1, intersect_genes)
    setdiff_list2 <- setdiff(list2, intersect_genes)
    
    results[[paste(list1_name, "_and_", list2_name, "_Intersection", sep = "")]] <- intersect_genes
    results[[paste(list1_name, "_minus_", list2_name, "_SetDiff", sep = "")]] <- setdiff_list1
    results[[paste(list2_name, "_minus_", list1_name, "_SetDiff", sep = "")]] <- setdiff_list2
  }
}


# Combine the results into a single data frame
combined_results <- data.frame(Genes = unique(unlist(results)))

# Add columns for each pair of gene lists
for (result_name in names(results)) {
  combined_results[[result_name]] <- combined_results$Genes %in% results[[result_name]]
}

# Write the combined_results data frame to a CSV file
write.csv(combined_results, "combined_venn_UP.DOWN_results.csv", row.names = FALSE)




# Find the union of all genes
all_genes <- unique(unlist(gene_lists))

# Create a data frame with the results
results_df <- data.frame(Genes = all_genes)

# Populate the data frame with presence indicators
for (list_name in names(gene_lists)) {
  results_df[[list_name]] <- as.integer(results_df$Genes %in% gene_lists[[list_name]])
}

# Save the results as a .txt table - This one is binary, so could be used for machine learning? 
write.table(results_df, "results_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

####################################################################################################################################
####################################### Make venn diagrams of overlapping significant up/down GO pathways #######################################################
####################################################################################################################################

dir.create("GO-Venns")
setwd("GO-Venns")

A1.up.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/A1/up/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
A1.up <- A1.up.GO$Description
A1.down.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/A1/down/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
A1.down <- A1.down.GO$Description
A2.up.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/A2/up/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
A2.up <- A2.up.GO$Description
A2.down.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/A2/down/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
A2.down <- A2.down.GO$Description
T2DM.up.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/T2DM/up/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
T2DM.up <- T2DM.up.GO$Description
T2DM.down.GO <-read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.1/functional-enrichment_v2/T2DM/down/GOEnrichmentAnalysis_results.txt",sep = "\t",header = T)
T2DM.down <- T2DM.down.GO$Description

library(VennDiagram)



venn_data <- list(
  A1.up = A1.up,
    A1.down = A1.down,
    A2.up = A2.up,
    A2.down = A2.down
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1.up", "A1.down", "A2.up", "A2.down"),
  filename = NULL
)
grid.draw(venn.plot)

png(file=paste0("A1.A2_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

venn_data <- list(
    A2.up = A2.up,
    A2.down = A2.down,
    T2DM.up = T2DM.up,
    T2DM.down = T2DM.down
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A2.up", "A2.down", "T2DM.up", "T2DM.down"),
  filename = NULL)

grid.draw(venn.plot)


png(file=paste0("A2.T2DM_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

venn_data <- list(
  A1.up = A1.up,
    A1.down = A1.down,
    T2DM.up = T2DM.up,
    T2DM.down = T2DM.down
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1.up", "A1.down","T2DM.up", "T2DM.down"),
  filename = NULL
)

grid.draw(venn.plot)

png(file=paste0("A1.T2DM_Up.Down_Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

# Create a list of gene lists
gene_lists <- list(
    A1.up = A1.up,
    A1.down = A1.down,
    A2.up = A2.up,
    A2.down = A2.down,
    T2DM.up = T2DM.up,
    T2DM.down = T2DM.down
)

# Perform intersections and set differences
results <- list()
for (i in 1:(length(gene_lists) - 1)) {
  for (j in (i + 1):length(gene_lists)) {
    list1_name <- names(gene_lists)[i]
    list2_name <- names(gene_lists)[j]
    list1 <- gene_lists[[i]]
    list2 <- gene_lists[[j]]
    
    intersect_genes <- intersect(list1, list2)
    setdiff_list1 <- setdiff(list1, intersect_genes)
    setdiff_list2 <- setdiff(list2, intersect_genes)
    
    results[[paste(list1_name, "_and_", list2_name, "_Intersection", sep = "")]] <- intersect_genes
    results[[paste(list1_name, "_minus_", list2_name, "_SetDiff", sep = "")]] <- setdiff_list1
    results[[paste(list2_name, "_minus_", list1_name, "_SetDiff", sep = "")]] <- setdiff_list2
  }
}

# Write the results to individual CSV files
for (result_name in names(results)) {
  result_genes <- results[[result_name]]
  result_df <- data.frame(Genes = result_genes)
  write.csv(result_df, paste(result_name, ".csv", sep = ""), row.names = FALSE)
}


# Perform intersections and set differences and make a combined results table
results <- list()
for (i in 1:(length(gene_lists) - 1)) {
  for (j in (i + 1):length(gene_lists)) {
    list1_name <- names(gene_lists)[i]
    list2_name <- names(gene_lists)[j]
    list1 <- gene_lists[[i]]
    list2 <- gene_lists[[j]]
    
    intersect_genes <- intersect(list1, list2)
    setdiff_list1 <- setdiff(list1, intersect_genes)
    setdiff_list2 <- setdiff(list2, intersect_genes)
    
    results[[paste(list1_name, "_and_", list2_name, "_Intersection", sep = "")]] <- intersect_genes
    results[[paste(list1_name, "_minus_", list2_name, "_SetDiff", sep = "")]] <- setdiff_list1
    results[[paste(list2_name, "_minus_", list1_name, "_SetDiff", sep = "")]] <- setdiff_list2
  }
}


# Combine the results into a single data frame
combined_results <- data.frame(Genes = unique(unlist(results)))

# Add columns for each pair of gene lists
for (result_name in names(results)) {
  combined_results[[result_name]] <- combined_results$Genes %in% results[[result_name]]
}

# Write the combined_results data frame to a CSV file
write.csv(combined_results, "combined_venn_UP.DOWN_results.csv", row.names = FALSE)




# Find the union of all genes
all_genes <- unique(unlist(gene_lists))

# Create a data frame with the results
results_df <- data.frame(Genes = all_genes)

# Populate the data frame with presence indicators
for (list_name in names(gene_lists)) {
  results_df[[list_name]] <- as.integer(results_df$Genes %in% gene_lists[[list_name]])
}

# Save the results as a .txt table - This one is binary, so could be used for machine learning? 
write.table(results_df, "results_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

setwd("..")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

## come back and make category visualization of pathways unique or conserved among signatures










