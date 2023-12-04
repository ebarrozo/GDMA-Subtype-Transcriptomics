## GDM-Atlas_v2.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine






############################################ ############################################  ############################################
############################################ Subset Yang et al., dataset that includes atlas metadata ############################################
############################################  ############################################  ############################################

set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)
load("human-rcpa-integrated-umap_v3.RData")

setwd("/home/ebarrozo/gdm.placenta/results/atlas")

dir.create("yang_subset")
setwd("yang_subset")

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
Idents(seurat.object) <- "etal"
levels(seurat.object)

seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap")

Idents(seurat.object) <- "etal"
yang.object <- subset(x = seurat.object, idents = c("Yang-2021"))
ncol(yang.object)  # 37848
rm(seurat.object)

Idents(yang.object) <- "orig.ident"
levels(yang.object)	# "C1" "C2" "C3" "C4" "P1" "P2" "G1" "G2"


Idents(yang.object) <- "orig.ident"
levels(yang.object)
new.metadata <- c("Control", "Control", "Elderly", "Elderly", "PE", "PE", "A1", "A2")
names(new.metadata) <- levels(yang.object)
yang.object <- RenameIdents(yang.object, new.metadata)
yang.object$Condition <- Idents(yang.object)
Idents(yang.object) <- "Condition"


p1 <- DimPlot(yang.object, reduction = "umap", group.by = "orig.ident", label= "TRUE", repel=TRUE, raster=T) + labs(title = NULL, color='Clinical-Group')
p1
p2 <- DimPlot(yang.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=T) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(yang.object, reduction = "umap", group.by = "Condition", label= "TRUE", repel=TRUE, raster=T) + labs(title = NULL, color='Condition')
p3
p5 <- DimPlot(yang.object, reduction = "umap", group.by = "Phase", raster=FALSE) + labs(title = NULL, color='Cell Cycle Phase')
p7 <- FeaturePlot(yang.object, features = 'percent.mt')

p4 <- DimPlot(yang.object, reduction = "umap", group.by = "Type", label= T, repel=TRUE, raster=T) + labs(title = NULL, color='Type')
p4

umap.combined <- p1 + p4 + p3 + p5
ggsave("yang.object.UMAP_integrated.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")


## make these umaps pretty
library(ggpubr)   # install.packages("ggpubr")
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(yang.object, reduction = "umap", group.by = "Platform", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Platform')
p1
p2 <- DimPlot(yang.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(yang.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p3
p4 <- DimPlot(yang.object, reduction = "umap", group.by = "Condition", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p4
p5 <- DimPlot(yang.object, reduction = "umap", group.by = "Phase", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cell Cycle Phase')
p7 <- FeaturePlot(yang.object, features = 'percent.mt') + labs(title = NULL, color='percent.mt')

umap.combined <- p3 + p4 + p5 + p7
ggsave("integrated.UMAP_final_colors.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

Idents(yang.object) <- "Condition"
yang.object <- subset(x = yang.object, idents = c("Control", "A1", "A2"))
ncol(yang.object)  # 19324

## make these umaps pretty
library(ggpubr)   # install.packages("ggpubr")
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(yang.object, reduction = "umap", group.by = "Platform", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Platform')
p1
p2 <- DimPlot(yang.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(yang.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p3
p4 <- DimPlot(yang.object, reduction = "umap", group.by = "Condition", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p4
p5 <- DimPlot(yang.object, reduction = "umap", group.by = "Phase", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cell Cycle Phase')
p7 <- FeaturePlot(yang.object, features = 'percent.mt') + labs(title = NULL, color='percent.mt')

umap.combined <- p3 + p4 + p5 + p7
ggsave("subset.condition.UMAP_final_colors.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")



dir.create("visualization")
setwd("visualization")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

## make these umaps pretty
library(ggpubr)   # install.packages("ggpubr")
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(yang.object, reduction = "umap", group.by = "Platform", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Platform')
p1
p2 <- DimPlot(yang.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(yang.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p3
p4 <- DimPlot(yang.object, reduction = "umap", group.by = "Condition", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p4
p5 <- DimPlot(yang.object, reduction = "umap", group.by = "Phase", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cell Cycle Phase')
p7 <- FeaturePlot(yang.object, features = 'percent.mt') + labs(title = NULL, color='percent.mt')

umap.combined <- p3 + p4 + p5 + p7
ggsave("subset.condition.UMAP_final_colors.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(yang.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(yang.object, features = 'nCount_RNA')
p7 <- FeaturePlot(yang.object, features = 'percent.mt')
p8 <- FeaturePlot(yang.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

setwd("..")

## Scale genes for DE analysis
top3k <- yang.object@assays[["SCT"]]@var.features
yang.object<- SCTransform(yang.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=F, do.scale=T, conserve.memory = F, verbose = TRUE)

####################  Perform  DE analysis by Type 
dir.create("DE_Type")
setwd("DE_Type")
Idents(yang.object) <- "Type"
de_markers <- FindAllMarkers(yang.object, features = intersect(rownames(yang.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_byType_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-Type_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top25_markers-TypebyCondtion_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_TypebyCondtion_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(yang.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top5.markers-Type-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top3.markers-Type-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top25.markers-TypebyCondition_DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

DefaultAssay(yang.object) <- "SCT"
Idents(yang.object) <- "Type"
cluster.averages <- AverageExpression(yang.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Type-top3_mc_heatmap_Type.Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, assay="SCT")
ggsave("DE_Type-top3_FC_heatmap_Type.Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")
####################  Perform  DE analysis by Condition 
dir.create("DE_Condition")
setwd("DE_Condition")
Idents(yang.object) <- "Condition"
de_markers <- FindAllMarkers(yang.object, features = intersect(rownames(yang.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_byCondition_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)

top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = F)
ggsave("integrated_top25_markers-Condition_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(yang.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-ConditionbyType_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_ConditionbyType_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top25.markers-Condition-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top15.markers-Condition-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(yang.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top25.markers-ConditionbyType_DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

DefaultAssay(yang.object) <- "SCT"
Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
cluster.averages <- AverageExpression(yang.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Condition-top25_mc_heatmap_Condition.Condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Condition-top25_FC_heatmap_Condition.Condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
Idents(yang.object) <- "Type"
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Condition-top25_mc_heatmap_Condition.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Condition-top25_FC_heatmap_Condition.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")


setwd("..")


dir.create("recluster.by.top3k")
setwd("recluster.by.top3k")
## run clustering on data with hpv genes
# viral.skin.object <- SCTransform(viral.skin.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=TRUE, do.scale=F, conserve.memory = TRUE, verbose = TRUE)
DefaultAssay(yang.object) <- "SCT"
top3k <- yang.object@assays[["SCT"]]@var.features

yang.object <- ScaleData(yang.object, features = top3k, assay = "SCT", vars.to.regress = c("CC.Difference", "percent.mt", "percent.ribo"))
yang.object<- RunPCA(yang.object, features = top3k, assay = "SCT", slot = "scale.data")
yang.object <- FindNeighbors(yang.object, reduction = "pca", dims = 1:5)
yang.object <- FindClusters(yang.object, res = 0.6, verbose = T)  # res=0.4=  clusters
yang.object <- RunUMAP(yang.object, reduction = "pca", dims = 1:5)
## scale all.genes for logfc vis and DE
yang.object <- ScaleData(yang.object, features = all.genes, assay = "SCT")

dir.create("visualization")
setwd("visualization")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

## make these umaps pretty
library(ggpubr)   # install.packages("ggpubr")
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(yang.object, reduction = "umap", group.by = "Platform", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Platform')
p1
p2 <- DimPlot(yang.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cluster')
p2
p3 <- DimPlot(yang.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=FALSE, cols = mypal3) + labs(title = NULL, color='Type')
p3
p4 <- DimPlot(yang.object, reduction = "umap", group.by = "Condition", label= "F", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Condition')
p4
p5 <- DimPlot(yang.object, reduction = "umap", group.by = "Phase", raster=FALSE, cols = mypal3) + labs(title = NULL, color='Cell Cycle Phase')
p7 <- FeaturePlot(yang.object, features = 'percent.mt') + labs(title = NULL, color='percent.mt')

umap.combined <- p2 + p3 + p4 + p5
ggsave("subset.condition.UMAP_final_colors.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(yang.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(yang.object, features = 'nCount_RNA')
p7 <- FeaturePlot(yang.object, features = 'percent.mt')
p8 <- FeaturePlot(yang.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")

setwd("..")


## Scale genes for DE analysis
# top3k <- yang.object@assays[["SCT"]]@var.features
# yang.object<- SCTransform(yang.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=F, do.scale=T, conserve.memory = F, verbose = TRUE)


####################  Perform  DE analysis by seurat_clusters 
dir.create("DE_seurat_clusters")
setwd("DE_seurat_clusters")
Idents(yang.object) <- "seurat_clusters"
de_markers <- FindAllMarkers(yang.object, features = intersect(rownames(yang.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_byseurat_clusters_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-seurat_clusters_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(yang.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-seurat_clustersbyType_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_seurat_clustersbyType_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, raster = T)
ggsave("integrated_top25_markers-seurat_clustersbyCondition_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(yang.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_seurat_clustersbyCondition_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(yang.object) <- "seurat_clusters"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top5.markers-seurat_clusters-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top3.markers-seurat_clusters-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(yang.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top25.markers-seurat_clustersbyCondition_DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(yang.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(yang.object, features = unique.top2)
png(file=paste0("yang.object_top3.markers-seurat_clustersbyType_DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(yang.object) <- "seurat_clusters"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
DefaultAssay(yang.object) <- "SCT"
Idents(yang.object) <- "seurat_clusters"
cluster.averages <- AverageExpression(yang.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_seurat_clusters-top3_mc_heatmap_seurat_clusters.seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, assay="SCT")
ggsave("DE_seurat_clusters-top3_FC_heatmap_seurat_clusters.seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")

## clean env and save image for monocle on local. 
DefaultAssay(yang.object) <- "SCT"
Idents(yang.object) <- "seurat_clusters"
# yang.object <- DietSeurat(yang.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap", assays= "SCT")

rm(cluster.averages)
rm(de_markers)
rm(df)
rm(df10)
rm(df11)
rm(df12)
rm(df13)
rm(df14)
rm(df15)
rm(df16)
rm(df17)
rm(df2)
rm(df3)
rm(df4)
rm(df5)
rm(df6)
rm(df7)
rm(df8)
rm(df9)
rm(feature.plot)
rm(final)
rm(final2)
rm(library.averages.heatmap)
rm(manifest)
rm(p1)
rm(p2)
rm(p3)
rm(p4)
rm(p5)
rm(p6)
rm(p7)
rm(p8)
rm(top5)
rm(top5.heatmap)
rm(umap.combined)

## clear all of the plots, then reload 
dev.off()


# save.image("yang.object_subset-reclustered.UMAP_v1.RData")

DefaultAssay(yang.object) <- "SCT"
Idents(yang.object) <- "Type"
write.csv(yang.object@meta.data,"yang.object_metadata.csv")
## delete irrelevant columns, add barcode to barcode column, and save as yang.object_metadata_clean.csv

write.csv(yang.object@meta.data,"yang.object_metadata.csv")
# write.csv(yang.object@assays[["RNA"]]@counts,"yang.object_rawcounts.csv") ## 53626 x 19324 :: Takes a very long time, 7.7GB
# write.csv(yang.object@assays[["SCT"]]@counts,"yang.object_sctcounts.csv")

write.table(yang.object@assays[["RNA"]]@counts, file='yang.object_rawcounts.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(yang.object@assays[["SCT"]]@counts, file='yang.object_SCTcounts.tsv', quote=FALSE, sep='\t', col.names = TRUE)

## proceed with ML using raw or SCT counts and associated metadata




############################################ ############################################  ############################################
############################################ Visualize top markers in the spatial datasets ############################################
############################################  ############################################  ############################################

set.seed(seed=1)
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)	# BiocManager::install("glmGamPoi")
library(ggpubr)

setwd("/home/ebarrozo/visium/results/seurat_human_v2")
load("/home/ebarrozo/visium/results/seurat_human_v2/annotated/human_spatial_data-integrated-annotated_v1.RData")
DefaultAssay(seurat.object) <- "Spatial"


DefaultAssay(seurat.object) <- "Spatial"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object) #  "S01"  "S03"  "S04"  "S15"  "S16"  "S17"  "S18"  "S19"  "S20"  "S21"  "S22" "S23a" "S23b" "S24"  "S25"  "S26"
new.metadata <- c("F", "F", "F", "F", 
	"M", "F", "M", "M", 
	"F", "F", "M", "F", 
	"F", "M","M","M")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$FetalSex <- Idents(seurat.object)
Idents(seurat.object) <- "FetalSex"


setwd("/home/ebarrozo/gdm.placenta/results/atlas")

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("NC1a", "NC1b", "NC1c", "ND1", "SP4", "SP1", "SP2", "SP3", "ND2", "ND3", "SP5", "HP1a", "HP1b", "NC2","NC3","NC4")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$SampleCode <- Idents(seurat.object)
Idents(seurat.object) <- "SampleCode"

Idents(seurat.object) <- "orig.ident"
levels(seurat.object)
new.metadata <- c("NegativeControl", "NegativeControl", "NegativeControl", "NotDetected", "SparsePositive", "SparsePositive", "SparsePositive", "SparsePositive", "NotDetected", "NotDetected", "SparsePositive", "HighPositive", "HighPositive", "NegativeControl","NegativeControl","NegativeControl")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$AnalysisCohort <- Idents(seurat.object)
Idents(seurat.object) <- "AnalysisCohort"

# setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3")
# load("human-integrated-umap_v3.RData")

Idents(seurat.object) <- "AnalysisCohort"
levels(seurat.object)
ncol(seurat.object)
	# 17927

seurat.object2 <- subset(x = seurat.object, idents = c("NegativeControl"))
ncol(seurat.object2) # 8481


setwd("/home/ebarrozo/gdm.placenta/results/atlas")
dir.create("TopFeature-SpatialPlots")
setwd("TopFeature-SpatialPlots")


Idents(seurat.object2) <- "orig.ident"
levels(seurat.object2)
new.metadata <- c("Villi", "Decidua", "Membranes",  "Parenchyma","Parenchyma","Parenchyma")
names(new.metadata) <- levels(seurat.object2)
seurat.object2 <- RenameIdents(seurat.object2, new.metadata)
seurat.object2$Region <- Idents(seurat.object2)
Idents(seurat.object2) <- "Region"


GDMA1.feature.list <- c("HES1", "PRSS22", "RMRP")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GDMA1.feature.list, images=image.list, ncol=6)
p5
ggsave("NegativeControl_SpatialPlots-cropped-GDM1-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

GDMA2.feature.list <- c("CD28", "CTSL", "ALB", "CD74")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GDMA2.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-GDM2-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

## path analysis to find top markers?
T2DM.feature.list <- c("CDK17", "CYB5B", "ITGB4", "CSF1R", "ALB", "LINC01002")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = T2DM.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-T2DM-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


GOI.feature.list <- c("HSP90AA1", "HSPA1A", "CD28", "CTSL", "CDK17", "CYB5B", "ITGB4", "CSF1R")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GOI.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-GOI-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

scGDM.feature.list <- c("SERPINE2", "TAC3", "AOC1", "CCL3", "CCL4", "HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = scGDM.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-scGDM-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

# Visualize top markers from volcano plots (figs 2-4)?
volcano.genes <- c("RMRP", "ALB", "PCDH15", "LINC01002")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = volcano.genes, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-volcano-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


# Visualize top markers from volcano plots (figs 2-4)?
FINAL.genes <- c("RMRP", "HES1", "PRSS22", "ALB", "CD28", "CTSL")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = FINAL.genes, images=image.list, ncol=6, crop=T)
ggsave("NegativeControl_SpatialPlots-cropped-FINAL-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 20, units = "in")
p5 <- SpatialPlot(seurat.object2, features = FINAL.genes, images=image.list, ncol=6, crop=F)
ggsave("NegativeControl_SpatialPlots-cropped-FINAL-topmarkers_uncropped.pdf", plot = p5, device = "pdf", width = 10, height = 20, units = "in")



A2.sig <- c("FSTL1", "ECH1", "APOE", "DSTN", "CD151", "HCFC1R1", "BCL11B", "ERGIC3", "PIP4K2A", "ARPC1B", "ENO1", "DCXR", "TPM2", "HAMP", "CLTB", "COL1A2", "MCAM", "PPIB", "CTSL", "GSTA3", "CCL21", "PLAC4", "RELN", "IGKV4-1")
A1.sig <- c("OAZ2", "HES1", "ING5", "CHMP1B", "TCF4", "PPFIBP1", "CTSK", "IGFBP5", "F13A1", "SASH1", "GP1BB", "AKAP12", "MPIG6B", "TMEM88", "LPP", "BID", "EXT2", "UBE2F", "ZNHIT1", "CTSV", "PRSS35", "PLAGL1", "LANCL3", "TPP1", "PAPSS1", "FOXP1", "HIST1H4C", "ABHD16A")
T2DM.sig <- c("TRGC2", "AP1S2", "ADI1", "MOB1B", "CCL14", "UQCRQ", "TNFAIP3", "CITED2", "HHEX", "NAPSA", "NFE2L3", "MRPL27", "F5", "MRPS6", "POLR3GL", "LTBP3", "ANKRD37", "CCDC50", "TCIM", "ITGB4", "ATP6V0B", "NECAP2", "NDUFA4", "MIDN", "TMEM100", "LGALSL", "CLDN5", "ERO1A", "HSD17B12", "LCN2", "IFT57", "NORAD", "CSF1R")

unique.type.genes <- union(A1.sig, A2.sig)
unique.type.genes <- union(unique.type.genes, T2DM.sig)
unique.type.genes <- unique(unique.type.genes)

Idents(seurat.object) <- "Region"
#seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = unique.type.genes, raster = T)
top5.heatmap 
ggsave("Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")

Idents(seurat.object) <- "Region"
#seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = A1.sig, raster = T)
ggsave("A1.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = A2.sig, raster = T)
ggsave("A2.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = T2DM.sig, raster = T)
ggsave("T2DM.sig.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")


Idents(atlas.object) <- "Region"
#atlas.object<- ScaleData(atlas.object, features = all.feature.list, assay = "SCT")
top5.heatmap <- DoHeatmap(atlas.object, features = A1.sig, raster = T)
ggsave("A1.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")
top5.heatmap <- DoHeatmap(atlas.object, features = A2.sig, raster = T)
ggsave("A2.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")
top5.heatmap <- DoHeatmap(atlas.object, features = T2DM.sig, raster = T)
ggsave("T2DM.sig.Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")



## Make heatmaps using gene list
Idents(seurat.object2) <- "Region"

all.feature.list <- union(GDMA1.feature.list, GDMA2.feature.list)
all.feature.list <- union(all.feature.list, T2DM.feature.list)
all.feature.list <- union(all.feature.list, GOI.feature.list)
all.feature.list <- union(all.feature.list, scGDM.feature.list)
all.feature.list <- unique(all.feature.list)
all.feature.list


Idents(seurat.object2) <- "Region"
seurat.object2<- ScaleData(seurat.object2, features = all.feature.list, assay = "Spatial")



############################################ ############################################  ############################################
############################################ Load full atlas for analysis: DE by site, condition, and type ############################################
############################################  ############################################  ############################################


set.seed(seed=1)
setwd("/home/ebarrozo/visium/results/seurat_human_v2/integration-v3/rcpa.integration")
library(dplyr)
library(Seurat)	 # Seurat_4.0.3
library(sctransform)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(glmGamPoi)
library(ggpubr)
load("human-rcpa-integrated-umap_v3.RData")

## Rename site with Region, no annotations
Idents(seurat.object) <- "Site"
levels(seurat.object)	# "PVBP" "DB"   "PV"   "CAM"  "BP"
new.metadata <- c("Parenchyma", "Decidua", "Villi",  "Membranes","BasalPlate")
names(new.metadata) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, new.metadata)
seurat.object$Region <- Idents(seurat.object)
Idents(seurat.object) <- "Region"

setwd("/home/ebarrozo/gdm.placenta/results/atlas")
dir.create("integrated.atlas")
setwd("integrated.atlas")


# seurat.object<- SCTransform(yang.object, method = "glmGamPoi", assay = "SCT", return.only.var.genes=T, do.scale=T, conserve.memory = T, verbose = TRUE)
DefaultAssay(seurat.object) <- "SCT"
seurat.object <- DietSeurat(seurat.object, counts=TRUE, data=TRUE, scale.data = FALSE, dimreducs="umap",assays="SCT")

PrepSCTFindMarkers(seurat.object, assay = "SCT", verbose = TRUE)
	# Given a merged object with multiple SCT models, this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
		# The counts slot of the SCT assay is replaced with recorrected counts and the data slot is replaced with log1p of recorrected counts.
		### Recorrecting SCT counts using minimum median counts: 3965
		# 1632 images present

top3k <- seurat.object@assays[["SCT"]]@var.features



Idents(seurat.object) <- "Region"

#seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
####################  Perform  DE analysis by Region 
dir.create("DE_Region")
setwd("DE_Region")
Idents(seurat.object) <- "Region"
# de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
# SCT assay does not contain median UMI information.Run 
	# `PrepSCTFindMarkers()` before running `FindMarkers()` or invoke `FindMarkers(recorrect_umi=FALSE)`.
de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301, recorrect_umi=FALSE)


write.table(de_markers, "integrated_DEGs_byRegion_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)

A2.sig <- c("FSTL1", "ECH1", "APOE", "DSTN", "CD151", "HCFC1R1", "BCL11B", "ERGIC3", "PIP4K2A", "ARPC1B", "ENO1", "DCXR", "TPM2", "HAMP", "CLTB", "COL1A2", "MCAM", "PPIB", "CTSL", "GSTA3", "CCL21", "PLAC4", "RELN", "IGKV4-1")
A1.sig <- c("OAZ2", "HES1", "ING5", "CHMP1B", "TCF4", "PPFIBP1", "CTSK", "IGFBP5", "F13A1", "SASH1", "GP1BB", "AKAP12", "MPIG6B", "TMEM88", "LPP", "BID", "EXT2", "UBE2F", "ZNHIT1", "CTSV", "PRSS35", "PLAGL1", "LANCL3", "TPP1", "PAPSS1", "FOXP1", "HIST1H4C", "ABHD16A")
T2DM.sig <- c("TRGC2", "AP1S2", "ADI1", "MOB1B", "CCL14", "UQCRQ", "TNFAIP3", "CITED2", "HHEX", "NAPSA", "NFE2L3", "MRPL27", "F5", "MRPS6", "POLR3GL", "LTBP3", "ANKRD37", "CCDC50", "TCIM", "ITGB4", "ATP6V0B", "NECAP2", "NDUFA4", "MIDN", "TMEM100", "LGALSL", "CLDN5", "ERO1A", "HSD17B12", "LCN2", "IFT57", "NORAD", "CSF1R")

unique.type.genes <- union(A1.sig, A2.sig)
unique.type.genes <- union(unique.type.genes, T2DM.sig)
unique.type.genes <- unique(unique.type.genes)

### DO NOT SCALE ALL GENES; ALMOST RAN OUT OF RAM
# seurat.object<- ScaleData(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT")
seurat.object<- ScaleData(seurat.object, features = unique.type.genes, assay = "SCT")


Idents(seurat.object) <- "Region"
#seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = unique.type.genes, raster = T)
ggsave("Unique.Type-Region_heatmap_logfc.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 3, units = "in")
ggsave("Unique.Type-Region_heatmap_logfc2.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 5, units = "in")
ggsave("Unique.Type-Region_heatmap_logfc3.pdf", plot = top5.heatmap, device = "pdf", width = 10, height = 8, units = "in")
ggsave("Unique.Type-Region_heatmap_logfc4.pdf", plot = top5.heatmap, device = "pdf", width = 20, height = 8, units = "in")

DefaultAssay(atlas.object) <- "SCT"
cluster.averages <- AverageExpression(atlas.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.type.genes)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.type.genes, raster = T, slot="counts", assay="SCT") + NoLegend()
ggsave("Unique.Type-Region_MC_MC.png", plot = library.averages.heatmap, device = "png", width = 5, height = 15, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.type.genes, raster = F, assay="SCT") + NoLegend()
ggsave("Unique.Type-Region_MC_logfc.png.png", plot = library.averages.heatmap, device = "png", width = 5, height = 15, units = "in")

library.averages.heatmap <- DoHeatmap(cluster.averages, features = A1.sig, raster = F, assay="SCT") + NoLegend()
ggsave("A1.sig_MC_logfc.png.png", plot = library.averages.heatmap, device = "png", width = 25, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = A2.sig, raster = F, assay="SCT") + NoLegend()
ggsave("A2.sig_MC_logfc.png.png", plot = library.averages.heatmap, device = "png", width = 25, height = 5, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = T2DM.sig, raster = F, assay="SCT") + NoLegend()
ggsave("T2DM.sig_MC_logfc.png.png", plot = library.averages.heatmap, device = "png", width = 25, height = 5, units = "in")





#   No requested features found in the scale.data slot for the SCT assay.
top5.heatmap 
# Zoom in and save a version on desktop as Unique.Type-Region_heatmap_logfc_legend.pdf

seurat.object<- ScaleData(seurat.object, features = top5$gene, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = F)
ggsave("integrated_top25_markers-Region_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-RegionbyType_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_RegionbyType_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(seurat.object) <- "Region"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top25.markers-Region-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top15.markers-Region-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(seurat.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top25.markers-RegionbyType_DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "Region"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Region-top25_mc_heatmap_Region.Region.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Region-top25_FC_heatmap_Region.Region.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
Idents(seurat.object) <- "Type"
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Region-top25_mc_heatmap_Region.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Region-top25_FC_heatmap_Region.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")


setwd("..")


####################  Perform  DE analysis by Type 
dir.create("DE_Type")
setwd("DE_Type")
Idents(seurat.object) <- "Type"
de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301, recorrect_umi=FALSE)
write.table(de_markers, "integrated_DEGs_byType_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

seurat.object<- ScaleData(seurat.object, features = top5$gene, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-Type_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = T)
ggsave("integrated_top25_markers-TypebyCondtion_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_TypebyCondtion_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(seurat.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top5.markers-Type-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top3.markers-Type-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(seurat.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top25.markers-TypebyCondition_DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Type-top3_mc_heatmap_Type.Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, assay="SCT")
ggsave("DE_Type-top3_FC_heatmap_Type.Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

setwd("..")
####################  Perform  DE analysis by Condition 
dir.create("DE_Condition")
setwd("DE_Condition")
Idents(seurat.object) <- "Condition"
de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_byCondition_pos-log2FC.txt", sep="\t")
top5 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
seurat.object<- ScaleData(seurat.object, features = top5$gene, assay = "SCT")

top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = F)
ggsave("integrated_top25_markers-Condition_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top25-markers_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")

Idents(seurat.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, raster = T)
ggsave("integrated_top5_markers-ConditionbyType_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")
top5.heatmap <- DoHeatmap(seurat.object, features = top5$gene, slot="counts", raster = T)
ggsave("integrated_top5-markers_ConditionbyType_counts_pos-log2fc_SCT.pdf", plot = top5.heatmap, device = "pdf", width = 7, height = 5, units = "in")


Idents(seurat.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top25.markers-Condition-DotPlot.png"),res=300, width=2500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top15.markers-Condition-DotPlot.png"),res=300, width=3000, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

Idents(seurat.object) <- "Type"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
feature.plot <- DotPlot(seurat.object, features = unique.top2)
png(file=paste0("seurat.object_top25.markers-ConditionbyType_DotPlot.png"),res=300, width=3500, height=1500)
feature.plot + theme(axis.text.x = element_text(angle = 90)) + scale_x_discrete(name ="Top Markers") + scale_y_discrete(name ="Transcription Profile")
dev.off()

DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "Condition"
top2 <- de_markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
top2 <- top2$gene
unique.top2 <- unique(top2)
seurat.object<- ScaleData(seurat.object, features = unique.top2, assay = "SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=unique.top2)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Condition-top25_mc_heatmap_Condition.Condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Condition-top25_FC_heatmap_Condition.Condition.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
Idents(seurat.object) <- "Type"
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = T, slot="counts", assay="SCT")
ggsave("DE_Condition-top25_mc_heatmap_Condition.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = unique.top2, raster = F, assay="SCT")
ggsave("DE_Condition-top25_FC_heatmap_Condition.Type.png", plot = library.averages.heatmap, device = "png", width = 5, height = 10, units = "in")


setwd("..")












dir.create("DE_Site")
setwd("DE_Site")
Idents(seurat.object) <- "Site"
levels(seurat.object)

# seurat.object<- ScaleData(seurat.object, features = all.genes, assay = "SCT")
# repSCTFindMarkers(seurat.object)
	## ran up to 400/440 gb ram so killing it for now. 
	top3k <- seurat.object@assays[["SCT"]]@var.features

de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_bySite_pos-log2FC.txt", sep="\t")

setwd("..")


Idents(seurat.object) <- "Type"
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=all.feature.list)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = TRUE, slot="counts", assay="SCT")
ggsave("GDM-gene.sets_mc_heatmap_Type-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = TRUE, assay="SCT")
ggsave("GDM-gene.sets_FC_heatmap_Type-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
ggsave("GDM-gene.sets_FC_heatmap_Type-avg_wider.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Type"
heatmap <- DoHeatmap(seurat.object, features = all.feature.list, raster = TRUE, slot="counts", assay="SCT")
ggsave("GDM-gene.sets_mc_heatmap_Type.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
heatmap <- DoHeatmap(seurat.object, features = all.feature.list, raster = TRUE, assay="SCT")
ggsave("GDM-gene.sets_FC_heatmap_Type.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")



Idents(seurat.object) <- "Site"
# seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
DefaultAssay(seurat.object) <- "SCT"

Idents(seurat.object) <- "Site"
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=all.feature.list)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = TRUE, slot="counts", assay="SCT")
ggsave("GDM-gene.sets_mc_heatmap_Site-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = TRUE, assay="SCT")
ggsave("GDM-gene.sets_FC_heatmap_Site-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
ggsave("GDM-gene.sets_FC_heatmap_Site-avg_wider.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Site"
heatmap <- DoHeatmap(seurat.object, features = all.feature.list, raster = TRUE, slot="counts", assay="SCT")
ggsave("GDM-gene.sets_mc_heatmap_Site.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
heatmap <- DoHeatmap(seurat.object, features = all.feature.list, raster = F, assay="SCT")
ggsave("GDM-gene.sets_FC_heatmap_Site.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")

library(ggpubr)
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

p1 <- DimPlot(seurat.object, reduction = "umap", group.by = "etal", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='et al.')
p1
ggsave("UMAP_XL-etal.pdf", plot = p1, device = "pdf", width = 8.5, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "seurat_clusters", label= "TRUE", repel=TRUE, raster=TRUE, cols = mypal3) + labs(title = NULL, color='Clusters')
p2
ggsave("UMAP-Clusters.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p2 <- DimPlot(seurat.object, reduction = "umap", group.by = "Type", label= "TRUE", repel=TRUE, raster=TRUE, cols = mypal3) + labs(title = NULL, color='Profile')
p2
ggsave("UMAP-Type.pdf", plot = p2, device = "pdf", width = 10, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Site", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Site')
ggsave("UMAP-Site.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p3 <- DimPlot(seurat.object, reduction = "umap", group.by = "Platform", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Platform')
ggsave("UMAP-Platform.pdf", plot = p3, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "Condition", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='Condition')
ggsave("UMAP-Condition.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

umap.combined <- p1 + p2 + p3 + p4
ggsave("integrated_UMAP-FINAL.pdf", plot = umap.combined, device = "pdf", width = 17, height = 12, units = "in")
ggsave("UMAP_XL.pdf", plot = umap.combined, device = "pdf", width = 22, height = 12, units = "in")



#################### 
#################### 
Idents(seurat.object) <- "Region"
feature.list <- c("HSP90AA1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSP90AA1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "HSP90AA1", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_HSP90AA1_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "Type"
feature.list <- c("HSP90AA1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSP90AA1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "HSP90AA1", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_HSP90AA1_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


Idents(seurat.object) <- "Region"
feature.list <- c("HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSPA1A.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "HSPA1A", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_HSPA1A_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "Type"
feature.list <- c("HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSPA1A.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "HSPA1A", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_HSPA1A_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


Idents(seurat.object) <- "Region"
feature.list <- c("CSH1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-CSH1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "CSH1", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_CSH1_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object) <- "Type"
feature.list <- c("CSH1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-CSH1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object, features = "CSH1", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_CSH1_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


