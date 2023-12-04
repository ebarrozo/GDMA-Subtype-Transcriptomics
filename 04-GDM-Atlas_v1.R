## GDM-Atlas_v1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## Analyze data on AagaardLab3
# 


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

feature.list <- c("HOXA10")

image.list <- c("S01", "S03.1", "S04.2", 
    "S24.13", "S25.14", "S26.15",
 "S15.3", "S20.8", "S21.9", 
 "S16.4", "S17.5", "S18.6", 
 "S19.7", "S22.10", "S23a.11", "S23b.12")
p5 <- SpatialPlot(atlas.object, features = feature.list, images=image.list, ncol=6)
ggsave("HOXA10-cropped-GDM1-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 20, units = "in")




GDMA1.feature.list <- c("HPSE", "KRT81", "SAA2")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GDMA1.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-GDM1-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

GDMA2.feature.list <- c("TUBA3E", "CCL14", "TAC3", "HOXD10", "CXCL9", "MT1H")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GDMA2.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-GDM2-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

T2DM.feature.list <- c("OMD", "PDZK1IP1", "SERPINA3", "NDP", "GSTA1", "AADAC", "FGB", "RXFP1", "SCARA5", "RBP4", "IGFBP1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = T2DM.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-T2DM-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")


GOI.feature.list <- c("CSH1", "EGFR", "PER1", "SOGA1", "IL2RB", "FOXO1", "SOD3", "DOCK5", "PIK3CB", "DIO2", "IGFBP1", "FAT2", "TAC3", "CSH2", "EIF6", "FN1", "IGF2", "SLC2A1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = GOI.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-GOI-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")

scGDM.feature.list <- c("HSP90AA1", "HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = scGDM.feature.list, images=image.list, ncol=6)
ggsave("NegativeControl_SpatialPlots-cropped-scGDM-topmarkers.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")



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


dir.create("DE_Site")
setwd("DE_Site")
Idents(seurat.object) <- "Site"

seurat.object<- ScaleData(seurat.object, features = all.genes, assay = "SCT")
PrepSCTFindMarkers(seurat.object)
de_markers <- FindAllMarkers(seurat.object, features = intersect(rownames(seurat.object), all.genes), assay = "SCT", only.pos = TRUE, min.pct = 0.20, logfc.threshold =  0.301)
write.table(de_markers, "integrated_DEGs_bySite_pos-log2FC.txt", sep="\t")


Idents(seurat.object2) <- "Region"
cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.genes)
write.csv(cluster.averages@assays[["SCT"]]@counts, file = "pseudobulk.viralmeancounts_byType.csv")

####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("GDM-gene.sets_mc_heatmap_region-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = FALSE, assay="Spatial")
ggsave("GDM-gene.sets_FC_heatmap_region-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
ggsave("GDM-gene.sets_FC_heatmap_region-avg_wider.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


heatmap <- DoHeatmap(seurat.object2, features = all.feature.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("GDM-gene.sets_mc_heatmap_Region.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
heatmap <- DoHeatmap(seurat.object2, features = all.feature.list, raster = FALSE, assay="Spatial")
ggsave("GDM-gene.sets_FC_heatmap_Region.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object2) <- "Type"
cluster.averages <- AverageExpression(seurat.object2, assays="Spatial", slot="counts", return.seurat = TRUE, features=all.feature.list)
####Use this to determine the sum of mean counts for t.lineages transcripts per cluster and rank t.lineages transcripts based on overall t.lineages transcription (especially for determining pseudotimes start point) 
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("GDM-gene.sets_mc_heatmap_Type-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
library.averages.heatmap <- DoHeatmap(cluster.averages, features = all.feature.list, raster = FALSE, assay="Spatial")
ggsave("GDM-gene.sets_FC_heatmap_Type-avg.png", plot = library.averages.heatmap, device = "png", width = 5, height = 8, units = "in")
ggsave("GDM-gene.sets_FC_heatmap_Type-avg_wider.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")


Idents(seurat.object2) <- "Type"
heatmap <- DoHeatmap(seurat.object2, features = all.feature.list, raster = FALSE, slot="counts", assay="Spatial")
ggsave("GDM-gene.sets_mc_heatmap_Type.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")
#heatmap of fold-change for t.lineages transcripts for each cluster
heatmap <- DoHeatmap(seurat.object2, features = all.feature.list, raster = FALSE, assay="Spatial")
ggsave("GDM-gene.sets_FC_heatmap_Type.png", plot = heatmap, device = "png", width = 11, height = 8, units = "in")


#################### 
#################### 
Idents(seurat.object2) <- "Region"
feature.list <- c("HSP90AA1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSP90AA1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "HSP90AA1", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_HSP90AA1_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object2) <- "Type"
feature.list <- c("HSP90AA1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSP90AA1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "HSP90AA1", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_HSP90AA1_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")


Idents(seurat.object2) <- "Region"
feature.list <- c("HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSPA1A.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "HSPA1A", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_HSPA1A_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object2) <- "Type"
feature.list <- c("HSPA1A")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-HSPA1A.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "HSPA1A", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_HSPA1A_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")



Idents(seurat.object2) <- "Region"
feature.list <- c("CSH1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-CSH1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "CSH1", assay="Spatial", group.by="Region")
library.averages.heatmap
ggsave("vlnplot_CSH1_Region.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")
Idents(seurat.object2) <- "Type"
feature.list <- c("CSH1")
image.list <- c("S01", "S03.1", "S04.2", "S24.13", "S25.14", "S26.15")
p5 <- SpatialPlot(seurat.object2, features = feature.list, images=image.list, ncol=3)
ggsave("NegativeControl_SpatialPlots-cropped-CSH1.pdf", plot = p5, device = "pdf", width = 10, height = 8, units = "in")
library.averages.heatmap <- VlnPlot(seurat.object2, features = "CSH1", assay="Spatial", group.by="Type")
library.averages.heatmap
ggsave("vlnplot_CSH1_Type.pdf", plot = library.averages.heatmap, device = "pdf", width = 10, height = 5, units = "in")

############# corr.plot in all and in subsets? 


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


Idents(seurat.object) <- "Type"
seurat.object<- ScaleData(seurat.object, features = all.feature.list, assay = "SCT")
DefaultAssay(seurat.object) <- "SCT"

setwd("/home/ebarrozo/gdm.placenta/results/atlas")
dir.create("TopFeature-IntegratedPlots_atlas")
setwd("TopFeature-IntegratedPlots_atlas")

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

p4 <- DimPlot(seurat.object, reduction = "umap", group.by = "FetalSex", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='FetalSex')
ggsave("UMAP-FetalSex.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")
p4 <- FeaturePlot(seurat.object, features="DDX3Y", reduction = "umap", cols = mypal3, raster=TRUE, label= "TRUE", repel=TRUE) + labs(title = NULL, color='DDX3Y')
ggsave("UMAP-DDX3Y.pdf", plot = p4, device = "pdf", width = 8.5, height = 6, units = "in")

#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_RNA')
p6 <- FeaturePlot(seurat.object, features = 'nCount_RNA')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-RNA.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- FeaturePlot(seurat.object, features = 'nFeature_Spatial'
p6 <- FeaturePlot(seurat.object, features = 'nCount_Spatial')
p7 <- FeaturePlot(seurat.object, features = 'percent.mt')
p8 <- FeaturePlot(seurat.object, features = 'percent.ribo')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots-Spatial.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")
#examine UMAPs with qc metrics 
p5 <- DimPlot(seurat.object, reduction = "umap", group.by = "Phase", cols = mypal3, raster=TRUE) + labs(title = NULL, color='Cell Cycle Phase')
p6 <- FeaturePlot(seurat.object, features = 'percent.Hemo')
p7 <- FeaturePlot(seurat.object, features = 'percent.Spike')
p8 <- FeaturePlot(seurat.object, features = 'percent.viral')
umap.combined <- CombinePlots(plots = list(p5, p6, p7, p8))
ggsave("integrated_UMAP_QCmetricsFeaturePlots2.pdf", plot = umap.combined, device = "pdf", width = 20, height = 12, units = "in")


Idents(seurat.object) <- "Platform"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Platform", label=TRUE)
ggsave("integrated_UMAP_splitby_platform-platform.pdf", plot = p2, device = "pdf", width = 9, height = 3, units = "in")
Idents(seurat.object) <- "Site"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Site", label=TRUE)
ggsave("integrated_UMAP_splitby_Site-Site.pdf", plot = p2, device = "pdf", width = 12, height = 3, units = "in")
Idents(seurat.object) <- "Condition"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Condition", label=TRUE)
ggsave("integrated_UMAP_splitby_Condition-Condition.pdf", plot = p2, device = "pdf", width = 12, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "etal", label=TRUE)
ggsave("integrated_UMAP_splitby_etal-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Condition", label=TRUE)
ggsave("integrated_UMAP_splitby_Condition-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")
Idents(seurat.object) <- "Type"
p2 <- DimPlot(seurat.object, reduction = "umap", split.by = "Site", label=TRUE)
ggsave("integrated_UMAP_splitby_Site-Type.pdf", plot = p2, device = "pdf", width = 20, height = 3, units = "in")


###### UMAP + Spatial DimPlots
DefaultAssay(seurat.object) <- "SCT"
plot2 <- SpatialFeaturePlot(seurat.object, features="SIGLEC1", ncol = 1) + patchwork::plot_layout(ncol = 1)
####### FOR SOME REASON, INTEGRATION ADDS IMAGES ARBITRARILY
ggsave("integrated_SpatialDimPlot_SIGLEC1_TEST.pdf", plot = plot2, device = "pdf", width = 8, height = 45, units = "in")
	## determine the correct imageID by copying and pasting them here

DefaultAssay(seurat.object) <- "SCT"
Idents(seurat.object) <- "etal"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_etal.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "Type"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_Type.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

Idents(seurat.object) <- "seurat_clusters"
seurat.object <- ScaleData(seurat.object, features = SARS.genes.receptors.hb.fetalsex.mac, assay="SCT")
cluster.averages <- AverageExpression(seurat.object, assays="SCT", slot="counts", return.seurat = TRUE, features=SARS.genes.receptors.hb.fetalsex.mac)
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, slot="counts", assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_mc_heatmap_seurat_clusters.png", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")
library.averages.heatmap <- DoHeatmap(cluster.averages, features = SARS.genes.receptors.hb.fetalsex.mac, raster = FALSE, assay="SCT")
ggsave("SARS.genes.receptors.hb.fetalsex.mac_FC_heatmap_seurat_clusters", plot = library.averages.heatmap, device = "png", width = 11, height = 8, units = "in")

