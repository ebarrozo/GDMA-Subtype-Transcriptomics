## gdm.placenta_bulk-RNA-seq_script_v1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## Analyze data on AagaardLab3
# 

# bulk RNA-seq from human placenta tissue: Control, GDMA1, and GDMA2, T2D
## metadata located in 'RNA seq_Clinical Validation_13 May 2016_v2 (2).xlsx' /Users/enricobarrozo/Library/CloudStorage/Box-Box/Melissa-GDM/RNA seq_Clinical Validation_13 May 2016_v2 (2).xlsx
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta
ls
    ## DM1 thru DM11, PE, stored as .bz2 e.g. DM1.read1.bz2 and DM1.read2.bz2 ; 22 total
    ## NCS5 thru 88, PE, stored as .bz2 e.g. NCS5.read1.bz2 and NCS5.read2.bz2 ; 54 total ; maybe ChIP? based on NCS designations in metadata excel sheet
    ## 6002 thru 6171, PE, stored as .bz2 e.g. 6002.read1.bz2 and 6002.read2.bz2 ; 51 total ; 6005 only has read1

################################################################################
##########################    Data Wrangling in Terminal      ######################################################
################################################################################
## in Terminal make a copy of the data in case there is an error
cd /home/ebarrozo/gdm.placenta/results
rm -r DESeq2_analysis_v1
mkdir DESeq2_analysis_v1
cp *.txt /home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1
cd /home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1

## fix filenames
rm NCS51abam_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
mv NCS60.paired_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt NCS60_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS60.paired_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS60.repair_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS62.repair_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt

## rm files that don't have metadata condition
rm X6030_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm X6151_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm X6153_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm DM5_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm DM9_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt

## add an X into the front of the filenames for 6002_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt samples
        mv 6002_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6002_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6003_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6003_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6022_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6022_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6025_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6025_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6026_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6026_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6033_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6033_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6035_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6035_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6039_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6039_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6040_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6040_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6041_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6041_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6064_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6064_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6070_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6070_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6073_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6073_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6076_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6076_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6077_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6077_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6090_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6090_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6105_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6105_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6157_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6157_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6162_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6162_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6170_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6170_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
        mv 6171_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt X6171_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
## done manually


## Run Sample_regex.sh in terminal 
  ## removes last 5 lines that are the counts summary
  ## Removes column A
for n in $(ls *.txt);
do
  echo "running $n";
  cat $n|cut -f2-3>$n.tmp.txt;
  head -n-5 $n.tmp.txt>$n;
  rm $n.tmp.txt;
done
  ## now there are 2 columns: genes and counts

## make a list of the file names we can call on later
ls *.txt > list.tsv

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 0_meta.R
## Go back to RStudio 
## Set the working directory
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1")

## load the libraries
library(tidyverse)
library(mosaic)
library(pheatmap)
library(vegan) # install.packages("vegan")
library(ggpubr)
library(ggsci)

## load the data using the list
df<-as_tibble(read.table("list.tsv",header = F,sep = "\t"))

# View the data
df  ## 53 samples
# write.table(x = df,file = "samples.tsv",sep = "\t",row.names = T)

## RNA-seq filenames    JMR30_HFD_1.trimmed_PE.fq.gz.duprm.sorted.htseq.counts.txt
## ChIP filenames     
## gdm.placenta         NCS62_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt

# Clip file names
df<-df%>%mutate(name=gsub("_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt","",V1))%>%select(-V1)
df  # 6002 

## Need to do something about the NAs in the replicate numbers
tally(~Sample,df)

data.frame(df)

df
## save the metadata
# write.table(df,"meta.tsv",sep = "\t",row.names = F)

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 1_import.R
#import the metadata df that we made in the previous script
getwd()
#meta<-as_tibble(read.table("meta.tsv",sep = "\t",header = T))

## from metadata using sampleID and class. Added an X in front of 6000 samples for R
meta<-as_tibble(read.table("/home/ebarrozo/gdm.placenta/docs/ERB_Metadata_slim_vX.tsv",sep = "\t",header = T))
meta
# Meta has SampleID and Class followed by X, X.1, thru X.21
    ## X thru X.1 are empty
## Remove those empty columns
keeps <- c("SampleID","Class")
meta <- meta[keeps]
## Need to do something about the NAs in the replicate numbers
tally(~SampleID,meta)
data.frame(meta)
## save the metadata
write.table(meta,"meta.tsv",sep = "\t",row.names = F)
meta<-as_tibble(read.table("meta.tsv",sep = "\t",header = T))
meta

#make a dataframe with a column that matches the names of all the count files in the results folder
df<-data.frame(file=list.files("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1"))
df<-df%>%filter(file!="list.tsv")
df<-df%>%filter(file!="meta.tsv")
df<-df%>%filter(file!="all-PCA")

df<-df%>%separate(col = file,into = c("Sample"),sep = "_",remove = F)
## View df
df

dim(df) ## 53 rows x 2

#merge the dfs together
# df <-full_join(df,meta)
df <-right_join(df,meta, by=c("Sample"="SampleID"))

## Confirm merge worked ; rows must match dims above
dim(df) #  [1] 53 x 3


#start up a dataframe that can be used as a matrix for joining
tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1/",df$file[1]),sep = "\t",header = F))

colnames(tmp)<-c("gene",df$Sample[1])
# colnames(tmp)<-c("gene",df$Sample[1])
tmp
df2<-tmp

unique(df$Sample)

#make a function that iteratively opens the file and merges the gene name a counts into a single matrix
jmac<-function(X)
{
  #import the file and annotate a column name
    tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1/",df$file[X]),sep = "\t",header = F))
  #rename the column names to match for merging
  colnames(tmp)<-c("gene",df$Sample[X])
  #merge the gene and counts dataframe with the larger file
  return(tmp)
}

#run the function across all the files
df3<-as_tibble(merge(x = df2,y = lapply(X = rep(2:length(df$Sample)),FUN = function(X) jmac(X))))
#remove the redundant gene columns
df4<-df3%>%select(gene,!contains("gene"))
#convert it to a matrix with the rownames being the gene
df5<-as.matrix(df4%>%select(-gene))
rownames(df5)<-df4$gene


dir.create("all-PCA")
setwd("all-PCA")
#write the count matrix to a table
 write.table(x = df5,file = "count_matrix.tsv",sep = "\t",row.names = T)
  ## YAY we have a counts matrix :)
################################################################################
################################################################################
##########################    DESeq2 Analysis and PCA Plot  ######################################################
################################################################################ EB modified from MJ 2_analysis.R


#convert the counts to a distance matrix
d<-vegdist(x = t(df5),method = "euclidean")
# Error in vegdist(x = t(df), method = "euclidean") : 
 # input data must be numeric

# Draw a plot for a non-vegan ordination (cmdscale).
dis <- vegdist(t(df5),method = "gower")
mds <- cmdscale(dis, eig = TRUE)

mds<-data.frame(mds$points)
mds
mds<-mds%>%mutate(Sample=rownames(mds))
mds
df$Sample
df<-full_join(df,mds)
df<-as_tibble(df)%>%
  mutate(Eig1=X1,Eig2=X2)%>%
  select(-c(X1,X2))%>%
  mutate(across(where(is.character),as.factor))


## mutate metadata table
dd_meta<-meta%>%
  # filter(Sample!="JMR16"|Tx!="HFRes2"|rep!=2)%>%  ### no need to remove samples from the metadata. They were removed upstream in this version. 
  #filter(Sample!="JMR30"|Tx!="HFD"|rep!=1)%>%
  mutate(across(where(is.character),factor))
dd_meta




################################################################################
ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Class",
          shape = "Class",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Class",
          shape = "Class",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")
a

## save PCA plot 
png(file=paste0("PCA_plot_default-Class.png"),
                res=300, 
                width=2000, 
                height=1500)
a
dev.off()

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Sample",
          shape = "Sample",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-Sample.png"),
                res=300, 
                width=2500, 
                height=2000)
a
dev.off()



# Add vertical and horizontal line to a ggscatterhist
plots <- ggscatterhist(df, x = "Eig2", y = "Eig1",
                       margin.params = list(theme= theme_pubr()),
                       color="Class",
                       palette = "aaas",print = FALSE)
plots$sp<-plots$sp+
  stat_density_2d(aes(colour = Class,size=0.05),contour = T,alpha=0.2,size=0.05,
                  contour_var = "count",n = 25)+
  xlim(c(-0.25,0.25))+
  ylim(-0.3,0.3)

plots$xplot<-plots$xplot+theme_pubr()
plots$yplot<-plots$yplot+theme_pubr()
plots

## save custom PCA plot 
png(file=paste0("PCA_custom-plot.png"),
                res=300, 
                width=2000, 
                height=1500)
plots
# ggpar(p = plots,legend.title = "Group") ## this command will change the legend from Class to whatever you want
dev.off()


### Make a custom plot 2
tally(~Class,df)
data.frame(df)
# customise plot
customised_plot <- 
  ggplot(df,mapping = aes(x = Eig2,y = Eig1, z=0, color=Class)) +
  geom_point()+
  geom_contour() +
scale_colour_brewer(palette = "Set1") +
  theme_pubr(legend = "bottom")+xlim(c(-0.25,0.25))+ylim(-0.3,0.3)+
  coord_fixed(ratio = 0.5, clip = "off")

 customised_plot

png(file=paste0("PCA_plot_simple.png"),
                res=300, 
                width=2000, 
                height=1500)
customised_plot
dev.off()

## End of DataWranglin_v2.R script, begin gdm.placenta_bulk-RNA-seq_script_v1.1R script
################################################################################
##########################    DESeq2 Load Object  ######################################################
################################################################################ EB modified from MJ help 7/20/2022
#Set the randomnness to a particular seed for reproducibility downstream
set.seed(seed=1)
# Load the required packages
library(DESeq2) # BiocManager::install("DESeq2")
library(dplyr)

## mutate metadata table
dd_meta<-meta%>%
  # filter(Sample!="JMR16"|Tx!="HFRes2"|rep!=2)%>%  ### no need to remove samples from the metadata. They were removed upstream in this version. 
  #filter(Sample!="JMR30"|Tx!="HFD"|rep!=1)%>%
  mutate(across(where(is.character),factor))
dd_meta

#add up all the duplicate transcripts
df4.1<-df4%>%
  group_by(gene)%>%
  summarise(across(everything(),sum))
df4.2<-df4.1 # %>%select(-JMR30_HFD) ## Sample removed upstream

df4.2
#convert 
df5.1<-as.matrix(df4.2%>%select(-gene))
df5.1
colSums(df5.1)

## Add gene as a column name
rownames(df5.1)<-df4.1$gene

## Finally, load data and metadata into DESeq2
deseq_counts <- DESeqDataSetFromMatrix(countData = df5.1, 
                                       colData =dd_meta,
                                       design = ~Class) 
getwd()
dds <-DESeq(deseq_counts,parallel = T)
## Make an elbowplot to determine how many PCs in a PCA = > 85% of variance. This will determine cutreen=
## https://support.bioconductor.org/p/83626/
     library(matrixStats)
     
     #How to get PCA plot?
     ##how to obtain d.deseq was described in DESeq2 manual
     
     cds=estimateDispersions(dds)
     vsd=varianceStabilizingTransformation(cds)
     plotPCA(vsd,intgroup=c("Class"))
     p1 <-      plotPCA(vsd,intgroup=c("Class"))
     p1
     #How to get PCA scree plot?
     
     ## calculate the variance for each gene
     rv <- rowVars(assay(vsd))

     ## select the ntop genes by variance
     select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
     
     ## perform a PCA on the data in assay(x) for the selected genes
     pca <- prcomp(t(assay(vsd)[select,]))
     
     ## the contribution to the total variance for each component
     percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
     
     ##plot the "percentVar"
     scree_plot=data.frame(percentVar)
     scree_plot[,2]<- c(1:24)
     p2 <-      scree_plot[,2]<- c(1:24)


     colnames(scree_plot)<-c("variance","component_number")
     ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")
     p3 <-      ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")

png(file=paste0("PCAplot-manual.png"),
                res=300, 
                width=2500, 
                height=1500)
p1
dev.off()
png(file=paste0("PCAplot-scree_plot-default.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

png(file=paste0("PCAplot-scree_plot-custom.png"),
                res=300, 
                width=2500, 
                height=1500)
p3
dev.off()

## elbow plot results -> use 4-5 PCs for >85% of variance
# save.image("gdm.placenta_data-v4.RData")

setwd("..")
###### DESeq2 will have 3 iterations ctrl vs GDMA1 or GDMA2 or T2DM
#########################################################################################################
##########################    Iteration I:    Control-vs-GDMA1Class Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1")
load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1/all-PCA/gdm.placenta_data-v4.RData")

dir.create("Control-vs-GDMA1Class")
setwd("Control-vs-GDMA1Class")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Class","Control","GDMA1"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the GDMA1 
  ## Wald hypothesis testing results in 0 significant genes. Coinsider LTR 
write.csv(res, file = "ControlClass_GDMA1_DESeq2_allResults-Wald.csv")


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
res <- results(dds, contrast=c("Class","Control","GDMA1"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
write.csv(res, file = "ControlClass_GDMA1_DESeq2_allResults-LTR.csv")
  ## LTR results in 13 up and 17 down sig genes; 30 total



########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$pvalue),]
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
    # 30 
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
  ## 1969
head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]
write.csv(topResults, file = "ControlClass_GDMA1_DESeq2_topResults.csv")
write.csv(resOrdered, file = "ControlClass_GDMA1_DESeq2_allResults-LTR-ordered.csv")

########## Generate a heatmap of the counts with only the Control and GDMA1 Class samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 30 sig genes
library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
res3
res4<-full_join(res3,df%>%select(Sample,Class)%>%distinct_all())
res4
res5<-res4%>%filter(Class%in%c("Control","GDMA1"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
sig2<-full_join(sig,df%>%select(Sample,Class)%>%distinct_all())
sig3<-sig2%>%filter(Class%in%c("Control","GDMA1"))%>%distinct_all()%>%
  group_by(Sample,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes)
top.genes.Control.GDMA1 <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = Sample,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 30 37
df_slim<-df%>%
  select(Sample,Class)%>%
  distinct_all()%>%
  filter(Class%in%c("Control","GDMA1"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(Class)
rownames(anot_col)<-df_slim$Sample

anot_col
library(RColorBrewer)

####### Default heatmap
pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")

p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")
png(file=paste0("Control-VS-GDMA1_counts-heatmap_row-scaled_notcut.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()


p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 3,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()
## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut3.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()

####### Heatmap cutting clusters
p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

## Save with the code
png(file=paste0("counts-heatmap_row-scaled-cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()
############# End Analysis with Michael 7/20/22

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("Class"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Class, shape=Class)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

png(file=paste0("PCAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

############ Plot histogram with 1 gene. Can be modified to compare any gene of interest.
png(file=paste0("ControlClass_GDMA1Class_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='Class')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlClass_GDMA1Class_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlClass_GDMA1Class_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlClass_GDMA1Class_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlClass_GDMA1Class_p-hist_smooth.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue[resOrdered$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()


# Make a basic volcano plot
#BiocManager::install("genefilter")
library("genefilter")
library(gplots)
topVarGenes <- head(order(-rowVars(assay(vst))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vst$GDMA1 ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$GDMA1,"-",vst$Class)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlClass_GDMA1Class_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
with(res, plot(log2FoldChange, -log10(pvalue))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


# Enhanced Volcano plot
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

png(file=paste0("ControlClass_GDMA1Class_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
# png(file=paste0("Top20-GDMA1Class_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
# pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) ## breaks my plots in RStudio
dev.off()


# png(file=paste0("Top20-GDMA1Class_vst-heatmap_rowsclust.png"),
                res=300, 
                width=2500, 
                height=1500)
# # pheatmap(assay(vst)[select,], cluster_rows=TRUE, show_rownames=TRUE)
dev.off()
## End

setwd("..")
#########################################################################################################
##########################    Iteration II:    Control-vs-GDMA2Class Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1")
load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1/all-PCA/gdm.placenta_data-v4.RData")

dir.create("Control-vs-GDMA2Class")
setwd("Control-vs-GDMA2Class")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Class","Control","GDMA2"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the GDMA2 
  ## Wald hypothesis testing results in 0 significant genes. Coinsider LTR 
write.csv(res, file = "ControlClass_GDMA2_DESeq2_allResults-Wald.csv")


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
res <- results(dds, contrast=c("Class","Control","GDMA2"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
write.csv(res, file = "ControlClass_GDMA2_DESeq2_allResults-LTR.csv")
  ## LTR results in 13 up and 17 down sig genes; 30 total



########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$pvalue),]
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
    # 30 
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
  ## 1969
head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]
write.csv(topResults, file = "ControlClass_GDMA2_DESeq2_topResults.csv")
write.csv(resOrdered, file = "ControlClass_GDMA2_DESeq2_allResults-LTR-ordered.csv")

########## Generate a heatmap of the counts with only the Control and GDMA2 Class samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 30 sig genes
library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
res3
res4<-full_join(res3,df%>%select(Sample,Class)%>%distinct_all())
res4
res5<-res4%>%filter(Class%in%c("Control","GDMA2"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
sig2<-full_join(sig,df%>%select(Sample,Class)%>%distinct_all())
sig3<-sig2%>%filter(Class%in%c("Control","GDMA2"))%>%distinct_all()%>%
  group_by(Sample,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes)
top.genes.Control.GDMA2 <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = Sample,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 30 37
df_slim<-df%>%
  select(Sample,Class)%>%
  distinct_all()%>%
  filter(Class%in%c("Control","GDMA2"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(Class)
rownames(anot_col)<-df_slim$Sample

anot_col
library(RColorBrewer)

####### Default heatmap
pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")

p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")
png(file=paste0("Control-VS-GDMA2_counts-heatmap_row-scaled_notcut.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()


p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 3,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()
## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut3.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()

####### Heatmap cutting clusters
p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

## Save with the code
png(file=paste0("counts-heatmap_row-scaled-cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()
############# End Analysis with Michael 7/20/22

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("Class"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Class, shape=Class)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

png(file=paste0("PCAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

############ Plot histogram with 1 gene. Can be modified to compare any gene of interest.
png(file=paste0("ControlClass_GDMA2Class_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='Class')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlClass_GDMA2Class_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlClass_GDMA2Class_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlClass_GDMA2Class_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlClass_GDMA2Class_p-hist_smooth.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue[resOrdered$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()


# Make a basic volcano plot
#BiocManager::install("genefilter")
library("genefilter")
library(gplots)
topVarGenes <- head(order(-rowVars(assay(vst))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vst$GDMA2 ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$GDMA2,"-",vst$Class)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlClass_GDMA2Class_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
with(res, plot(log2FoldChange, -log10(pvalue))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


# Enhanced Volcano plot
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

png(file=paste0("ControlClass_GDMA2Class_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
# png(file=paste0("Top20-GDMA2Class_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
# pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) ## breaks my plots in RStudio
dev.off()


# png(file=paste0("Top20-GDMA2Class_vst-heatmap_rowsclust.png"),
                res=300, 
                width=2500, 
                height=1500)
# # pheatmap(assay(vst)[select,], cluster_rows=TRUE, show_rownames=TRUE)
dev.off()
## End

setwd("..")
#########################################################################################################
##########################    Iteration III:    Control-vs-T2DMClass Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1")
load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v1/all-PCA/gdm.placenta_data-v4.RData")

dir.create("Control-vs-T2DMClass")
setwd("Control-vs-T2DMClass")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Class","Control","T2DM"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the T2DM 
  ## Wald hypothesis testing results in 0 significant genes. Coinsider LTR 
write.csv(res, file = "ControlClass_T2DM_DESeq2_allResults-Wald.csv")


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
res <- results(dds, contrast=c("Class","Control","T2DM"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
write.csv(res, file = "ControlClass_T2DM_DESeq2_allResults-LTR.csv")
  ## LTR results in 13 up and 17 down sig genes; 30 total



########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$pvalue),]
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
    # 30 
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
  ## 1969
head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]
write.csv(topResults, file = "ControlClass_T2DM_DESeq2_topResults.csv")
write.csv(resOrdered, file = "ControlClass_T2DM_DESeq2_allResults-LTR-ordered.csv")

########## Generate a heatmap of the counts with only the Control and T2DM Class samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 30 sig genes
library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
res3
res4<-full_join(res3,df%>%select(Sample,Class)%>%distinct_all())
res4
res5<-res4%>%filter(Class%in%c("Control","T2DM"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
sig2<-full_join(sig,df%>%select(Sample,Class)%>%distinct_all())
sig3<-sig2%>%filter(Class%in%c("Control","T2DM"))%>%distinct_all()%>%
  group_by(Sample,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes)
top.genes.Control.T2DM <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = Sample,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 30 37
df_slim<-df%>%
  select(Sample,Class)%>%
  distinct_all()%>%
  filter(Class%in%c("Control","T2DM"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(Class)
rownames(anot_col)<-df_slim$Sample

anot_col
library(RColorBrewer)

####### Default heatmap
pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")

p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")
png(file=paste0("Control-VS-T2DM_counts-heatmap_row-scaled_notcut.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()


p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 3,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()
## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut3.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()

####### Heatmap cutting clusters
p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

## Save with the code
png(file=paste0("counts-heatmap_row-scaled-cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()
############# End Analysis with Michael 7/20/22

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("Class"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Class, shape=Class)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

png(file=paste0("PCAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

############ Plot histogram with 1 gene. Can be modified to compare any gene of interest.
png(file=paste0("ControlClass_T2DMClass_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='Class')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlClass_T2DMClass_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlClass_T2DMClass_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlClass_T2DMClass_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlClass_T2DMClass_p-hist_smooth.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue[resOrdered$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()


# Make a basic volcano plot
#BiocManager::install("genefilter")
library("genefilter")
library(gplots)
topVarGenes <- head(order(-rowVars(assay(vst))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vst$T2DM ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$T2DM,"-",vst$Class)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlClass_T2DMClass_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
with(res, plot(log2FoldChange, -log10(pvalue))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()


# Enhanced Volcano plot
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

png(file=paste0("ControlClass_T2DMClass_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
# png(file=paste0("Top20-T2DMClass_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
# pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) ## breaks my plots in RStudio
dev.off()


# png(file=paste0("Top20-T2DMClass_vst-heatmap_rowsclust.png"),
                res=300, 
                width=2500, 
                height=1500)
# # pheatmap(assay(vst)[select,], cluster_rows=TRUE, show_rownames=TRUE)
dev.off()
## End

setwd("..")

#########################################################################################################
##########################    Iteration 4:    GDM-SigPooled Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 

## This iteration 8 requires iterations 1-3 above. 

dir.create("GDM-SigPooled")
setwd("GDM-SigPooled")

## Combine the list of genes significantly DE for GDM relative to all controls
GDM.genes <- union(top.genes.Control.GDMA1, top.genes.Control.GDMA2)
GDM.genes <- unique(GDM.genes)
GDM.genes
write.csv(GDM.genes, file = "merged-sig.GDM-genes.csv")
write.csv(top.genes.Control.GDMA1, file = "Control.GDMA1-sig-genes.csv")
write.csv(top.genes.Control.GDMA2, file = "Control.GDMA2-sig-genes.csv")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
res <- results(dds, alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## LTR results in 17 up and 13 down sig genes


########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA
resultsNames(de) # [1] "Intercept"                       "Class_GDM_vs_Control"         "Class_GDMRescue_vs_Control"  
        # [4] "Class_GDMA1_vs_Control"

resOrdered <- res[order(res$padj),]
res05 <- res
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]

########## Generate a heatmap of with all samples and pooled GDM significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
    ## 30 sig genes

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
res3
res4<-full_join(res3,df%>%select(Sample,Class)%>%distinct_all())
res4
res5<-res4%>%filter(Class%in%c("Control","GDMA1", "GDMA2", "T2DM"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "Sample",values_to = "count")
sig2<-full_join(sig,df%>%select(Sample,Class)%>%distinct_all())
sig3<-sig2%>%filter(Class%in%c("Control","GDMA1", "GDMA2", "T2DM"))%>%distinct_all()%>%
  group_by(Sample,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  



sig4<-sig3%>%
  pivot_wider(names_from = Sample,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))


h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-GDM.genes
dim(h)      # [1] 30 53
df_slim<-df%>%
  select(Sample,Class)%>%
  #distinct_all()%>%
  #filter(Class%in%c("Control","GDMA1", "GDMA2", "T2DM"))%>%
  distinct_all()
df_slim     # 53 × 2
df
anot_col<-data.frame(df_slim)%>%select(Class)
rownames(anot_col)<-df_slim$Sample

anot_col
library(RColorBrewer)

####### Default heatmap
pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")

p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 1,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1
## Save with the code
png(file=paste0("counts-heatmap_row-scaled_notcut.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()
p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 4,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_colcut4.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()

####### Heatmap cutting clusters
p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p1

## Save with the code
png(file=paste0("counts-heatmap_row-scaled-cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 5,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut5by5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 3,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut3by5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()



p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and colums
         annotation_col = anot_col,cutree_rows = 3,cutree_cols = 3,border_color = "black", ## Can change cutree_rows and colums
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cutcols3_cutrows3.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()


# Enhanced Volcano plot
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res2,
    lab = rownames(res2),
    x = 'log2FoldChange',
    y = 'pvalue')
p3
png(file=paste0("Enhancedvolcanoplot-sig.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

png(file=paste0("Enhancedvolcanoplot-all.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

p1 <- EnhancedVolcano(res,
  lab = rownames(res),
  x = 'log2FoldChange',
  y = 'pvalue',
  pCutoff = 10e-4,
  FCcutoff = 1.333,
  xlim = c(-5.5, 5.5),
  ylim = c(0, -log10(10e-12)),
  pointSize = 1.5,
  labSize = 4.5,
  title = 'DESeq2 results',
  subtitle = 'Differential expression',
  caption = 'FC cutoff, 1.333; p-value cutoff, 10e-4',
  legendPosition = "right",
  legendLabSize = 14,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.9,
  drawConnectors = TRUE,
  hline = c(10e-8),
  widthConnectors = 0.5)
p1

png(file=paste0("Enhancedvolcanoplot-default.png"),
                res=300, 
                width=3500, 
                height=2500)
p1
dev.off()


# Make a basic volcano plot
#BiocManager::install("genefilter")
library("genefilter")
library(gplots)
vst <- vst(dds, blind=FALSE)
topVarGenes <- head(order(-rowVars(assay(vst))),35)
topVarGenes
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vst$Control ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$Control,"-",vst$Class)

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)    # 31
#reset par
par(mfrow=c(1,1))

png(file=paste0("Basic_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
with(res, plot(log2FoldChange, -log10(pvalue)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

# save.image("gdm.placenta_data-final.RData")



####### 
sessionInfo()


