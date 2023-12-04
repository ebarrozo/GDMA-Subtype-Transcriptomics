## gdm.placenta_bulk-RNA-seq_script_v3.2.1.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## v3 - utilizes new diabetes subtype categories including unknowns
  ## v3.1 changes the dash in the metatable names to make sure DE analysis was correct. 
  ## v3.2 renames A1 and A2 to GDMA1 and GDMA2, and excludes the in between samples NCS27, NCS51, NCS88

## metadata located in 'DEIDENTIFIED_TABLE1.xlsx' /Users/enricobarrozo/Library/CloudStorage/Box-Box/AagaardLab/ManuscriptDrafts/GDM/AJOG/DEIDENTIFIED_TABLE1.xlsx

################################################################################
##########################    Data Wrangling in Terminal      ######################################################
################################################################################
## in Terminal make a copy of the data in case there is an error
cd /home/ebarrozo/gdm.placenta/results
rm -r DESeq2_analysis_v3.2
mkdir DESeq2_analysis_v3.2
cp *.txt /home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2
cd /home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2

## fix filenames
mv NCS51abam_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt NCS51_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
mv NCS60.paired_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt NCS60_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS60.paired_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS60.repair_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS62.repair_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt

## exclude DM7 because it is a type I diabetic
rm DM7_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
# 8/14/23- Adding back these samples into metadata because we have the metadata and sequencing data.
#cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/6030_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO
#cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/6151_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO
#cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/6153_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO
#cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/DM5_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO
#cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/DM9_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO
# cp /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta/DM7_1.trimmed_PE.fastq.gz /home/ebarrozo/gdm.placenta/data/GEO

## add an X into the front of the filenames for 6002_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt samples

rm NCS27_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS51_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
rm NCS88_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt


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
  ## 54 samples;; Matches Table 1

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 0_meta.R
## Go back to RStudio 
## Set the working directory
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")

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
df  ## 57 samples
# write.table(x = df,file = "samples.tsv",sep = "\t",row.names = T)

## RNA-seq filenames    JMR30_HFD_1.trimmed_PE.fq.gz.duprm.sorted.htseq.counts.txt
## ChIP filenames     
## gdm.placenta         NCS62_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt

# Clip file names
df<-df%>%mutate(name=gsub("_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt","",V1))%>%select(-V1)
df  # X6002 

## Need to do something about the NAs in the replicate numbers
tally(~name,df)

data.frame(df)

df
## save the metadata
write.table(df,"meta.tsv",sep = "\t",row.names = F)

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 1_import.R
#import the metadata df that we made in the previous script
getwd()
#meta<-as_tibble(read.table("meta.tsv",sep = "\t",header = T))


## update meta with deidentified sampleID, DiabetesSubtype, and fetalbiologicalsex
## from metadata using sampleID and DiabetesSubtype. Added an X in front of 6000 samples for R
meta<-as_tibble(read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/meta.tsv",sep = "\t",header = T))
meta
# Meta has SampleID and DiabetesSubtype followed by X, X.1, thru X.21
    ## X thru X.1 are empty
## Remove those empty columns
keeps <- c("SampleID","DiabetesSubtype", "FetusBiologicalSex1F2M")
meta <- meta[keeps]
## Need to do something about the NAs in the replicate numbers
tally(~SampleID,meta)
data.frame(meta)
## save the metadata
write.table(meta,"meta_final.tsv",sep = "\t",row.names = F)
\
#make a dataframe with a column that matches the names of all the count files in the results folder
df<-data.frame(file=list.files("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2"))
df<-df%>%filter(file!="list.tsv")
df<-df%>%filter(file!="meta_final.tsv")
df<-df%>%filter(file!="meta.tsv")
## exclude all-PCA if re-running this
df<-df%>%filter(file!="all-PCA")

df<-df%>%separate(col = file,into = c("SampleID"),sep = "_",remove = F)
## View df
df

dim(df) ## 54 rows x 2

#merge the dfs together
# df <-full_join(df,meta)
#df <-right_join(df,meta, by=c("Sample"="SampleID"))
df2 <-full_join(df,meta, by=c("SampleID"))
df<- df2
rm(df2)

## Confirm merge worked ; rows must match dims above
dim(df) #  [1]54 x 5


#start up a dataframe that can be used as a matrix for joining
tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/",df$file[1]),sep = "\t",header = F))

colnames(tmp)<-c("gene",df$SampleID[1])
# colnames(tmp)<-c("gene",df$Sample[1])
tmp
df2<-tmp

unique(df$SampleID)

#make a function that iteratively opens the file and merges the gene name a counts into a single matrix
jmac<-function(X)
{
  #import the file and annotate a column name
    tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/",df$file[X]),sep = "\t",header = F))
  #rename the column names to match for merging
  colnames(tmp)<-c("gene",df$SampleID[X])
  #merge the gene and counts dataframe with the larger file
  return(tmp)
}

#run the function across all the files
df3<-as_tibble(merge(x = df2,y = lapply(X = rep(2:length(df$SampleID)),FUN = function(X) jmac(X))))
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
  ## in terminal - copy over to the GEO folder
#  cp /home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/all-PCA/count_matrix.tsv /home/ebarrozo/gdm.placenta/data/GEO
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
mds<-mds%>%mutate(SampleID=rownames(mds))
mds
df$SampleID
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
mypal1 <-get_palette("ucscgb",50)
mypal2 <-get_palette("aaas",50)
mypal3 <-get_palette("igv",50)
mypal4 <-get_palette("npg",50)

ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
         #  shape = "DiabetesSubtype",
          ellipse = T,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")




## save PCA plot 
png(file=paste0("PCA_plot_default-DiabetesSubtype.png"),
                res=300, 
                width=2000, 
                height=1500)
ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
         #  shape = "DiabetesSubtype",
          ellipse = T,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = F,
          title = "PCA plot of sample distances by gene count")
dev.off()

png(file=paste0("PCA_plot_default-DiabetesSubtype.png"),
                res=300, 
                width=2000, 
                height=1500)
ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
          shape = "DiabetesSubtype",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "mypal3",rug = T,
          title = "PCA plot of sample distances by gene count")
dev.off()

png(file=paste0("PCA_plot_default-DiabetesSubtype_Elipse.png"),
                res=300, 
                width=2000, 
                height=1500)
ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
          shape = "DiabetesSubtype",
          ellipse = T,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "mypal3",rug = T,
          title = "PCA plot of sample distances by gene count")
dev.off()

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "SampleID",
          shape = "SampleID",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = F,palette = "mypal3",rug = F,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-Sample.png"),
                res=300, 
                width=2500, 
                height=3000)
a
dev.off()


df.unk<-df%>%filter(DiabetesSubtype!=c("Control"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("A1"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("A1.Control"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("A2"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("A2.T2DM"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("A1.T2DM"))
df.unk<-df.unk%>%filter(DiabetesSubtype!=c("T2DM"))

glimpse(df.unk)

a<-ggscatterhist(data = df.unk,
          x = "Eig2",
          y = "Eig1",
          color = "SampleID",
          shape = "SampleID",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = F,palette = "mypal3",rug = F,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-Sample_unk.png"),
                res=300, 
                width=2500, 
                height=3000)
a
dev.off()

a<-ggscatterhist(data = df.unk,
          x = "Eig2",
          y = "Eig1",
          color = "SampleID",
          shape = "SampleID",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = F,palette = "mypal3",rug = F,
           label = "SampleID", repel = TRUE,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-Sample_unk2.png"),
                res=300, 
                width=2500, 
                height=3000)
a
dev.off()


a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
          shape = "DiabetesSubtype",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = F,palette = "mypal3",rug = F,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-DiabetesSubtype_v2.png"),
                res=300, 
                width=1500, 
                height=1500)
a
dev.off()

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "DiabetesSubtype",
          shape = "DiabetesSubtype",
          ellipse = T,
          ellipse.alpha = 0,
          ellipse.type = F,palette = "mypal3",rug = F,
          title = "PCA plot of sample distances by gene count")

## save PCA plot 
png(file=paste0("PCA_plot_default-DiabetesSubtype_v2-elipse.png"),
                res=300, 
                width=1500, 
                height=1500)
a
dev.off()


# Add vertical and horizontal line to a ggscatterhist
plots <- ggscatterhist(df, x = "Eig2", y = "Eig1",
                       margin.params = list(theme= theme_pubr()),
                       color="DiabetesSubtype",
                       palette = "mypal3",print = FALSE)
plots$sp<-plots$sp+
  stat_density_2d(aes(colour = DiabetesSubtype,size=0.05),contour = T,alpha=0.2,size=0.05,
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
# ggpar(p = plots,legend.title = "Group") ## this command will change the legend from DiabetesSubtype to whatever you want
dev.off()


### Make a custom plot 2
tally(~DiabetesSubtype,df)
data.frame(df)
# customise plot
customised_plot <- 
  ggplot(df,mapping = aes(x = Eig2,y = Eig1, z=0, color=DiabetesSubtype)) +
  geom_point()+
  geom_contour() +
scale_colour_brewer(palette = "Set1") +
  theme_pubr(legend = "bottom")+xlim(c(-0.1,0.2))+ylim(-0.2,0.2)+
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
                                       design = ~DiabetesSubtype) 
## note- R will not like dash in the A1-T2DM column names

getwd()
dds <-DESeq(deseq_counts,parallel = T)
## Make an elbowplot to determine how many PCs in a PCA = > 85% of variance. This will determine cutreen=
## https://support.bioconductor.org/p/83626/
     library(matrixStats)
     
     #How to get PCA plot?
     ##how to obtain d.deseq was described in DESeq2 manual
     
     cds=estimateDispersions(dds)
     vsd=varianceStabilizingTransformation(cds)
     plotPCA(vsd,intgroup=c("DiabetesSubtype"))
     p1 <-      plotPCA(vsd,intgroup=c("DiabetesSubtype"))
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
     scree_plot[,2]<- c(1:57)
     p2 <-      scree_plot[,2]<- c(1:57)


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
# save.image("gdm.placenta_data-v3.RData")

setwd("..")
###### DESeq2 will have 3 iterations ctrl vs GDMA1 or GDMA2 or T2DM
#########################################################################################################
##########################    Iteration I:    Control-vs-GDMA1DiabetesSubtype Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")
# load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/gdm.placenta_data-v3.RData")

dir.create("Control-vs-GDMA1DiabetesSubtype")
setwd("Control-vs-GDMA1DiabetesSubtype")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("DiabetesSubtype","Control","A1"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the A1 condition
  ## Wald hypothesis testing results in:
LFC > 0 (up)       : 4287, 12%
LFC < 0 (down)     : 2467, 7.1%
outliers [1]       : 0, 0%
low counts [2]     : 9377, 27%
(mean count < 1)
write.csv(res, file = "ControlDiabetesSubtype_A1_DESeq2_allResults-Wald.csv")
  ## Same results as before when I had dashes 
    ## same results with cooksCutoff=T

res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$padj),] # order genes by significance
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
res05
    # 6754
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
res05
  ## 9582

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)

head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 2, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < -2, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue','padj')]


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
# dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","A1"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
# summary(res)
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7725, 22%
LFC < 0 (down)     : 5022, 14%
outliers [1]       : 0, 0%
low counts [2]     : 8037, 23%
(mean count < 0)
# write.csv(res, file = "ControlDiabetesSubtype_A1_DESeq2_allResults-LTR.csv")
  ## Same as before, also the same with cooksCutoff=T


######################################################################################### running manual heatmaps for now ###################################################################
########## Generate a heatmap of the counts with only the Control and GDMA1 DiabetesSubtype samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 6,754 sig genes
A1.genes.wald <- res2$gene

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
  ## df4 is 36,621 × 58
  ## res2 is 6,754 × 7
    # A tibble: 385,662 × 3 That doesn't seem right?
res3
res4<-full_join(res3,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
res4 # 385,662 × 4
res5<-res4%>%filter(DiabetesSubtype%in%c("Control","A1"))%>%distinct_all()
res5 #108,128 × 4

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
sig2<-full_join(sig,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
sig3<-sig2%>%filter(DiabetesSubtype%in%c("Control","A1"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes) # 6754 - same as A1.genes.wald above
top.genes.Control.A1 <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = SampleID,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 6754   16
df_slim<-df%>%
  select(SampleID,DiabetesSubtype)%>%
  distinct_all()%>%
  filter(DiabetesSubtype%in%c("Control","A1"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(DiabetesSubtype)
rownames(anot_col)<-df_slim$SampleID

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


#########################################################################################  manual heatmaps done now ###################################################################

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("DiabetesSubtype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype, shape=DiabetesSubtype)) +

p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype)) +
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
png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='DiabetesSubtype')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_p-hist_smooth.png"),
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
sidecols <- c("grey","dodgerblue")[ vst$A1 ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$A1,"-",vst$DiabetesSubtype)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
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

png(file=paste0("ControlDiabetesSubtype_GDMA1DiabetesSubtype_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
png(file=paste0("Top20-GDMA1DiabetesSubtype_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) 
dev.off()

setwd("..")
#########################################################################################################
##########################    Iteration II:    Control-vs-A2DiabetesSubtype Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")
# load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/gdm.placenta_data-v3.RData")
dir.create("Control-vs-A2DiabetesSubtype")
setwd("Control-vs-A2DiabetesSubtype")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("DiabetesSubtype","Control","A2"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the A2 condition
  ## Wald hypothesis testing results in:
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1678, 4.8%
LFC < 0 (down)     : 1637, 4.7%
outliers [1]       : 0, 0%
low counts [2]     : 10047, 29%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
write.csv(res, file = "ControlDiabetesSubtype_A2_DESeq2_allResults-Wald.csv")
  ## Same results as before when I had dashes 
    ## same results with cooksCutoff=T

res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$padj),] # order genes by significance
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
res05
    # 6754
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
res05
  ## 9582

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)

head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 2, ][1:n,], resOrdered[ resOrdered[,'log2FoldChange'] < -2, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue','padj')]


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
# dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","A2"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
# summary(res)
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7725, 22%
LFC < 0 (down)     : 5022, 14%
outliers [1]       : 0, 0%
low counts [2]     : 8037, 23%
(mean count < 0)
# write.csv(res, file = "ControlDiabetesSubtype_A2_DESeq2_allResults-LTR.csv")
  ## Same as before, also the same with cooksCutoff=T


######################################################################################### running manual heatmaps for now ###################################################################
########## Generate a heatmap of the counts with only the Control and A2 DiabetesSubtype samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 6,754 sig genes
A2.genes.wald <- res2$gene

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
  ## df4 is 36,621 × 58
  ## res2 is 6,754 × 7
    # A tibble: 385,662 × 3 That doesn't seem right?
res3
res4<-full_join(res3,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
res4 # 385,662 × 4
res5<-res4%>%filter(DiabetesSubtype%in%c("Control","A2"))%>%distinct_all()
res5 #108,128 × 4

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
sig2<-full_join(sig,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
sig3<-sig2%>%filter(DiabetesSubtype%in%c("Control","A2"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes) # 6754 - same as A2.genes.wald above
top.genes.Control.A2 <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = SampleID,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 6754   16
df_slim<-df%>%
  select(SampleID,DiabetesSubtype)%>%
  distinct_all()%>%
  filter(DiabetesSubtype%in%c("Control","A2"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(DiabetesSubtype)
rownames(anot_col)<-df_slim$SampleID

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
png(file=paste0("Control-VS-A2_counts-heatmap_row-scaled_notcut.png"),
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


#########################################################################################  manual heatmaps done now ###################################################################

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("DiabetesSubtype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype, shape=DiabetesSubtype)) +

p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype)) +
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
png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='DiabetesSubtype')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_p-hist_smooth.png"),
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
sidecols <- c("grey","dodgerblue")[ vst$A2 ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$A2,"-",vst$DiabetesSubtype)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
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

png(file=paste0("ControlDiabetesSubtype_A2DiabetesSubtype_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
png(file=paste0("Top20-A2DiabetesSubtype_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) 
dev.off()

setwd("..")
#########################################################################################################
##########################    Iteration III:    Control-vs-T2DMDiabetesSubtype Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")
# load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/gdm.placenta_data-v3.RData")
dir.create("Control-vs-T2DMDiabetesSubtype")
setwd("Control-vs-T2DMDiabetesSubtype")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("DiabetesSubtype","Control","T2DM"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the T2DM condition
  ## Wald hypothesis testing results in:
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3382, 9.8%
LFC < 0 (down)     : 2607, 7.5%
outliers [1]       : 0, 0%
low counts [2]     : 9377, 27%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
write.csv(res, file = "ControlDiabetesSubtype_T2DM_DESeq2_allResults-Wald.csv")
  ## Same results as before when I had dashes 
    ## same results with cooksCutoff=T

res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$padj),] # order genes by significance
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
res05
    # 6754
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
res05
  ## 9582

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)

head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 2, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < -2, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue','padj')]


# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
# dds <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
# res <- results(dds, contrast=c("DiabetesSubtype","Control","T2DM"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
# summary(res)
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7725, 22%
LFC < 0 (down)     : 5022, 14%
outliers [1]       : 0, 0%
low counts [2]     : 8037, 23%
(mean count < 0)
# write.csv(res, file = "ControlDiabetesSubtype_T2DM_DESeq2_allResults-LTR.csv")
  ## Same as before, also the same with cooksCutoff=T


######################################################################################### running manual heatmaps for now ###################################################################
########## Generate a heatmap of the counts with only the Control and T2DM DiabetesSubtype samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 6,754 sig genes
T2DM.genes.wald <- res2$gene

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
  ## df4 is 36,621 × 58
  ## res2 is 6,754 × 7
    # A tibble: 385,662 × 3 That doesn't seem right?
res3
res4<-full_join(res3,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
res4 # 385,662 × 4
res5<-res4%>%filter(DiabetesSubtype%in%c("Control","T2DM"))%>%distinct_all()
res5 #108,128 × 4

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
sig2<-full_join(sig,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
sig3<-sig2%>%filter(DiabetesSubtype%in%c("Control","T2DM"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes) # 6754 - same as T2DM.genes.wald above
top.genes.Control.T2DM <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = SampleID,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 6754   16
df_slim<-df%>%
  select(SampleID,DiabetesSubtype)%>%
  distinct_all()%>%
  filter(DiabetesSubtype%in%c("Control","T2DM"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(DiabetesSubtype)
rownames(anot_col)<-df_slim$SampleID

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


#########################################################################################  manual heatmaps done now ###################################################################

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("DiabetesSubtype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype, shape=DiabetesSubtype)) +

p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype)) +
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
png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='DiabetesSubtype')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_p-hist_smooth.png"),
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
colnames(mat) <- paste0(vst$T2DM,"-",vst$DiabetesSubtype)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
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

png(file=paste0("ControlDiabetesSubtype_T2DMDiabetesSubtype_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
png(file=paste0("Top20-T2DMDiabetesSubtype_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) 
dev.off()

setwd("..")


#########################################################################################################
##########################    Iteration IV:    Control-vs-UNKDiabetesSubtype Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")
# load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/gdm.placenta_data-v3.RData")
dir.create("Control-vs-UNKDiabetesSubtype")
setwd("Control-vs-UNKDiabetesSubtype")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("DiabetesSubtype","Control","UNK"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the UNK condition
  ## Wald hypothesis testing results in:
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 8224, 24%
LFC < 0 (down)     : 4177, 12%
outliers [1]       : 0, 0%
low counts [2]     : 6698, 19%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
write.csv(res, file = "ControlDiabetesSubtype_UNK_DESeq2_allResults-Wald.csv")


res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$padj),] # order genes by significance
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
res05
    # 12401
res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
res05
  ## 14355

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)

head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 2, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < -2, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue','padj')]




######################################################################################### running manual heatmaps for now ###################################################################
########## Generate a heatmap of the counts with only the Control and UNK DiabetesSubtype samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2
  ## 6,754 sig genes
UNK.genes.wald <- res2$gene

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
  ## df4 is 36,621 × 58
  ## res2 is 6,754 × 7
    # A tibble: 385,662 × 3 That doesn't seem right?
res3
res4<-full_join(res3,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
res4 # 385,662 × 4
res5<-res4%>%filter(DiabetesSubtype%in%c("Control","UNK"))%>%distinct_all()
res5 #108,128 × 4

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
sig2<-full_join(sig,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
sig3<-sig2%>%filter(DiabetesSubtype%in%c("Control","UNK"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

top.genes <- sig3$gene
top.genes <- unique(top.genes) # 6754 - same as UNK.genes.wald above
top.genes.Control.UNK <- top.genes

sig4<-sig3%>%
  pivot_wider(names_from = SampleID,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 6754   16
df_slim<-df%>%
  select(SampleID,DiabetesSubtype)%>%
  distinct_all()%>%
  filter(DiabetesSubtype%in%c("Control","UNK"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(DiabetesSubtype)
rownames(anot_col)<-df_slim$SampleID

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
png(file=paste0("Control-VS-UNK_counts-heatmap_row-scaled_notcut.png"),
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


#########################################################################################  manual heatmaps done now ###################################################################

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("DiabetesSubtype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype, shape=DiabetesSubtype)) +

p2 <- ggplot(pcaData, aes(PC1, PC2, color=DiabetesSubtype)) +
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
png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='DiabetesSubtype')
dev.off()

### Shows dispersion of all data
png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_p-hist_smooth.png"),
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
sidecols <- c("grey","dodgerblue")[ vst$UNK ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$UNK,"-",vst$DiabetesSubtype)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
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

png(file=paste0("ControlDiabetesSubtype_UNKDiabetesSubtype_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
png(file=paste0("Top20-UNKDiabetesSubtype_vst-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=TRUE) 
dev.off()

setwd("..")

#########################################################################################################
##########################    Iteration IV, V, and VI:    Control-vs-A1-2TDMDiabetesSubtype Analysis       ###########################################
#########################################################################################################
## NVM WE ONLY HAVE 1 SAMPLE MARKED AS A1-T2DM, 1 marked as A2-T2DM, and 2 as A1-Control
#########################################################################################################
##########################    Iteration final:    GDM-SigPooled Analysis       ###########################################
#########################################################################################################
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2")
dir.create("GDM-SigPooled")
setwd("GDM-SigPooled")

## Combine the list of genes significantly DE for GDM relative to all controls

A1.genes.wald
A2.genes.wald
T2DM.genes.wald
GDM.genes <- union(A1.genes.wald, A2.genes.wald)
GDM.genes <- union(GDM.genes, T2DM.genes.wald)
GDM.genes <- unique(GDM.genes)
GDM.genes # 8749
write.csv(GDM.genes, file = "merged-sig.GDM-genes.csv")

A1A2.genes <- union(A1.genes.wald, A2.genes.wald)
A1A2.genes <- unique(A1A2.genes)
T2DM.signature <- setdiff(T2DM.genes.wald, A1A2.genes)
T2DM.signature # 1520
write.csv(T2DM.signature, file = "T2DM.signature.csv")


T2DMA2.genes <- union(T2DM.genes.wald, A2.genes.wald)
T2DMA2.genes <- unique(T2DMA2.genes)
A1.signature <- setdiff(A1.genes.wald, T2DMA2.genes)
A1.signature #2052
write.csv(A1.signature, file = "A1.signature.csv")

T2DMA1.genes <- union(T2DM.genes.wald, A1.genes.wald)
T2DMA1.genes <- unique(T2DMA1.genes)
A2.signature <- setdiff(A2.genes.wald, T2DMA1.genes)
A2.signature # 267
write.csv(A2.signature, file = "A2.signature.csv")


UNK.signature <- setdiff(UNK.genes.wald, GDM.genes)
UNK.signature # 6001
write.csv(UNK.signature.signature, file = "UNK.signature.csv")


library(VennDiagram) # install.packages("VennDiagram")
venn_data <- list(
  A1 = A1.genes.wald,
  A2 = A2.genes.wald,
  T2DM = T2DM.genes.wald
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1", "A2", "T2DM"),
  filename = NULL  # Set this to a filename if you want to save the plot as an image file
)

grid.draw(venn.plot)

png(file=paste0("Venn.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()



# Calculate overlapping and unique genes
overlap_1_2 <- intersect(A1.genes.wald, A2.genes.wald)
write.csv(overlap_1_2, "overlap_A1_A2.csv", row.names = F)
overlap_1_3 <- intersect(A1.genes.wald, T2DM.genes.wald)
write.csv(overlap_1_3, "overlap_A1_T2DM.csv", row.names = F)
overlap_2_3 <- intersect(A2.genes.wald, T2DM.genes.wald)
write.csv(overlap_2_3, "overlap_A2_T2DM.csv", row.names = F)
overlap_all <- intersect(overlap_1_2, T2DM.genes.wald)
write.csv(overlap_all, "overlap_A1_A2_T2DM.csv", row.names = F)

unique_1 <- setdiff(A1.genes.wald, c(overlap_1_2, overlap_1_3))
write.csv(unique_1, "unique_A1.csv", row.names = F)
unique_2 <- setdiff(A2.genes.wald, c(overlap_1_2, overlap_2_3))
write.csv(unique_2, "unique_A2.csv", row.names = F)
unique_3 <- setdiff(T2DM.genes.wald, c(overlap_1_3, overlap_2_3))
write.csv(unique_3, "unique_T2DM.csv", row.names = F)


# Create a data frame with the results
overlap_table <- data.frame(
  Overlap_A1_A2 = paste(overlap_1_2, collapse = ", "),
  Overlap_A1_T2DM = paste(overlap_1_3, collapse = ", "),
  Overlap_A2_T2DM = paste(overlap_2_3, collapse = ", "),
  Overlap_All = paste(overlap_all, collapse = ", "),
  Unique_A1 = paste(unique_1, collapse = ", "),
  Unique_A2 = paste(unique_2, collapse = ", "),
  Unique_T2DM = paste(unique_3, collapse = ", ")
)

# Transpose the overlap_table to have columns as rows and rows as columns
overlap_table <- as.data.frame(t(overlap_table))

# Write the overlap_table data frame to a CSV file
write.csv(overlap_table, "overlap_table.csv", row.names = TRUE)





library(VennDiagram) # install.packages("VennDiagram")
venn_data <- list(
  A1 = A1.genes.wald,
  A2 = A2.genes.wald,
  T2DM = T2DM.genes.wald, 
  UNK = UNK.genes.wald
)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("A1", "A2", "T2DM", "UNK"),
  filename = NULL  # Set this to a filename if you want to save the plot as an image file
)

grid.draw(venn.plot)

png(file=paste0("Venn_UNK.png"),
                res=300, 
                width=1500, 
                height=1500)
grid.draw(venn.plot)
dev.off()

# Calculate overlapping and unique genes
overlap_1_2 <- intersect(A1.genes.wald, A2.genes.wald)
write.csv(overlap_1_2, "overlap_A1_A2-unk.csv", row.names = F)
overlap_1_3 <- intersect(A1.genes.wald, T2DM.genes.wald)
write.csv(overlap_1_3, "overlap_A1_T2DM-unk.csv", row.names = F)
overlap_1_4 <- intersect(A1.genes.wald, UNK.genes.wald)
write.csv(overlap_1_4, "overlap_A1_UNK-unk.csv", row.names = F)
overlap_2_3 <- intersect(A2.genes.wald, T2DM.genes.wald)
write.csv(overlap_2_3, "overlap_A2_T2DM-unk.csv", row.names = F)
overlap_2_4 <- intersect(A2.genes.wald, UNK.genes.wald)
write.csv(overlap_2_4, "overlap_A2_UNK-unk.csv", row.names = F)
overlap_3_4 <- intersect(T2DM.genes.wald, UNK.genes.wald)
write.csv(overlap_3_4, "overlap_T2DM_UNK-unk.csv", row.names = F)
overlap_all <- intersect(overlap_1_2, overlap_1_3)
overlap_all <- intersect(overlap_all, overlap_1_4)
overlap_all <- intersect(overlap_all, overlap_2_3)
overlap_all <- intersect(overlap_all, overlap_2_4)
overlap_all <- intersect(overlap_all, overlap_3_4)

write.csv(overlap_all, "overlap_all-unk.csv", row.names = F)

  ## overlap_all is empty
unique_1 <- setdiff(A1.genes.wald, c(overlap_1_2, overlap_1_3, overlap_1_4))
write.csv(unique_1, "unique_A1-unk.csv", row.names = F)
unique_2 <- setdiff(A2.genes.wald, c(overlap_1_2, overlap_2_3, overlap_2_4))
write.csv(unique_2, "unique_A2-unk.csv", row.names = F)

unique_3 <- setdiff(T2DM.genes.wald, c(overlap_1_3, overlap_2_3, overlap_3_4))
write.csv(unique_3, "unique_T2DM-unk.csv", row.names = F)

unique_4 <- setdiff(UNK.genes.wald, c(overlap_1_4, overlap_2_4, overlap_3_4))
write.csv(unique_4, "unique_UNK-unk.csv", row.names = F)


# Create a data frame with the results
overlap_table <- data.frame(
  Overlap_A1_A2 = paste(overlap_1_2, collapse = ", "),
  Overlap_A1_T2DM = paste(overlap_1_3, collapse = ", "),
  Overlap_A1_UNK = paste(overlap_1_4, collapse = ", "),
  Overlap_A2_T2DM = paste(overlap_2_3, collapse = ", "),
  Overlap_A2_UNK = paste(overlap_2_4, collapse = ", "),
  Overlap_T2DM_UNK = paste(overlap_3_4, collapse = ", "),
 # Overlap_All = paste(overlap_all, collapse = ", "),
  Unique_A1 = paste(unique_1, collapse = ", "),
  Unique_A2 = paste(unique_2, collapse = ", "),
  Unique_T2DM = paste(unique_3, collapse = ", "),
  Unique_UNK = paste(unique_4, collapse = ", ")
)

# Transpose the overlap_table to have columns as rows and rows as columns
overlap_table <- as.data.frame(t(overlap_table))

# Write the overlap_table data frame to a CSV file
write.csv(overlap_table, "overlap_table_unk.csv", row.names = TRUE)






## add a Diabetes column with Yes or No
meta<-as_tibble(read.table("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/meta.tsv",sep = "\t",header = T))
keeps <- c("SampleID","DiabetesSubtype", "FetusBiologicalSex1F2M", "Diabetes")
meta <- meta[keeps]
tally(~SampleID,meta)
data.frame(meta)
write.table(meta,"meta_final.tsv",sep = "\t",row.names = F)
meta<-as_tibble(read.table("meta_final.tsv",sep = "\t",header = T))

dd_meta<-meta%>%
  # filter(Sample!="JMR16"|Tx!="HFRes2"|rep!=2)%>%  ### no need to remove samples from the metadata. They were removed upstream in this version. 
  #filter(Sample!="JMR30"|Tx!="HFD"|rep!=1)%>%
  mutate(across(where(is.character),factor))
dd_meta

## Finally, load data and metadata into DESeq2
deseq_counts <- DESeqDataSetFromMatrix(countData = df5.1, 
                                       colData =dd_meta,
                                       design = ~Diabetes) 
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Diabetes","Yes","No"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
              ## DE values for yes, relative to No controls
                
summary(res) ## Filtering adj. p < 0.05
out of 34648 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3240, 9.4%
LFC < 0 (down)     : 5121, 15%
outliers [1]       : 0, 0%
low counts [2]     : 6698, 19%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
write.csv(res, file = "Diabetes_YesrelativetoNo_DESeq2_allResults-Wald.csv")


########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA


resOrdered <- res[order(res$padj),]
res05 <- res
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
# 8361
head(resOrdered)
log2 fold change (MLE): Diabetes Yes vs No 
Wald test p-value: Diabetes Yes vs No 
DataFrame with 6 rows and 6 columns


n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue','padj')]

########## Generate a heatmap of with all samples and pooled GDM significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
    ## 8,361 sig genes

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
res3
res4<-full_join(res3,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
res4
res5<-res4%>%filter(DiabetesSubtype%in%c("Control", "A1.Control", "A1", "A2", "A2.T2DM", "T2DM"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
sig2<-full_join(sig,df%>%select(SampleID,DiabetesSubtype)%>%distinct_all())
sig3<-sig2%>%filter(DiabetesSubtype%in%c("Control", "A1.Control", "A1", "A2", "A2.T2DM", "T2DM"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  

## stopping here 8/15 12:06 pm

sig4<-sig3%>%
  pivot_wider(names_from = SampleID,values_from = count) # %>%select(-JMR30)%>%filter(!is.na(gene))


h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-GDM.genes
dim(h)      # [1] 30 53
df_slim<-df%>%
  select(Sample,DiabetesSubtype)%>%
  #distinct_all()%>%
  #filter(DiabetesSubtype%in%c("Control","GDMA1", "GDMA2", "T2DM"))%>%
  distinct_all()
df_slim     # 53 × 2
df
anot_col<-data.frame(df_slim)%>%select(DiabetesSubtype)
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
colnames(mat) <- paste0(vst$Control,"-",vst$DiabetesSubtype)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)    # 31
#reset par
par(mfrow=c(1,1))

png(file=paste0("Basic_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,5)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()
sessionInfo()

# save.image("gdm.placenta_data-final.RData")



####### 
sessionInfo()

setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/GDM-SigPooled")
load("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/GDM-SigPooled/gdm.placenta_data-final.RData")


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
dir.create("functional-enrichment_v1")
setwd("/home/ebarrozo/gdm.placenta/results/DESeq2_analysis_v3.2/functional-enrichment_v1")

# deseq_counts <- DESeqDataSetFromMatrix(countData = df5.1, 
 #                                      colData =dd_meta,
 #                                      design = ~Diabetes) 
# dds <-DESeq(deseq_counts,parallel = T)
# res <- results(dds, contrast=c("Diabetes","Yes","No"), alpha=0.05, cooksCutoff=FALSE)

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = df5.1, 
                                      colData =dd_meta,
                                      design = ~Diabetes)
view(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=" normalized_counts.txt", sep="\t", quote=F, col.names=NA)

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


plotPCA(rld, intgroup="DiabetesSubtype", ntop=3000, ellipse=T)

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
anot_col<-data.frame(df)%>%select(DiabetesSubtype)
rownames(anot_col)<-df$SampleID

pheatmap(mat = rld_cor, annotation_col = anot_col, cutree_rows = 6,cutree_cols = 5,border_color = "black")

## Save with the code
png(file=paste0("Heatmap_DiabetesSubtype-SampleID_pairwisecorrelation-cutmodules.png"),res=300, width=3000, height=2000)
pheatmap(mat = rld_cor, annotation_col = anot_col,cutree_rows = 3,cutree_cols = 7,border_color = "black")
dev.off()



# dds <- DESeq(dds)
dds <-DESeq(dds,parallel = T)


plotDispEsts(dds)

# res <- results(dds)
res <- results(dds, contrast=c("Diabetes","Yes","No"))
summary(res)
out of 34648 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3835, 11%
LFC < 0 (down)     : 6656, 19%
outliers [1]       : 0, 0%
low counts [2]     : 6698, 19%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

plotMA(res, ylim=c(-2,2)) # The genes that are significantly DE are colored to be easily identified.


############################################## Functional analysis with clusterProfiler ##############################################
# Over-representation analysis with clusterProfiler


# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
  ## Add fold-change filter?
# signif_res <- res[res$padj < 0.05 & !is.na(res$padj) & res$log2FoldChange > 2 & !is.na(res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))

# Run GO enrichment analysis
ego <- enrichGO(gene = signif_genes, universe = all_genes,
keyType = "ENSEMBL",
OrgDb = org.Hs.eg.db,
ont = "BP",
pAdjustMethod = "BH",
qvalueCutoff = 0.05,
readable = TRUE)

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

write.table(cluster_summary, file="GOEnrichmentAnalysis_DiabetesYesNo_BH.txt", sep="\t", quote=F, col.names=NA)

# Visualizing clusterProfiler results
# The dotplot shows the number of genes associated with the first 50 terms (size) and the p-adjusted values for these terms (color). This plot displays the top 50 genes by gene ratio (# genes related to GO term / total number of sig genes), not p-adjusted value.
dotplot(ego, showCategory=10)


png(file=paste0("enrichGO-dotplot-top10.png"),res=300, width=2000, height=1200)
dotplot(ego, showCategory=10)
dev.off()

png(file=paste0("enrichGO-dotplot-top25.png"),res=300, width=2500, height=3500)
dotplot(ego, showCategory=25)
dev.off()

emapplot(ego, showCategory=50)
  # Results in an error. 

## troubleshooting with https://github.com/YuLab-SMU/enrichplot/issues/79
library(enrichplot) # BiocManager::install("enrichplot")

x2 <- pairwise_termsim(ego) 
emapplot(x2, showCategory=10)

emapplot(x2, cex_category=1.5)

png(file=paste0("enrichGO-emaplot-top25.png"),res=300, width=3000, height=3000)
emapplot(x2, cex_category=1.5) 
dev.off()


# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)
## 8361 genes

signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
  ##  2753 genes log(2)=0.301
signif_genes <- as.character(rownames(signif_res))


signif_res_lFC <- signif_res$log2FoldChange
cnetplot(ego,
categorySize="pvalue",
showCategory = 5,
foldChange= signif_res_lFC, vertex.label.font=6)

signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
  ##  2753 genes log(2)=0.301
signif_genes <- as.character(rownames(signif_res))

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

signif_res <- signif_res[signif_res$log2FoldChange > 1 & !is.na(signif_res$log2FoldChange), ]
  ##  201 genes log(X)=1
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

signif_res <- signif_res[signif_res$log2FoldChange > 2 & !is.na(signif_res$log2FoldChange), ]
  ##  14 genes 
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

######################################################  too many genes for the cneplot- let's run again with the top 201
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
  ##  22753 genes log(2)=0.301
  ## or filter by the top 200
# signif_res <- names(signif_res)[1:200]
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

## subset the top 200 genes for visualization
signif_res <- signif_res[1:200,]
signif_res_lFC <- signif_res$log2FoldChange

cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE) 

categorys <- c("macromolecule methylation", "methylation",
                   "histone modification", "protein methylation", "protein alkylation")
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 

##### Plot the top 200 genes in the "methylation" node categories
png(file=paste0("cnetplot-methylation-top200_circle.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = TRUE, colorEdge = TRUE, showCategory = categorys) 
dev.off()
png(file=paste0("cnetplot-methylation-top200.png"),res=300, width=3000, height=3000)
cnetplot(ego, foldChange=signif_res_lFC, circular = F, colorEdge = TRUE, showCategory = categorys) 
dev.off()


png(file=paste0("treeplot_q.05.png"),res=300, width=3000, height=3000)
treeplot(x2, hclust_method = "average")
dev.off()
# while (!is.null(dev.list()))  dev.off()


https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
#library(upsetplot) # BiocManager::install("upsetplot")
#upsetplot(x2)

# ridgeplot(x2)


# Optional: gene set enrichment analysis (GSEA) using clusterProfiler and Pathview
library(biomaRt) # BiocManager::install("biomaRt")


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters="ensembl_gene_id",
attributes=c("ensembl_gene_id", "entrezgene_id"), values= all_genes,
mart=mart)
indNA = which(is.na(genes$entrezgene_id))
genes_noNA <- genes[-indNA,]
indnodup = which(duplicated(genes_noNA$ entrezgene_id) == F)
genes_noNA_nodup <- genes_noNA[indnodup,]
lFC <- res$log2FoldChange[-indNA]
lFC <- lFC[indnodup]
names(lFC) <- genes_noNA_nodup$entrezgene_id
# Sort fold changes in decreasing order
lFC <- sort(lFC, decreasing = TRUE)

gseaKEGG <- gseKEGG(geneList = lFC,
organism = "hsa", nPerm = 1000, # default number permutations
minGSSize = 5, # minimum gene set size
pvalueCutoff = 0.1, # padj cutoff value
verbose = FALSE)
# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result



signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
# signif_res <- signif_res[signif_res$log2FoldChange > 0.301 & !is.na(signif_res$log2FoldChange), ]
signif_genes <- as.character(rownames(signif_res))
signif_res_lFC <- signif_res$log2FoldChange

names(signif_res_lFC) <- row.names(signif_res)

signif_res_lFC <- sort(signif_res_lFC, decreasing = TRUE)




mart <- useDataset("hsapiens_gene_ensembl", useMart("SYMBOL"))
genes <- getBM(filters="ensembl_gene_id",
attributes=c("ensembl_gene_id", "entrezgene_id"), values= all_genes,
mart=mart)
indNA = which(is.na(genes$entrezgene_id))
genes_noNA <- genes[-indNA,]
indnodup = which(duplicated(genes_noNA$ entrezgene_id) == F)
genes_noNA_nodup <- genes_noNA[indnodup,]
lFC <- res$log2FoldChange[-indNA]
lFC <- lFC[indnodup]
names(lFC) <- genes_noNA_nodup$entrezgene_id
# Sort fold changes in decreasing order
lFC <- sort(lFC, decreasing = TRUE)



gseaKEGG <- gseKEGG(geneList = signif_genes,
organism = "hsa", nPerm = 1000, # default number permutations
minGSSize = 5, # minimum gene set size
pvalueCutoff = 0.1, # padj cutoff value
verbose = T)
# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

# http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
# https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_2021_June/Markdowns/12_Gene_set_testing.pdf
# https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html


## see GO_analysis_v1.R for running this on each subset







