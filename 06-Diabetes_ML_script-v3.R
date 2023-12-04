## Diabetes_ML_script-v3.R
## 8.29.2023

## Michael Jochum, Ph.D. &
# Enrico R. Barrozo, Ph.D. 
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## see import_wrangle_filter_ml_caret.R for ML and data wrangling source, gdm.placenta_bulk-RNA-seq_script_v3.1.R for bulk RNA-seq, or GDM-Atlas_v2.R for single-cell 

## load required packages

library(tidyverse)
library(pheatmap)
library(corrgram)
library(mosaic)
library(vegan)
library(ggpubr)
mypal3 <-get_palette("igv",50)
library(stringr)
library(data.table)
library(caret)
library(caretEnsemble)
library(ggpubr)
library(doParallel)   # install.packages("doParallel", repos="http://R-Forge.R-project.org")
library(AppliedPredictiveModeling) # install.packages("AppliedPredictiveModeling")


setwd("/home/ebarrozo/gdm.placenta/results")
#################################################################################
#################### DATA IMPORT ###############################################
#################################################################################
dir.create("/home/ebarrozo/gdm.placenta/results/ML_v3")
setwd("/home/ebarrozo/gdm.placenta/results/ML_v3")

sc <- readRDS("/home/ebarrozo/gdm.placenta/results/sc.RDS")
df3_notimpute_long_clean <- readRDS("/home/ebarrozo/gdm.placenta/results/df3_notimpute_long_clean.RDS")

df3_notimpute_long_clean%>%arrange(desc(count))



library(ggpubr)
gghistogram(data = df3_notimpute_long_clean,
            x = "count",
            color = "project",bins = 100,
            facet.by = "target")+xscale("log10")

# install.packages("caret", dependencies = T,suggests=T)#save the rds
# install.packages("caretEnsemble", dependencies = T)#save the rds

unique(df3_notimpute_long_clean$target)
df3_notimpute_wide_clean<-df3_notimpute_long_clean%>%
  pivot_wider(names_from = Gene,values_from = count,values_fill = 0)
df3_notimpute_wide_clean
# saveRDS(df3_notimpute_wide_clean,"df3_notimpute_wide_clean.RDS")

df3_notimpute_wide_clean <- readRDS("/home/ebarrozo/gdm.placenta/results/df3_notimpute_wide_clean.RDS")


#########################################
#DATA VISUALIZATION
#########################################
library(caret)
library(caretEnsemble)  # install.packages("caretEnsemble", dependencies = T)
df3_notimpute_wide_clean
A1<-df3_notimpute_wide_clean%>%
  filter(target=="A1")%>%
  mutate(perinatal_GTT=FALSE,
         Meds=FALSE,
         postnatal_GTT=TRUE)
A2<-df3_notimpute_wide_clean%>%
  filter(target=="A2")%>%
  mutate(perinatal_GTT=FALSE,
         Meds=TRUE,
         postnatal_GTT=TRUE)
T2DM<-df3_notimpute_wide_clean%>%
  filter(target=="T2DM")%>%
  mutate(perinatal_GTT=FALSE,
         Meds=TRUE,
         postnatal_GTT=FALSE)
Control<-df3_notimpute_wide_clean%>%
  filter(target=="Control")%>%
  mutate(perinatal_GTT=TRUE,
         Meds=FALSE,
         postnatal_GTT=TRUE)
UNK<-df3_notimpute_wide_clean%>%
  filter(target=="UNK")%>%
  mutate(perinatal_GTT=NA,
         Meds=NA,
         postnatal_GTT=NA)

annot<-full_join(A1,A2)
annot<-full_join(annot,T2DM)
annot<-full_join(annot,Control)
annot<-full_join(annot,UNK)
annot
# annot<-full_join(annot,UNK)


################################################################################ ################################################################################ ################################################################################
################################################################################  iscontrol Model ################################################################################
################################################################################ ################################################################################ ################################################################################
setwd("/home/ebarrozo/gdm.placenta/results/ML_v3")
dir.create("iscontrol")
setwd("iscontrol")

iscontrol<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target!="UNK")%>%
  mutate(target=as.factor(target=="Control"))

tally(annot$target)
#  A1      A2 Control    T2DM     UNK 
#   6      10      13       5      24


unk_preproc<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target=="UNK")%>%
  mutate(target=as.factor(target=="Control"))



prrproc<-preProcess(x = iscontrol,
                    method = c("center","scale", "corr", "zv", "nzv"),
                    verbose = T,
                    outcome ="target")

iscontrol_preproc<-predict(prrproc,iscontrol)
colnames(iscontrol_preproc)


notunk<-annot%>%filter(target!="UNK")
notunk
rownames(iscontrol_preproc)<-notunk$SampleID
rownames(iscontrol_preproc)


unk_preproc<-predict(prrproc,unk_preproc)
UNK
rownames(unk_preproc)<-UNK$SampleID
rownames(unk_preproc)
################################################################################
################################################################################
# MACHINE LEARNING TIME
################################################################################
################################################################################
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(iscontrol_preproc$target, p = .8, 
                                  list = FALSE, 
                                  times = 1)


head(trainIndex)


iscontrol_preproc
Train<-slice_sample(.data = iscontrol_preproc,prop = 0.8,weight_by = target,replace = F)
tally(Train$target)
# FALSE  TRUE 
#   16    11
Test<-iscontrol_preproc%>%
  mutate(SampleID=rownames(iscontrol_preproc))%>%
  filter(!SampleID%in%rownames(Train))%>%
  select(-SampleID)


fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 5,
  allowParallel = T,
  ## repeated ten times
  repeats = 3)

library(doParallel)   # install.packages("doParallel", repos="http://R-Forge.R-project.org")

cl <- makePSOCKcluster(48)
registerDoParallel(cl)

## All subsequent models are then run in parallel
# model <- train(y ~ ., data = training, method = "rf")
Train2<-data.frame(Train)
rownames(Train2)<-rownames(Train)
rownames(Train2)
set.seed(314)
rfFit1 <- train(target ~ ., data = Train2, 
                 method = "rf", 
                 tuneLength = 25,
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = TRUE)
## When you are done:
# stopCluster(cl)
res<-rfFit1


set.seed(314)
gbmFit2 <- train(target ~ ., data = Train2, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneLength=2)
gbmFit2
# Something is wrong; all the Accuracy metric values are missing:


################################################################################
# MODEL EVALUATION
rf_res<-predict(rfFit1, newdata = data.frame(Test))
Test$pred<-rf_res
Test_res<-Test%>%select(target,pred)
Test_res
write.table(Test_res, file="isControl_TestPrediction_key.txt", sep="\t", col.names=T, row.names=F)

rf_res_prob<-predict(rfFit1, newdata = data.frame(Test), type = "prob")
rf_res_prob


rf_res_prob<-rf_res_prob%>%mutate(SampleID=rownames(Test))
rf_res_prob$pred<-rf_res
colnames(rf_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isControlPrediction")
write.table(rf_res_prob, file="isControl_TestPrediction_results.txt", sep="\t", col.names=T, row.names=F)

unk_res<-predict(rfFit1, newdata = data.frame(unk_preproc))
unk_res_prob<-predict(rfFit1, newdata = data.frame(unk_preproc),type="prob")
unk_res_prob<-unk_res_prob%>%mutate(SampleID=rownames(unk_preproc))
unk_res_prob$pred<-unk_res
colnames(unk_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isControlPrediction")
ggscatter(unk_res_prob,x = "ProbabilityTrue",y = "ProbabilityFalse",color = "isControlPrediction")
unk_res_prob
class(unk_res_prob)

write.table(unk_res_prob, file="isControl_UNKPrediction_results.txt", sep="\t", col.names=T, row.names=F)


postResample(pred = Test$pred,obs = Test$target)
#Accuracy    Kappa 
#       1        1
################################################################################ pca_plots.R

annot
library(vegan)
annot_mat<-as.matrix(annot%>%select(where(is.numeric)))
rownames(annot_mat)<-annot$SampleID
dist<-vegdist(annot_mat,method = "euclidean")
distmat<-as.matrix(dist)
distmat

hdist<-hclust(dist,method = "ward.D2")
hdist
ordiplot()

##################
## Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.mds
dune.dis <- dist
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, annot_mat, expand = TRUE)
dune.mds$points

annot3<-annot%>%
  mutate(PCA1=dune.mds$points[,1],
         PCA2=dune.mds$points[,2])
annot3
annot3_slim<-annot3%>%select(SampleID,target,project,PCA1,PCA2)%>%distinct_all()%>%mutate(target2=target=="UNK")
keep<-annot3_slim%>%filter(target2=="TRUE")%>%mutate(target2=SampleID)
other<-annot3_slim%>%filter(target2!="TRUE")%>%mutate(target2=NA)
annot4_slim<-full_join(keep,other)

my_pal<-pal_d3(palette = "cat")
ggscatter(annot4_slim,x = "PCA1",
          y = "PCA2",
          # label = "target2",
          color = "target",
          palette = "startrek",
          rug = F,
          shape = "project",repel = T,
          ggtheme = theme_pubr(legend = "right"),
          group="target")#+
  # geom_label(mapping = aes(label=target2),position = "jitter",)#+
  
  # geom_smooth(method = "glm",mapping = aes(group=target,fill=NA,color=target))+scale_color_aaas()
          # ellipse.alpha = 0,ellipse.level = 0.95,
          # ellipse.type = "norm",
          # ellipse = T)

# pl <- ordiplot(dune.mds, type = "none")
# points(pl, "sites", pch=21, col="red", bg="yellow")
# text(pl, "species", col="blue", cex=0.9)
# fwrite(Train2,"train.tsv",sep = "\t",row.names = T)
# fwrite(data.frame(Test),"test.tsv",sep = "\t",row.names = T)
# saveRDS(rfFit1,"iscontrol_randomForest.RDS")
summary(rfFit1)
rfFit1$bestTune
plot(rfFit1$results)

png(file=paste0("isControl_rfFit1_Results.png"),
                res=300, 
                width=2000, 
                height=1500)
plot(rfFit1$results)
dev.off()


dim(rfFit1$results)
dim(Train)  # 27 x 13452

resamples(rfFit1)
rfFit1$results


################################################################################ stacked_matrix.R
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))
top<-varImp(rfFit1)
# rf variable importance 
# only 20 most important variables shown (out of 13452)
top<-top$importance%>%arrange(desc(Overall))
top<-top%>%mutate(Gene=rownames(top))%>%mutate(Gene=gsub("\\.","-",Gene))

write.table(top, file="isControl_topvariables.txt", sep="\t", col.names=T, row.names=F)


keep<-top[1:4,]
keep


annot
rownames(top)
annot_preproc<-predict(prrproc,annot)
features<-data.frame(annot_preproc[,c("HSPB2","AC110597.1","AC139768.1","LCN1","target")])
as_tibble(features)
features<-features%>%mutate(target=as.factor(target))
as_tibble(features)
f<-featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
f$layout

png(file=paste0("isControl_topvariables.png"),
                res=300, 
                width=2000, 
                height=1500)
featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
dev.off()


################################################################################
setwd("..")
## When you are done:
# stopCluster(cl)

################################################################################ ################################################################################ ################################################################################
################################################################################  isA1 Model ################################################################################
################################################################################ ################################################################################ ################################################################################
setwd("/home/ebarrozo/gdm.placenta/results/ML_v3")
dir.create("isA1")
setwd("isA1")

isA1<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target!="UNK")%>%
  mutate(target=as.factor(target=="A1"))

tally(annot$target)
  #   A1      A2 Control    T2DM     UNK 
  #    6      10      13       5      24 


unk_preproc<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target=="UNK")%>%
  mutate(target=as.factor(target=="A1"))



prrproc<-preProcess(x = isA1,
                    method = c("center","scale", "corr", "zv", "nzv"),
                    verbose = T,
                    outcome ="target")

isA1_preproc<-predict(prrproc,isA1)
colnames(isA1_preproc)


notunk<-annot%>%filter(target!="UNK")
notunk
rownames(isA1_preproc)<-notunk$SampleID
rownames(isA1_preproc)


unk_preproc<-predict(prrproc,unk_preproc)
UNK
rownames(unk_preproc)<-UNK$SampleID
rownames(unk_preproc)
################################################################################
################################################################################
# MACHINE LEARNING TIME
################################################################################
################################################################################
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(isA1_preproc$target, p = .8, 
                                  list = FALSE, 
                                  times = 1)


head(trainIndex)


isA1_preproc
Train<-slice_sample(.data = isA1_preproc,prop = 0.8,weight_by = target,replace = F)
tally(Train$target)
# FALSE  TRUE 
#   21    6
Test<-isA1_preproc%>%
  mutate(SampleID=rownames(isA1_preproc))%>%
  filter(!SampleID%in%rownames(Train))%>%
  select(-SampleID)


fitA1 <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 5,
  allowParallel = T,
  ## repeated ten times
  repeats = 3)

library(doParallel)   # install.packages("doParallel", repos="http://R-Forge.R-project.org")

cl <- makePSOCKcluster(48)
registerDoParallel(cl)

## All subsequent models are then run in parallel
# model <- train(y ~ ., data = training, method = "rf")
Train2<-data.frame(Train)
rownames(Train2)<-rownames(Train)
rownames(Train2)
set.seed(314)
rfFit1 <- train(target ~ ., data = Train2, 
                 method = "rf", 
                 tuneLength = 25,
                 trA1 = fitA1,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = TRUE)
## When you are done:
# stopCluster(cl)
res<-rfFit1


set.seed(314)
gbmFit2 <- train(target ~ ., data = Train2, 
                 method = "gbm", 
                 trA1 = fitA1, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneLength=2)
gbmFit2
# Something is wrong; all the Accuracy metric values are missing:


################################################################################
# MODEL EVALUATION
rf_res<-predict(rfFit1, newdata = data.frame(Test))
Test$pred<-rf_res
Test_res<-Test%>%select(target,pred)
Test_res
write.table(Test_res, file="isA1_TestPrediction_key.txt", sep="\t", col.names=T, row.names=F)

rf_res_prob<-predict(rfFit1, newdata = data.frame(Test), type = "prob")
rf_res_prob


rf_res_prob<-rf_res_prob%>%mutate(SampleID=rownames(Test))
rf_res_prob$pred<-rf_res
colnames(rf_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isA1Prediction")
rf_res_prob
write.table(rf_res_prob, file="isA1_TestPrediction_results.txt", sep="\t", col.names=T, row.names=F)

unk_res<-predict(rfFit1, newdata = data.frame(unk_preproc))
unk_res_prob<-predict(rfFit1, newdata = data.frame(unk_preproc),type="prob")
unk_res_prob<-unk_res_prob%>%mutate(SampleID=rownames(unk_preproc))
unk_res_prob$pred<-unk_res
colnames(unk_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isA1Prediction")
unk_res_prob
ggscatter(unk_res_prob,x = "ProbabilityTrue",y = "ProbabilityFalse",color = "isA1Prediction")
unk_res_prob
class(unk_res_prob)

write.table(unk_res_prob, file="isA1_UNKPrediction_results.txt", sep="\t", col.names=T, row.names=F)


postResample(pred = Test$pred,obs = Test$target)
#Accuracy    Kappa 
#      1       NA
################################################################################ pca_plots.R

annot
library(vegan)
annot_mat<-as.matrix(annot%>%select(where(is.numeric)))
rownames(annot_mat)<-annot$SampleID
dist<-vegdist(annot_mat,method = "euclidean")
distmat<-as.matrix(dist)
distmat

hdist<-hclust(dist,method = "ward.D2")
hdist
ordiplot()

##################
## Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.mds
dune.dis <- dist
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, annot_mat, expand = TRUE)
dune.mds$points

annot3<-annot%>%
  mutate(PCA1=dune.mds$points[,1],
         PCA2=dune.mds$points[,2])
annot3
annot3_slim<-annot3%>%select(SampleID,target,project,PCA1,PCA2)%>%distinct_all()%>%mutate(target2=target=="UNK")
keep<-annot3_slim%>%filter(target2=="TRUE")%>%mutate(target2=SampleID)
other<-annot3_slim%>%filter(target2!="TRUE")%>%mutate(target2=NA)
annot4_slim<-full_join(keep,other)

my_pal<-pal_d3(palette = "cat")
ggscatter(annot4_slim,x = "PCA1",
          y = "PCA2",
          # label = "target2",
          color = "target",
          palette = "startrek",
          rug = F,
          shape = "project",repel = T,
          ggtheme = theme_pubr(legend = "right"),
          group="target")#+
  # geom_label(mapping = aes(label=target2),position = "jitter",)#+
  
  # geom_smooth(method = "glm",mapping = aes(group=target,fill=NA,color=target))+scale_color_aaas()
          # ellipse.alpha = 0,ellipse.level = 0.95,
          # ellipse.type = "norm",
          # ellipse = T)

# pl <- ordiplot(dune.mds, type = "none")
# points(pl, "sites", pch=21, col="red", bg="yellow")
# text(pl, "species", col="blue", cex=0.9)
# fwrite(Train2,"train.tsv",sep = "\t",row.names = T)
# fwrite(data.frame(Test),"test.tsv",sep = "\t",row.names = T)
# saveRDS(rfFit1,"isA1_randomForest.RDS")
summary(rfFit1)
rfFit1$bestTune
plot(rfFit1$results)

png(file=paste0("isA1_rfFit1_Results.png"),
                res=300, 
                width=2000, 
                height=1500)
plot(rfFit1$results)
dev.off()


dim(rfFit1$results)
dim(Train)  # 27 x 13453

resamples(rfFit1)
rfFit1$results


################################################################################ stacked_matrix.R
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))
top<-varImp(rfFit1)
top
# rf variable importance 
# only 20 most important variables shown (out of 13452)
top<-top$importance%>%arrange(desc(Overall))
top<-top%>%mutate(Gene=rownames(top))%>%mutate(Gene=gsub("\\.","-",Gene))

write.table(top, file="isA1_topvariables.txt", sep="\t", col.names=T, row.names=F)


keep<-top[1:4,]
keep


annot
rownames(top)
annot_preproc<-predict(prrproc,annot)
features<-data.frame(annot_preproc[,c("SMAD5-AS1","ENTPD8","CTSG","SIGLEC15","target")]) ## convert . to -
as_tibble(features)
features<-features%>%mutate(target=as.factor(target))
as_tibble(features)
f<-featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
f$layout
f
png(file=paste0("isA1_topvariables.png"),
                res=300, 
                width=2000, 
                height=1500)
featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
dev.off()


################################################################################
setwd("..")
## When you are done:
# stopCluster(cl)
################################################################################ ################################################################################ ################################################################################
################################################################################  isA2 Model ################################################################################
################################################################################ ################################################################################ ################################################################################
setwd("/home/ebarrozo/gdm.placenta/results/ML_v3")
dir.create("isA2")
setwd("isA2")

isA2<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target!="UNK")%>%
  mutate(target=as.factor(target=="A2"))

tally(annot$target)
  #   A1      A2 Control    T2DM     UNK 
  #    6      10      13       5      24 


unk_preproc<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target=="UNK")%>%
  mutate(target=as.factor(target=="A2"))



prrproc<-preProcess(x = isA2,
                    method = c("center","scale", "corr", "zv", "nzv"),
                    verbose = T,
                    outcome ="target")

isA2_preproc<-predict(prrproc,isA2)
colnames(isA2_preproc)


notunk<-annot%>%filter(target!="UNK")
notunk
rownames(isA2_preproc)<-notunk$SampleID
rownames(isA2_preproc)


unk_preproc<-predict(prrproc,unk_preproc)
UNK
rownames(unk_preproc)<-UNK$SampleID
rownames(unk_preproc)
################################################################################
################################################################################
# MACHINE LEARNING TIME
################################################################################
################################################################################
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(isA2_preproc$target, p = .8, 
                                  list = FALSE, 
                                  times = 1)


head(trainIndex)


isA2_preproc
Train<-slice_sample(.data = isA2_preproc,prop = 0.8,weight_by = target,replace = F)
tally(Train$target)
# FALSE  TRUE 
#   18    9
Test<-isA2_preproc%>%
  mutate(SampleID=rownames(isA2_preproc))%>%
  filter(!SampleID%in%rownames(Train))%>%
  select(-SampleID)


fitA2 <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 5,
  allowParallel = T,
  ## repeated ten times
  repeats = 3)

library(doParallel)   # install.packages("doParallel", repos="http://R-Forge.R-project.org")

cl <- makePSOCKcluster(48)
registerDoParallel(cl)

## All subsequent models are then run in parallel
# model <- train(y ~ ., data = training, method = "rf")
Train2<-data.frame(Train)
rownames(Train2)<-rownames(Train)
rownames(Train2)
set.seed(314)
rfFit1 <- train(target ~ ., data = Train2, 
                 method = "rf", 
                 tuneLength = 25,
                 trA2 = fitA2,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = TRUE)
## When you are done:
# stopCluster(cl)
res<-rfFit1


set.seed(314)
gbmFit2 <- train(target ~ ., data = Train2, 
                 method = "gbm", 
                 trA2 = fitA2, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneLength=2)
gbmFit2
# Something is wrong; all the Accuracy metric values are missing:


################################################################################
# MODEL EVALUATION
rf_res<-predict(rfFit1, newdata = data.frame(Test))
Test$pred<-rf_res
Test_res<-Test%>%select(target,pred)
Test_res
write.table(Test_res, file="isA2_TestPrediction_key.txt", sep="\t", col.names=T, row.names=F)

rf_res_prob<-predict(rfFit1, newdata = data.frame(Test), type = "prob")
rf_res_prob


rf_res_prob<-rf_res_prob%>%mutate(SampleID=rownames(Test))
rf_res_prob$pred<-rf_res
colnames(rf_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isA2Prediction")
rf_res_prob
write.table(rf_res_prob, file="isA2_TestPrediction_results.txt", sep="\t", col.names=T, row.names=F)

unk_res<-predict(rfFit1, newdata = data.frame(unk_preproc))
unk_res_prob<-predict(rfFit1, newdata = data.frame(unk_preproc),type="prob")
unk_res_prob<-unk_res_prob%>%mutate(SampleID=rownames(unk_preproc))
unk_res_prob$pred<-unk_res
colnames(unk_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isA2Prediction")
unk_res_prob
ggscatter(unk_res_prob,x = "ProbabilityTrue",y = "ProbabilityFalse",color = "isA2Prediction")
unk_res_prob
class(unk_res_prob)

write.table(unk_res_prob, file="isA2_UNKPrediction_results.txt", sep="\t", col.names=T, row.names=F)


postResample(pred = Test$pred,obs = Test$target)
#Accuracy    Kappa 
#      0.8571429       0.0000000
################################################################################ pca_plots.R

annot
library(vegan)
annot_mat<-as.matrix(annot%>%select(where(is.numeric)))
rownames(annot_mat)<-annot$SampleID
dist<-vegdist(annot_mat,method = "euclidean")
distmat<-as.matrix(dist)
distmat

hdist<-hclust(dist,method = "ward.D2")
hdist
ordiplot()

##################
## Draw a plot for a non-vegan ordination (cmdscale).
data(dune)
dune.mds
dune.dis <- dist
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, annot_mat, expand = TRUE)
dune.mds$points

annot3<-annot%>%
  mutate(PCA2=dune.mds$points[,1],
         PCA2=dune.mds$points[,2])
annot3
annot3_slim<-annot3%>%select(SampleID,target,project,PCA2,PCA2)%>%distinct_all()%>%mutate(target2=target=="UNK")
keep<-annot3_slim%>%filter(target2=="TRUE")%>%mutate(target2=SampleID)
other<-annot3_slim%>%filter(target2!="TRUE")%>%mutate(target2=NA)
annot4_slim<-full_join(keep,other)

my_pal<-pal_d3(palette = "cat")
ggscatter(annot4_slim,x = "PCA2",
          y = "PCA2",
          # label = "target2",
          color = "target",
          palette = "startrek",
          rug = F,
          shape = "project",repel = T,
          ggtheme = theme_pubr(legend = "right"),
          group="target")#+
  # geom_label(mapping = aes(label=target2),position = "jitter",)#+
  
  # geom_smooth(method = "glm",mapping = aes(group=target,fill=NA,color=target))+scale_color_aaas()
          # ellipse.alpha = 0,ellipse.level = 0.95,
          # ellipse.type = "norm",
          # ellipse = T)

# pl <- ordiplot(dune.mds, type = "none")
# points(pl, "sites", pch=21, col="red", bg="yellow")
# text(pl, "species", col="blue", cex=0.9)
# fwrite(Train2,"train.tsv",sep = "\t",row.names = T)
# fwrite(data.frame(Test),"test.tsv",sep = "\t",row.names = T)
# saveRDS(rfFit1,"isA2_randomForest.RDS")
summary(rfFit1)
rfFit1$bestTune
plot(rfFit1$results)

png(file=paste0("isA2_rfFit1_Results.png"),
                res=300, 
                width=2000, 
                height=1500)
plot(rfFit1$results)
dev.off()


dim(rfFit1$results)
dim(Train)  # 27 x 13453

resamples(rfFit1)
rfFit1$results


################################################################################ stacked_matrix.R
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))
top<-varImp(rfFit1)
top
# rf variable importance 
# only 20 most important variables shown (out of 13452)
top<-top$importance%>%arrange(desc(Overall))
top<-top%>%mutate(Gene=rownames(top))%>%mutate(Gene=gsub("\\.","-",Gene))

write.table(top, file="isA2_topvariables.txt", sep="\t", col.names=T, row.names=F)


keep<-top[1:4,]
keep


annot
rownames(top)
annot_preproc<-predict(prrproc,annot)
features<-data.frame(annot_preproc[,c("CLMAT3","UBE2Q1-AS1","ZC3H11A","POU1F1","target")]) ## convert . to - ;; FAM66B POU1F1 SLC34A3 
as_tibble(features)
features<-features%>%mutate(target=as.factor(target))
as_tibble(features)
f<-featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
f$layout
f
png(file=paste0("isA2_topvariables.png"),
                res=300, 
                width=2000, 
                height=1500)
featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
dev.off()


################################################################################
setwd("..")
## When you are done:
# stopCluster(cl)
################################################################################ ################################################################################ ################################################################################
################################################################################  isT2DM Model ################################################################################
################################################################################ ################################################################################ ################################################################################
setwd("/home/ebarrozo/gdm.placenta/results/ML_v3")
dir.create("isT2DM")
setwd("isT2DM")

isT2DM<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target!="UNK")%>%
  mutate(target=as.factor(target=="T2DM"))

tally(annot$target)
  #   A1      A2 Control    T2DM     UNK 
  #    6      10      13       5      24 


unk_preproc<-annot%>%
  select(-c(SampleID,project,perinatal_GTT,Meds,postnatal_GTT))%>%
  filter(target=="UNK")%>%
  mutate(target=as.factor(target=="T2DM"))



prrproc<-preProcess(x = isT2DM,
                    method = c("center","scale", "corr", "zv", "nzv"),
                    verbose = T,
                    outcome ="target")

isT2DM_preproc<-predict(prrproc,isT2DM)
colnames(isT2DM_preproc)


notunk<-annot%>%filter(target!="UNK")
notunk
rownames(isT2DM_preproc)<-notunk$SampleID
rownames(isT2DM_preproc)


unk_preproc<-predict(prrproc,unk_preproc)
UNK
rownames(unk_preproc)<-UNK$SampleID
rownames(unk_preproc)
################################################################################
################################################################################
# MACHINE LEARNING TIME
################################################################################
################################################################################
library(caret)
set.seed(3456)
trainIndex <- createDataPartition(isT2DM_preproc$target, p = .8, 
                                  list = FALSE, 
                                  times = 1)


head(trainIndex)


isT2DM_preproc
Train<-slice_sample(.data = isT2DM_preproc,prop = 0.8,weight_by = target,replace = F)
tally(Train$target)
# FALSE  TRUE 
#   22    5
Test<-isT2DM_preproc%>%
  mutate(SampleID=rownames(isT2DM_preproc))%>%
  filter(!SampleID%in%rownames(Train))%>%
  select(-SampleID)


fitT2DM <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 5,
  allowParallel = T,
  ## repeated ten times
  repeats = 3)

library(doParallel)   # install.packages("doParallel", repos="http://R-Forge.R-project.org")

cl <- makePSOCKcluster(48)
registerDoParallel(cl)

## All subsequent models are then run in parallel
# model <- train(y ~ ., data = training, method = "rf")
Train2<-data.frame(Train)
rownames(Train2)<-rownames(Train)
rownames(Train2)
set.seed(314)
rfFit1 <- train(target ~ ., data = Train2, 
                 method = "rf", 
                 tuneLength = 25,
                 trT2DM = fitT2DM,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = TRUE)
## When you are done:
# stopCluster(cl)
res<-rfFit1


set.seed(314)
gbmFit2 <- train(target ~ ., data = Train2, 
                 method = "gbm", 
                 trT2DM = fitT2DM, 
                 verbose = FALSE, 
                 ## Now specify the exact models 
                 ## to evaluate:
                 tuneLength=2)
gbmFit2
# Something is wrong; all the Accuracy metric values are missing:


################################################################################
# MODEL EVALUATION
rf_res<-predict(rfFit1, newdata = data.frame(Test))
Test$pred<-rf_res
Test_res<-Test%>%select(target,pred)
Test_res
write.table(Test_res, file="isT2DM_TestPrediction_key.txt", sep="\t", col.names=T, row.names=F)

rf_res_prob<-predict(rfFit1, newdata = data.frame(Test), type = "prob")
rf_res_prob


rf_res_prob<-rf_res_prob%>%mutate(SampleID=rownames(Test))
rf_res_prob$pred<-rf_res
colnames(rf_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isT2DMPrediction")
rf_res_prob
write.table(rf_res_prob, file="isT2DM_TestPrediction_results.txt", sep="\t", col.names=T, row.names=F)

unk_res<-predict(rfFit1, newdata = data.frame(unk_preproc))
unk_res_prob<-predict(rfFit1, newdata = data.frame(unk_preproc),type="prob")
unk_res_prob<-unk_res_prob%>%mutate(SampleID=rownames(unk_preproc))
unk_res_prob$pred<-unk_res
colnames(unk_res_prob)<-c("ProbabilityFalse","ProbabilityTrue","SampleID","isT2DMPrediction")
unk_res_prob
ggscatter(unk_res_prob,x = "ProbabilityTrue",y = "ProbabilityFalse",color = "isT2DMPrediction")
unk_res_prob
class(unk_res_prob)

write.table(unk_res_prob, file="isT2DM_UNKPrediction_results.txt", sep="\t", col.names=T, row.names=F)


postResample(pred = Test$pred,obs = Test$target)
#Accuracy    Kappa 
#      1       NA


# pl <- ordiplot(dune.mds, type = "none")
# points(pl, "sites", pch=21, col="red", bg="yellow")
# text(pl, "species", col="blue", cex=0.9)
# fwrite(Train2,"train.tsv",sep = "\t",row.names = T)
# fwrite(data.frame(Test),"test.tsv",sep = "\t",row.names = T)
# saveRDS(rfFit1,"isT2DM_randomForest.RDS")
summary(rfFit1)
rfFit1$bestTune
plot(rfFit1$results)

png(file=paste0("isT2DM_rfFit1_Results.png"),
                res=300, 
                width=2000, 
                height=1500)
plot(rfFit1$results)
dev.off()


dim(rfFit1$results)
dim(Train)  # 27 x 13453

resamples(rfFit1)
rfFit1$results


################################################################################ stacked_matrix.R
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)
library(caret)
featurePlot(x = iris[, 1:4], 
            y = iris$Species, 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 3))
top<-varImp(rfFit1)
top
# rf variable importance 
# only 20 most important variables shown (out of 13452)
top<-top$importance%>%arrange(desc(Overall))
top<-top%>%mutate(Gene=rownames(top))%>%mutate(Gene=gsub("\\.","-",Gene))

write.table(top, file="isT2DM_topvariables.txt", sep="\t", col.names=T, row.names=F)


keep<-top[1:4,]
keep


annot
rownames(top)
annot_preproc<-predict(prrproc,annot)
features<-data.frame(annot_preproc[,c("REG3G","THEMIS","MYRIP","RCVRN","target")]) ## convert . to -
as_tibble(features)
features<-features%>%mutate(target=as.factor(target))
as_tibble(features)
f<-featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
f$layout
f
png(file=paste0("isT2DM_topvariables.png"),
                res=300, 
                width=2000, 
                height=1500)
featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
dev.off()


################################################################################
setwd("..")
## When you are done:
stopCluster(cl)




annot
rownames(top)
annot_preproc<-predict(prrproc,annot)
features<-data.frame(annot_preproc[,c("HSPB2","SMAD5-AS1","CLMAT3","REG3G","target")]) ## convert . to - ; HSPB2 SMAD5-AS1 CLMAT3 REG3G
as_tibble(features)
features<-features%>%mutate(target=as.factor(target))
as_tibble(features)
f<-featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
f$layout
f
png(file=paste0("topvariables.png"),
                res=300, 
                width=2000, 
                height=1500)
featurePlot(x = features[,1:4],
            y = features$target,plot = "pairs",
            auto.key = list(columns = 5))
dev.off()

################################################################################ pca_plots.R

annot
library(vegan)
annot_mat<-as.matrix(annot%>%select(where(is.numeric)))
rownames(annot_mat)<-annot$SampleID
dist<-vegdist(annot_mat,method = "euclidean")
distmat<-as.matrix(dist)
distmat

hdist<-hclust(dist,method = "ward.D2")
hdist
# ordiplot()

##################
## Draw a plot for a non-vegan ordination (cmdscale).
data(dune)

dune.dis <- dist
dune.mds <- cmdscale(dune.dis, eig = TRUE)
dune.mds$species <- wascores(dune.mds$points, annot_mat, expand = TRUE)
dune.mds$points

annot3<-annot%>%
  mutate(PCT2DM=dune.mds$points[,1],
         PCT2DM=dune.mds$points[,2])
annot3
annot3_slim<-annot3%>%select(SampleID,target,project,PCT2DM,PCT2DM)%>%distinct_all()%>%mutate(target2=target=="UNK")
keep<-annot3_slim%>%filter(target2=="TRUE")%>%mutate(target2=SampleID)
other<-annot3_slim%>%filter(target2!="TRUE")%>%mutate(target2=NA)
annot4_slim<-full_join(keep,other)

# my_pal<-pal_d3(palette = "cat")
ggscatter(annot4_slim,x = "PCT2DM",
          y = "PCT2DM",
          # label = "target2",
          color = "target",
          palette = "startrek",
          rug = F,
          shape = "project",repel = T,
          ggtheme = theme_pubr(legend = "right"),
          group="target")#+
  # geom_label(mapping = aes(label=target2),position = "jitter",)#+
  
  # geom_smooth(method = "glm",mapping = aes(group=target,fill=NA,color=target))+scale_color_aaas()
          # ellipse.alpha = 0,ellipse.level = 0.95,
          # ellipse.type = "norm",
          # ellipse = T)

################################################################################ ################################################################################ ################################################################################
################################################################################ ################################################################################ ################################################################################
################################################################################ ################################################################################ ################################################################################
reannotate 7 samples with predictions 
CIs for each test
unsupervised clustering, clustering genes 
visualizing top markers
overlap of top markers and DE analysis


