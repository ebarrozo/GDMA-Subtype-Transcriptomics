# GDMA-Subtype-Transcriptomics

Study design. This case-cohort study was designed to rigorously examine placental gene expression differences between gestational and pre-gestational diabetics. To accomplish this, we used a two phase study in which participants were recruited from labor and delivery at Ben Taub General Hospital in Houston, Texas is accordance with procedures and protocols approved by the Institutional Review Board of Baylor College of Medicine.  (IRB number H-28623, approved 04-26-2011). An initial discovery cohort of 65 subject’s samples were collected, consisting of 40 control, 6 GDMA1, 11 GDMA2, 1 Type 1 DM, and 7 Type 2 DM placentae, and analyzed by whole transcriptome sequencing. The Type 1 DM was excluded from the current analysis. The second phase was designed for independent validation of first phase results and included 145 placentas from non-overlapping subjects (29 control, 41 GDMA1, 43 GDMA2, and 32 T2DM samples). The cohorts were distinct subjects and samples, without individual subject overlap (Fig. S1). Extensive clinical metadata was collected from the electronic medical record and curated in a secured database by trained research personnel. Selected cases were routinely audited and adjudicated by a board-certified maternal-fetal medicine physician scientist to ensure diagnostic accuracy. 

Diagnostic criteria for cohort designation. Diagnosis of gestational diabetes was based on uniform established institutional criteria using the Carpenter-Coustan glucose tolerance test (GTT), with a screening glucose challenge test (GCT) value of greater than or equal to 140 mg/dl defining a positive screen. Pre-gestational diabetes was classified by one or more of the following: established diagnosis of diabetes prior to pregnancy, positive GTT, or elevated hemoglobin A1c early in pregnancy (HbA1c > 6.4). Hypertensive disorders of pregnancy, including both gestational hypertension and preeclampsia, were defined using the American College of Obstetricians and Gynecologists’ classification.61 

Sample collection and processing. All samples were collected by personnel trained in perinatal and placental pathology under strict uniform protocol. Briefly, following standard obstetrical practice the placenta was delivered and immediately passed to trained personnel in a sterile clean container. In a separate room, two samples were collected from midway between the cord insertion and placental margin by incision through the fetal surface into the parenchyma, but not to the maternal surface. All samples were collected within one hour of delivery under clean and sterile conditions as detailed above, placed on dry ice in sterile closed vials, transported to the laboratory, and stored at −80°C until mRNA extraction.

Transcriptomics. RNA was extracted from placental tissue as previously described.37 The Machery Nagel Nucleospin II kit was used to extract RNA from each sample. Samples were stored at -80°C. Each sample was analyzed for quality control on the Agilent Bioanalyzer, with a minimum RNA Integrity Number (RIN) of 4.0 accepted for transcriptomic analysis. For each library, mRNA was purified from 10µg of total RNA using the DynaBeads mRNA Purification Kit (Invitrogen) and fragmented using the RNA Fragmentation Reagents (Ambion). Double-stranded cDNA was synthesized from fragmented mRNA using the Superscript Double-Stranded cDNA Synthesis Kit (Invitrogen) and Random Hexamer primers (50ng/ul, Invitrogen). DNA sequencing libraries were generated from the cDNA according to manufacturer’s protocol and cluster generation and sequencing was performed on the Illumina cBot station and Illumina Hiseq 2000.
Sequenced reads were mapped to the human reference genome (GRCh38) with HISAT262 and assigned to genomic features with StringTie.63 Differential expression was determined using tximport and DESeq2,64,65 and genes with a false discovery rate (FDR) multiple testing corrected value <0.05 were considered differentially expressed. Factors that were examined for an association with differential expression were birth weight, BMI at delivery, chronic hypertension, delivery weight, diabetic state, diabetes medication, ethnicity, gestational age at birth, gravidity, height, maternal age, parity, preeclampsia or gestational hypertension, pregnant weight change, and pre-pregnancy/1st trimester weight. Consensus classes were assigned using ConsensusClusterPlus66 and weighted gene co-expression network analysis was performed using WGCNA.67 Gene set enrichment and ontology analysis was performed with EGSEA and clusterProfiler.68–70 Cytoscape software was used for network visualization.71

Sequencing quality validation and verification. From 64 samples, 4.4 billion reads were generated, of which 83% mapped to the reference genome. Estimation of experimental power of our dataset indicates a >80% power to detect genes with greater than two-fold change in expression at p<0.01 (Fig. S3a).28 After quantile normalization, all samples exhibited consistent expression and distribution as measured by gene abundance and density (Figures S3b, S3c) respectively. The accuracy of our analysis pipeline was verified and differential expression was measured between placental samples from males and females (Fig. S4a). Of 53 differentially expressed genes detected in this comparison, 36 were located on the Y chromosome and displayed higher expression in the male samples (Fig. S4b). Eight of the genes were located on the X chromosome, of which six were more highly expressed in females, including the X Inactive Specific Transcript (XIST) gene.

GEO Contents: The bulk RNA-seq data have been deposited to the Gene Expression Omnibus (GEO accession GSEXXXXXX).

# Repository Contents

[The reference GRCh38.p13 was downloaded from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build.] 

(A) Bash script for pre-processing and quality filtering reads, trimming adapters/barcodes, alignment using STAR, deduplication, and counts matrix generation.

(B) R script used for differential expression analysis.


# Workflow
Reads were pre-processed, and quality filtered using FastqQC (v0.11.9). 

Barcodes and adapters were trimmed using Trimmomatic (v0.33). 

Reads were aligned to the human  transcriptome (GRCh38.p13) using STAR (v2.7.8). 

PCR duplicate reads were removed using Picard (v2.24.0). 

Unique reads were counted using HTseq (v0.11.1). 

Counts were used for differential expression analysis in R (v4.0.2) using DEseq2 (v3.12). 
