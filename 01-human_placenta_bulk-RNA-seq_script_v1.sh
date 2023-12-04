## human_placenta_bulk-RNA-seq_script_v1.sh

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## Analyze data on AagaardLab3



# make gdm.placenta directory and tree
cd
mkdir gdm.placenta
cd gdm.placenta
mkdir scripts
mkdir data
mkdir docs
mkdir results
cd results
mkdir qc
cd ..

# bulk RNA-seq from human placenta tissue: Ctrl, GDMA2, Type II
## metadata located in 'RNA seq_Clinical Validation_13 May 2016_v2 (2).xlsx' /Users/enricobarrozo/Library/CloudStorage/Box-Box/Melissa-GDM/RNA seq_Clinical Validation_13 May 2016_v2 (2).xlsx
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta
ls
	## DM1 thru DM11, PE, stored as .bz2 e.g. DM1.read1.bz2 and DM1.read2.bz2 ; 22 total
	## NCS5 thru 88, PE, stored as .bz2 e.g. NCS5.read1.bz2 and NCS5.read2.bz2 ; 54 total ; maybe ChIP? based on NCS designations in metadata excel sheet
	## 6002 thru 6171, PE, stored as .bz2 e.g. 6002.read1.bz2 and 6002.read2.bz2 ; 51 total ; 6005 only has read1


## run fastqc to determine quality, read length, and adapter/index retention
fastqc DM1.read1.bz2 DM1.read2.bz2 --threads 32 --outdir /home/ebarrozo/gdm.placenta/results/qc
	## 100 bp, no adapters detected
		## start generate genome with sjdb 99 and trimmomatic

############################################################################################################
####################################  Use the human reference from 10x genomics refdata-gex-GRCh38-2020-A  #############################################
############################################################################################################

cd /home/ebarrozo/gdm.placenta/docs
## Generate STAR genome index; make sure sjdbOverhang matches read length minus 1, as determined by fastqc outputs b4 and after trimming
mkdir STAR
cd STAR
STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/gdm.placenta/docs/STAR \
--genomeFastaFiles /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--sjdbGTFfile /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--sjdbOverhang 99

## see STAR manual 2.7.0a https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

############################################################################################################
####################################  Convert .bz2 to .gz and run trimmomatic  #############################################
############################################################################################################
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta
parallel --dryrun 'bzcat {} | gzip -c > {.}.fastq.gz' ::: *bz2

parallel 'bzcat {} | gzip -c > {.}.fastq.gz' ::: *bz2

## Run Trimmomatic PE loop
for file in *.read1.fastq.gz; do
withpath="${file}"
filename=${withpath##*/}
base="${filename%.read1.fastq.gz}"
echo "${base}"
trimmomatic PE -threads 40 "${base}".read1.fastq.gz "${base}".read2.fastq.gz "${base}"_1.trimmed_PE.fastq.gz "${base}"_1.trimmed_SE.fastq.gz "${base}"_2.trimmed_PE.fastq.gz "${base}"_2.trimmed_SE.fastq.gz SLIDINGWINDOW:4:15 MINLEN:50 ILLUMINACLIP:/media/jochum00/Aagaard_Raid2/ebarrozo/adapters.fa:2:30:10
done

## see if adding the adapters found in the first fastqc run fixed the problem
fastqc JDM1.trimmed_PE.fastq.gz --threads 32 --outdir /home/ebarrozo/gdm.placenta/results/qc

## Remove the files that won't be used
for f in *SE.fastq.gz; do
rm $f 
done
for f in *2.trimmed_PE.fastq.gz; do
rm $f 
done

## Remove the .gz files converted from the original .bz2 files after trimming is complete
for f in *.read1.fastq.gz; do
rm $f 
done
for f in *.read2.fastq.gz; do
rm $f 
done

############################################################################################################
####################################  Use STAR to align to human reference (refdata-gex-GRCh38-2020-A)  #############################################
############################################################################################################
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta

## Align to GCA_003118495.1 using STAR with the standard ENCODE options (see pg 8 in the STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
for f in *_1.trimmed_PE.fastq.gz; do
STAR --runMode alignReads --runThreadN 40 \
--genomeDir /home/ebarrozo/gdm.placenta/docs/STAR \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--sjdbGTFfile /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
--sjdbOverhang 99 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn <(gunzip -c $f) \
--outFileNamePrefix ${f%.trimmed.fastq}
done

# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file JMC11_HFRes_K4_trimmed.fastq.gz_STARtmp//BAMsort/19/46
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.

# Convert logs to txt files and save to local
for f in *Log.final.out; do
cat $f > ${f%.out}.txt
done


############################################################################################################
####################################  Remove PCR duplicates using picard rmdup  #############################################
############################################################################################################
## Remove PCR duplicates
for f in *Aligned.sortedByCoord.out.bam; do
picard MarkDuplicates \
INPUT=$f \
OUTPUT=${f%Aligned.sortedByCoord.out.bam}.duprm.bam \
METRICS_FILE=${f%Aligned.sortedByCoord.out.bam}.met \
REMOVE_DUPLICATES=true
done

############################################################################################################
####################################  Use htseq2 to count transcripts  #############################################
############################################################################################################
## Make counts files for biological replicates
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta

for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done

## htseq requires gene_id in the each row of the .gtf file
cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/
wc -l genes.gtf
	# 2765974 genes total
# grep "exon" genes.gtf > filtered.genes.gtf
# wc -l filtered.genes.gtf
	## 1025984
# grep "CDS" genes.gtf > filtered.genes.gtf
# wc -l filtered.genes.gtf
	## 867547
grep "gene_id" genes.gtf > filtered.genes.gtf
wc -l filtered.genes.gtf	
	## 2765969 lines with gene_id

cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta

## removed -n which sorts by name not position
for f in *.duprm.bam; do
samtools sort -@ 48 $f -o ${f%.bam}.sorted.bam
done

for f in *duprm.sorted.bam; do
samtools index $f -@ 48
done

## added -r pos to specify alignmed sorted by position
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no -r pos --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/filtered.genes.gtf > ${f%.bam}.htseq.counts.txt
done


for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/gdm.placenta/results
done



################### for IGV

cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta
for f in *.sorted.bam; do
bedtools bamtobed -i $f > ${f%.bam}.bed
done

#cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta
#for f in *.sorted.bam; do
#	samtools index $f 
#done

# mkdir /home/ebarrozo/gdm.placenta/results/IGV
for f in *.sorted.bed; do
cp $f /home/ebarrozo/gdm.placenta/results/IGV
done

for f in *.bai; do
cp $f /home/ebarrozo/gdm.placenta/results/IGV
done

###################################################################################################################################### 
######################################	Remove files to save space	#############################################################
######################################################################################################################################
for f in *.gz.duprm.bam; do
	rm $f
done

for f in *Aligned.sortedByCoord.out.bam; do
	rm $f
done

for f in *.fastq.gz; do
	rm $f
done