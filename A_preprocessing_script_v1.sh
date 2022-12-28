## preprocessing_script_v1.sh

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


## count_script_v1.sh
############################################################################################################
####################################  Use STAR to align to human reference (refdata-gex-GRCh38-2020-A)  #############################################
############################################################################################################
cd /media/jochum00/Aagaard_Raid2/ebarrozo/'New directory'/Barrozo/placenta


############################################################################################################
####################################  Use htseq2 to count transcripts  #############################################
############################################################################################################
## Make counts files for biological replicates
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