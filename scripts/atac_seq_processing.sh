#!/bin/bash -i
# version 1.05

# Set error handling
set -e

# Track execution time
res1=$(date +%s.%N)

### ATAC-seq analysis script
### Suitable for paired-end data

# Define input file 
i=$1

# Define project directories
home_projects=$HOME/projects/def-ablai2
project_dir=$HOME/scratch/ATAC_OCT
final_project_dir=$home_projects/ATAC-seq-ramya
process_report_file=$final_project_dir/ATACseq_01.process_report.$i.out
temp_directory=$HOME/scratch/
genome_dir=$home_projects/Genomes/mm9/

# Load required modules
module load nixpkgs/16.09 intel/2018.3
module load fastp/0.20.0 star/2.7.0a picard/2.18.9 samtools/1.9 fastqc/0.11.9
module load gcc/7.3.0 gsl/2.5
module load r/3.6.0

# Create directories if they don't exist
mkdir -p $project_dir/fastQC_reports
mkdir -p $project_dir/fastp_reports
mkdir -p $project_dir/fastq_files
mkdir -p $project_dir/bam_and_bai_files
mkdir -p $project_dir/STAR_reports
mkdir -p $project_dir/picard_reports
mkdir -p $project_dir/ataqv_reports
mkdir -p $project_dir/bedgraph_bigwig_files
mkdir -p $project_dir/fingerprinting_reports
mkdir -p $project_dir/peak_calling_MACS2

echo "Starting computation for $i" >> $process_report_file

# Load required modules
module load nixpkgs/16.09 intel/2018.3
module load fastp/0.20.0 star/2.7.0a picard/2.18.9 samtools/1.9 fastqc/0.11.9

#Concatenate the reads coming from the two lanes into one file 
echo Concatenating input files for $i >>$process_report_file

cat $ori_fastq/$i"_L001_R1_001.fastq.gz" \
$ori_fastq/$i"_L002_R1_001.fastq.gz" \
>$project_dir/fastq_files/$i"_L0012_R1_001.fastq.gz"

cat $ori_fastq/$i"_L001_R2_001.fastq.gz" \
$ori_fastq/$i"_L002_R2_001.fastq.gz" \
>$project_dir/fastq_files/$i"_L0012_R2_001.fastq.gz"

# FastQC prior to trimming
echo "FastQC before trimming for $i" >> $process_report_file

# FastQC Read1
zcat $project_dir/fastq_files/$i"_L0012_R1_001.fastq.gz" |\
fastqc stdin \
-o $project_dir/fastQC_reports/ \
-f fastq \
--contaminants $home_projects/Genomes/fastQC_sequences/contaminant_list.txt \
--adapters $home_projects/Genomes/fastQC_sequences/adapter_list.txt \
--limits $home_projects/Genomes/fastQC_sequences/limits.txt \
-d $temp_directory \
>>$process_report_file 2>&1

# FastQC Read2
zcat $project_dir/fastq_files/$i"_L0012_R2_001.fastq.gz" |\
fastqc stdin \
-o $project_dir/fastQC_reports/ \
-f fastq \
--contaminants $home_projects/Genomes/fastQC_sequences/contaminant_list.txt \
--adapters $home_projects/Genomes/fastQC_sequences/adapter_list.txt \
--limits $home_projects/Genomes/fastQC_sequences/limits.txt \
-d $temp_directory \
>>$process_report_file 2>&1

### Remove low quality sections of reads and remove sequencing adapter contamination
# Note: when a command is very long, we can break it over multiple lines as this helps with readability of the code
echo "Trimming reads with fastp for $i" >> $process_report_file

fastp \
--in1 $project_dir/fastq_files/$i"_L0012_R1_001.fastq.gz" \
--in2 $project_dir/fastq_files/$i"_L0012_R2_001.fastq.gz" \
--thread 16 \
--length_required 25 \
--adapter_sequence CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
--adapter_sequence_r2 CTGTCTCTTATACACATCTGACGCTGCCGACGA \
--cut_tail \
--cut_tail_window_size 4 \
--cut_tail_mean_quality 20 \
--disable_quality_filtering \
--overrepresentation_analysis \
--overrepresentation_sampling 10 \
--json $project_dir/fastp_reports/$i.fastp_report.json \
--html $project_dir/fastp_reports/$i.fastp_report.html \
--out1 $project_dir/fastq_files/$i.trimmed.R1.fastq \
--out2 $project_dir/fastq_files/$i.trimmed.R2.fastq \
>>$process_report_file 2>&1

# FastQC after trimming
echo "FastQC after trimming for $i" >> $process_report_file

fastqc \
-o $project_dir/fastQC_reports/ \
-f fastq \
--contaminants $home_projects/Genomes/fastQC_sequences/contaminant_list.txt \
--adapters $home_projects/Genomes/fastQC_sequences/adapter_list.txt \
--limits $home_projects/Genomes/fastQC_sequences/limits.txt \
-d $temp_directory \
$project_dir/fastq_files/$i.trimmed.R1.fastq $project_dir/fastq_files/$i.trimmed.R2.fastq \
>>$process_report_file 2>&1

# Map reads with STAR
echo "Mapping reads with STAR for $i" >> $process_report_file

# Set the parameters so that the maximum INTRON size is 1 bp, and the maximum FRAGMENT size is 2000 bp
STAR \
--runThreadN 16 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--alignMatesGapMax 2000 \
--outFilterMatchNminOverLread 0.10 \
--outFilterScoreMinOverLread 0.10 \
--genomeDir $genome_dir/star_index_noAnnot/ \
--readFilesIn $project_dir/fastq_files/$i.trimmed.R1.fastq $project_dir/fastq_files/$i.trimmed.R2.fastq \
--outSAMtype BAM SortedByCoordinate \
--outSAMmapqUnique 40 \
--outFilterMultimapNmax 1 \
--outFileNamePrefix $project_dir/bam_and_bai_files/$i.star_index_noAnnot. \
>>$process_report_file 2>&1

# Move STAR report files
mv $project_dir/bam_and_bai_files/$i.star_index_noAnnot.Log.out $project_dir/STAR_reports/
mv $project_dir/bam_and_bai_files/$i.star_index_noAnnot.Log.final.out $project_dir/STAR_reports/
mv $project_dir/bam_and_bai_files/$i.star_index_noAnnot.Log.progress.out $project_dir/STAR_reports/

# Index BAM file
echo "Indexing after STAR for $i" >> $process_report_file

samtools index -@ 16 $project_dir/bam_and_bai_files/$i.star_index_noAnnot.Aligned.sortedByCoord.out.bam \
$project_dir/bam_and_bai_files/$i.star_index_noAnnot.Aligned.sortedByCoord.out.bam.bai \
>>$process_report_file 2>&1

# Mark duplicates with Picard
echo "Marking duplicates with Picard for $i" >> $process_report_file

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
VERBOSITY=WARNING \
I=$project_dir/bam_and_bai_files/$i.star_index_noAnnot.Aligned.sortedByCoord.out.bam \
O=$project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam \
M=$project_dir/picard_reports/$i.star_index_noAnnot.marked_dup_metrics.txt \
TMP_DIR=$temp_directory \
>>$process_report_file 2>&1

# Index BAM file after marking duplicates
echo "Indexing after Picard for $i" >> $process_report_file

samtools index -@ 16 $project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam \
$final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam.bai \
>>$process_report_file 2>&1

# Collect insert size metrics

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
I=$project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam \
O=$project_dir/picard_reports/$i.star_index_noAnnot.markedDups.insert_metrics_table.txt \
H=$project_dir/picard_reports/$i.star_index_noAnnot.markedDups.insert_metrics_histogram.pdf \
HISTOGRAM_WIDTH=2000 \
>>$process_report_file 2>&1

# Peak calling with MACS2
module purge

# Load necessary modules for MACS2
module load python/3.7.0 scipy-stack nixpkgs/16.09 gcc/7.3.0 intel/2018.3 samtools/1.9 bedtools/2.29.2

# Prepare BAM files for MACS2 in BED format
echo "Preparing BAM files for MACS2 in BED format for $i" >> $process_report_file

# Convert BAM to BED for first and second mates
# samtools with -f 64 restricts the output to only those reads (mates) that are FIRST IN PAIR, while RETAINing those reads marked as duplicates. 

samtools view -hf 64 - | samtools view -hF 1024 - $final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam |\
bedtools bamtobed -i stdin >$final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.first_noDups.bed

#samtools with -f 128 restricts the output to only those reads (mates) that are SECOND IN PAIR,RETAIN those reads marked as duplicates. 
samtools view -hf 128 - | samtools view -hF 1024 $final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam |\
bedtools bamtobed -i stdin >$final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.second_noDups.bed

# Convert BAM to BED for both mates
# To exclude supplementary alignments: -hF 1024 
samtools view -hF 1024 - $final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.bam |\
bedtools bamtobed -i stdin >$final_project_dir/bam_and_bai_files/$i.star_index_noAnnot.markedDups.both_noDups.bed

# Perform MACS2 peak calling
echo "MACS2 peak calling for $i" >> $process_report_file

samples='Blais_002 Blais_003 Blais_004 Blais_006 Blais_007 Blais_008 Blais_010 Blais_011 Blais_012'
# Define file name extension
term=star_index_noAnnot.markedDups.bam
# Define file name extension of files after sorting
termF=star_index_noAnnot.markedDups.first_noDups.bed
termS=star_index_noAnnot.markedDups.second_noDups.bed
termB=star_index_noAnnot.markedDups.both_noDups.bed
term_final=star_index_noAnnot.markedDups.both_noDups.sorted.bed

# Define location where bam files are stored
bam_loc=$HOME/projects/def-ablai2/ATAC-seq-ramya/bam_and_bai_files/

module load nixpkgs/16.09 intel/2018.3
module load samtools/1.9 bedtools/2.29.2

for i in $samples
do
# note that splitting reads like this with flags 64 and 128 will not remove reads marked as duplicates.
samtools view -hf 64 $bam_loc/$i.$term | samtools view -hF 1024 - | bedtools bamtobed -i stdin >$bam_loc/$i.$termF
samtools view -hf 128 $bam_loc/$i.$term | samtools view -hF 1024 - | bedtools bamtobed -i stdin >$bam_loc/$i.$termS
done

# MACS2 peak calling specifying input in bed format

module load python/3.7.0
module load scipy-stack
module load nixpkgs/16.09 gcc/7.3.0
source ~/projects/def-ablai2/env_MACS2/bin/activate

macs_out=$HOME/projects/def-ablai2/ATAC-seq-ramya/Peak_calling/MACS2

for i in $samples
do
macs2 callpeak \
-t $bam_loc/$i.$termF $bam_loc/$i.$termS \
-f BED \
-g mm \
--nomodel --shift -37 --extsize 73 --keep-dup all \
--outdir $macs_out \
--name $i.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73
done

#Remove blacklisted regions from MACS2 output narrowPeak files
module load bedtools/2.29.2

for i in $samples
do
bedtools intersect -v \
-a $macs_out/$i.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73_peaks.narrowPeak \
-b ~/projects/def-ablai2/Genomes/mm9/mm9-blacklist.bed \
>$macs_out/$i.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak

done

# Convert the narrowpeak file into a bed file that can be used by MSPC (improves peak calling when we have replicates)
# Rank by decreasing q value and place the q value in the 5th column

sort -k9,9nr \
$macs_out/$i.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak |\
awk -F "\t" '{printf($1"\t"$2"\t"$3"\t"$4"\t"$9"\n")}' - \
>$macs_out/$i.MACS2_bed_wDups.s37_e73.noBL.mspcBED

# Run MSPC to aggregate replicates
echo "MSPC peak calling for $i" >> $process_report_file
module load nixpkgs/16.09
module load dotnet-core/3.0.0

dotnet ~/projects/def-ablai2/mspc/mspc.dll \
-i $macs_out/Blais_002.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_003.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_004.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_006.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_007.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_008.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_010.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_011.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
$macs_out/Blais_012.star_index_noAnnot.markedDups.MACS2_bed_noDups.s37_e73.noBL.narrowPeak \
-r bio -w 1e-5 -s 1e-10 -d 6 -c 6 --parser ~/projects/def-ablai2/mspc/mspc_macs2_config_pvalue.json \
--output $macs_out/A485_new_MSPC_output_w-5_s-10_c6_parser_pval_noDups/

# Use the MACS2 output, filtered out of black-list peaks and processed with MSPC with parameters w-5_s-10_c6_parser_pval_noDups

#Rename output file
cp $macs_out/A485_new_MSPC_output_w-5_s-10_c6_parser_pval_noDups/ConsensusPeaks.bed \
$macs_out/A485_new_MSPC_output_w-5_s-10_c6_parser_pval_noDups/C2C12_GM_A485.ATAC_cons_peaks_noDups.mspcBED

# Use R to format the consensus file into a correct BED file
# Take the percentile rank of xSqrd values, with maximum of 1000

module load nixpkgs/16.09  intel/2018.3 gcc/7.3.0
module load r/3.6.0

Rscript -e '
  data <- read.table("~/projects/def-ablai2/ATAC-seq-ramya/Peak_calling/MACS2/A485_new_MSPC_output_w-5_s-10_c6_parser_pval_noDups/ConsensusPeaks_mspc_peaks.txt", sep="\t", header=TRUE)

  perc.rank <- function(x) round(1000*(trunc(rank(x))/length(x)), digits=0)
  data <- within(data, score <- perc.rank(xSqrd))
  data$strand <- rep(".", times=nrow(data))
  data2 <- data[,c(1:4,9,10)]
  write.table(data2, file="~/projects/def-ablai2/ATAC-seq-ramya/Peak_calling/MACS2/A485_new_MSPC_output_w-5_s-10_c6_parser_pval_noDups/C2C12_GM_A485.ATAC_cons_peaks_noDups.BED",  sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
'

### Keep track of time
# This records the time when the execution ended
res2=$(date +%s.%N)

# This calcualtes the difference in time between start and end
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

### Calculate time difference and report it at the end.
echo "#################################" >>$process_report_file
echo Subject: Finished $i >>$process_report_file
printf "Total runtime (days-hours-minutes-seconds-nanosec): %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds >>$process_report_file