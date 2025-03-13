# Start interactive session by logging on a compute node
# salloc --account=def-ablai2 --time=6:0:0 --x11 --mem-per-cpu=10Gb -n 10 -N 1

# load modules needed to run R
# module load nixpkgs/16.09  gcc/7.3.0 gsl/2.5
# module load r/3.6.0
# start R
R
# Load libraries
library(ATACseqQC)
library(GenomicAlignments)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm9)
library(magrittr)
library(dplyr)
library(Biostrings)
library(Rsubread)
library(rtracklayer)

# Set working directory
setwd("~/projects/def-ablai2/2021-03-31_ATAC-Ramya/3,6,24h_Analysis")

# Load file describing ATAC-seq samples
atac_samples <- as.data.frame(read.table("ATAC_C2C12_A485.samples.merged.txt", 
                                         sep = "\t", header = TRUE))

# Create sample names with A-485 dosage
atac_samples$SampleName <- paste(atac_samples$Sample, atac_samples$Treatment.uM., sep = "_")

# Define extra columns for narrowPeak file import
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", 
                          qValue = "numeric", peak = "integer")

# Path to Six1 CUT&Tag narrowPeak file
path_to_file_co8 <- "~/projects/def-ablai2/CnT-seq-Ramya/RB_AS_qValue-co8_Six1_ctrl_macs2_bampe_noDups.noBL.noIgG_peaks.narrowPeak"

# Import narrowPeak file
raw_peaks <- import.bed(path_to_file_co8, extraCols = extraCols_narrowPeak)

# Remove peaks in blacklist or on chrM
# Assuming blkList is defined elsewhere in your code
clean_peaks <- raw_peaks[!raw_peaks %over% blkList & !seqnames(raw_peaks) %in% "chrM"]

# Export cleaned peaks to BED file
export.bed(clean_peaks, 'CnT_Six1_co8_combined_bampe_nodups_c2.bed')

# Convert peaks to data frame for featureCounts
regionsToCount <- data.frame(
  GeneID = paste("ID", seqnames(clean_peaks), start(clean_peaks), end(clean_peaks), sep = "_"),
  Chr = seqnames(clean_peaks),
  Start = start(clean_peaks),
  End = end(clean_peaks),
  Strand = strand(clean_peaks)
)

# Specify A-485 Treatment ATAC-seq BAM files to count signal at Six1 CUT&Tag peaks
bamsToCount <- as.character(atac_samples$Path)

# Perform feature counting, ignoring duplicates
fcResults.noDups <- featureCounts(
  bamsToCount, 
  annot.ext = regionsToCount, 
  isPairedEnd = FALSE, 
  autosort = TRUE, 
  countMultiMappingReads = FALSE, 
  read2pos = 5, 
  ignoreDup = TRUE, 
  nthreads = 10
)

# Extract counts
myCounts_noDups <- fcResults.noDups$counts

# Rename columns with sample names
colnames(myCounts_noDups) <- as.character(atac_samples$SampleName)

# Write counts to file
write.table(myCounts_noDups, 
            file = "A-485_ATAC_counts.read2pos5_bothEnds_noDups.Six1.peaks.txt")
