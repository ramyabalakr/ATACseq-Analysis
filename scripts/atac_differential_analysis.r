#Author: Ramya Balakrishnan

#Load libraries
library(GenomicAlignments)
library(GenomicRanges)
library("DESeq2")
library(edgeR)
library(BSgenome.Mmusculus.UCSC.mm9)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggfortify)
library(ggrepel)
library("pheatmap")
library(ggpubr)
library("gridExtra")
library("RUVSeq")
library("rtracklayer")
library(reshape2)

# Set working directory
setwd("C:/Users/balak/Desktop/ATACseq_Analysis")

#Load file describing CnT samples

atac_samples <- as.data.frame(read.table("ATAC_C2C12_A485.samples.merged.txt", sep="\t", header=TRUE))

# Reorder samples so order is DMSO, 3h, 6h. 
# Note: Blais_116 should be excluded from analysis
atac_samples <- atac_samples[c(1,5,9,13,14,15,18,21,24,2,6,10,3,7,11,4,8,12,16,19,22,17,23),]

#Sample metadata 

atac_samples$batch <- c(rep("July_2020", times=15),rep("March_2021",times=8))
atac_samples$ShortName <- paste(atac_samples$drug_name, atac_samples$drug_dose_uM, atac_samples$duration.hr, atac_samples$Replicate, sep="|")
atac_samples$DrugDose <- paste(atac_samples$drug_name, atac_samples$drug_dose_uM, sep="_")
atac_samples$DrugDoseDuration <- paste(atac_samples$DrugDose, atac_samples$duration.hr, sep="_")
atac_samples$Treatment_uM_duration <- paste(atac_samples$Treatment_uM, atac_samples$duration.hr, sep="_")
atac_samples$ShortName <- paste(atac_samples$Treatment, atac_samples$Replicate_number, atac_samples$Replicate_date, atac_samples$Replicate, sep="_")

#Load feature counts data; containing the number of sequencing reads that overlap with ATAC-seq peaks

featureCounts_tables <- read.table("A-485_ATAC_counts.read2pos5_bothEnds_noDups.Six1.peaks.txt")

#Again, reorder samples so order is DMSO, 3h, 6h. 
# Note: Blais_116 should be excluded from analysis
featureCounts_tables <- featureCounts_tables[,c(1,5,9,13,14,15,18,21,24,2,6,10,3,7,11,4,8,12,16,19,22,17,23)]

rownames(atac_samples) <- make.unique(atac_samples$ShortName)
colnames(featureCounts_tables) <- atac_samples$SampleName

# Transform featureCounts into an edgeR object.
# readDGE is a function from the edgeR package.
# The 'align_counts' edgeR object stores all the information regarding our experiment

align_counts <- DGEList(featureCounts_tables)

# Confirm that original sample information matches what is reported inside the edgeR object
# This should report "TRUE"
all(rownames(align_counts$samples) == rownames(atac_samples))

# If TRUE, we can add extra information about the samples to the edgeR object
align_counts$samples <- cbind(align_counts$samples, atac_samples)

#Confirm information was entered correctly. 
#Verify that the "$samples" section of the edgeR object is in the same order as the columns of the "$counts" section of the object.

rownames(align_counts$samples) <- colnames(align_counts$counts)
all(rownames(align_counts$samples) == colnames(align_counts$counts))

#Filter out low-quality peaks with counts < 20

idx_good_expression <- rowSums(align_counts$counts > 20) >= 0
table(idx_good_expression)
align_counts_filter <- align_counts[idx_good_expression, ]
	
#Print out the number of peaks that were retained and removed 
cat("Initial peaks:", nrow(align_counts$counts), " | Retained peaks:", nrow(align_counts_filter$counts), " | Removed peaks:", nrow(align_counts$counts) - nrow(align_counts_filter$counts), "\n")

#Generate RLE and PCA Plots to Evaluate Optimal Normalization Methods

#Create R function for PCA plot creation 
								plotPCA1234 <- function (object, pc_to_plot=c(1,2), labels = FALSE, isLog = FALSE, ...)
								{
								if (!isLog) {
								  Y <- apply(log(object + 1), 1, function(y) scale(y, center = TRUE, scale = FALSE))
								}
								else {
								  Y <- apply(object, 1, function(y) scale(y, center = TRUE, scale = FALSE))
								}
								s <- svd(Y)
								percent <- s$d^2/sum(s$d^2) * 100
								labs <- sapply(seq_along(percent), function(i) {
								  paste("PC ", i, " (", round(percent[i], 2), "%)", sep = "")
								})
								if (labels) {
								  plot(s$u[, pc_to_plot[1]], s$u[, pc_to_plot[2]], type = "n", ..., xlab = labs[pc_to_plot[1]], ylab = labs[pc_to_plot[2]])
								  text(s$u[, pc_to_plot[1]], s$u[, pc_to_plot[2]], labels = colnames(object), ...)
								}
								else {
								  plot(s$u[, pc_to_plot[1]], s$u[, pc_to_plot[2]], ..., xlab = labs[pc_to_plot[1]], ylab = labs[pc_to_plot[2]])
								}
								}


# Define the design matrix to include the effect of DrugDose
design_a = model.matrix(~0 + DrugDose, data = droplevels(atac_samples))

# Filter the count data to only include peaks with good expression
y.TMM <- align_counts[idx_good_expression, ]

# Calculate normalization factors using the TMM method
y.TMM <- calcNormFactors(y.TMM, method="TMM")

# Set the design matrix for normalization and dispersion estimation
design_for_normtest <- design_a
		
# Ensure that the row names of atac_samples match the sample names		
rownames(atac_samples) <- atac_samples$SampleName

# Estimate common dispersion for the GLM model using the specified design
y_a.TMM <- estimateGLMCommonDisp(y.TMM, design_for_normtest)
# Estimate tag-wise dispersion for the GLM model
y_a.TMM <- estimateGLMTagwiseDisp(y_a.TMM, design_for_normtest)
# Fit the GLM model to the data using the specified design
fit_y_a.TMM <- glmQLFit(y_a.TMM, design_for_normtest)
# Calculate counts per million (CPM) for the normalized data
y_a_cpm.TMM <- cpm(y_a.TMM, normalized.lib.sizes = TRUE, log=FALSE)
# Create a new SeqExpressionSet object with the normalized counts and sample information
y_a_set.TMM <- newSeqExpressionSet(counts=y_a.TMM$counts, phenoData = atac_samples)
# Scale the counts for visualization purposes
multiplier <- 1e-6*(mean(y_a.TMM$samples$lib.size))
normCounts(y_a_set.TMM) <- round(multiplier*y_a_cpm.TMM, digits=0)
# Calculate residuals from the GLM fit
res_y_a.TMM <- residuals(fit_y_a.TMM, type="deviance")

# Define the range of RUVr k values
for (k in 1:3) {
  assign(paste0("y_a_set.TMM_ruvrk", k), 
         RUVr(y_a_set.TMM, rownames(y_a_set.TMM), k = k, res_y_a.TMM))
}

# Open a PDF file for saving the plots

pdf("PCA_plots.test_designs+ruvr.pdf", height = 26, width = 30,  title="edgeR different normalization methods")
# Set the layout for the plots
# First number is the number of rows, second is the number of columns			
par(mfrow=c(5,6), oma = c(0, 0, 8, 0))
# Define colors and shapes for the plots
treatment_cols <- c("grey","darkcyan","blue3","deepskyblue","dimgrey","pink","red")
PCA_cols <- c(rep(treatment_cols[1], times=6),rep(treatment_cols[2], times=3),rep(treatment_cols[3], times=3),rep(treatment_cols[4], times=3),rep(treatment_cols[5], times=3),rep(treatment_cols[6], times=3))
PCA_shapes <- c(16,17,18,25,5,3,16,17,18,16,17,18,16,17,18,16,17,18,16,17,18)

# Define legends as strings to be evaluated later
rle_legend <- 'legend("topright", legend=c("DMSO_24h", "A485 0.3", "A485 1.0", "A485 3.0"), col=treatment_cols, pch=c(15,15,15), cex=1.0)'
pca_legend <- 'legend("topleft", legend=c("DMSO_24h", "A485 0.3", "A485 1.0", "A485 3.0", "A485 1.0 3h", "A485 1.0 6h", "Rep_1", "Rep_2", "Rep_3","Rep_4","Rep_5","Rep_6"), col=c(treatment_cols, "black", "black", "black"), pch=c(15,15,15, 15, 15, 15, 15, 16, 17, 18, 25, 5, 3), cex=1.0)'


# Define parameters
pc_combinations <- list(c(1,2), c(2,3))
ruvrk_values <- c(1, 2, 3)

# Loop through combinations and ruvrk values
for (k in ruvrk_values) {
  for (pc in pc_combinations) {
    plotPCA1234(normCounts(get(paste0("y_a_set.TMM_ruvrk", k))), 
                col=PCA_cols, cex=2, 
                main=paste("TMM, design a ruvrk", k, ", PCA", pc[1], "+", pc[2]), 
                isLog=FALSE, pc_to_plot=pc, pch=PCA_shapes, labels=FALSE)
  }
}

# Close the PDF file
dev.off()			

# Convert atac_samples to a data frame for easier manipulation
metaData <- data.frame(atac_samples)

# Get RUVr-adjusted expression values for all genes using k=1
table_log2_ruvrk1 <- log2(normCounts(y_a_set.TMM_ruvrk1))

# Median center the expression estimates to reduce bias
table_log2_ruvrk1_medCenter <- table_log2_ruvrk1 - apply(table_log2_ruvrk1, 1, median)
m <- table_log2_ruvrk1_medCenter

# Load cluster information from a saved R data file
load(file='cluster_info_k=3.Rda')
p300_cluster_info <- cluster_info

# Define gap positions for clustering based on cluster sizes
cluster_gaps <- as.integer(list(982,3418))

# Reorder the peak information to match the p300 cluster information
new_df <- m[ order(row.names(p300_cluster_info)), ]

# Combine the expression data with cluster information and reorder based on cluster ranks
m2 <- cbind(m, p300_cluster_info$k.means.cluster)
o <- order(m2[, ncol(m2)])
m3 <- m[o, ]

# Set the title for the heatmap plot
plot_title <- "ATAC-seq signal at Six1 peaks\nsiCtrl (Ctrl) vs siSix1 conditions\nRemoval of batch identified by RUVr k=1"

# Create a color palette for the heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("deepskyblue2", "black", "yellow"))(paletteLength)
breaksmin <- -1
breaksmax <- 1
myBreaks <- c(seq(from=breaksmin, to=breaksmax, by=(breaksmax-breaksmin)/paletteLength))
myColorAdjusted <- c(myColor[1], myColor, myColor[paletteLength])
myBreaksAdjusted <- c(-2, myBreaks, 2)

# Prepare data for the heatmap
annot.df <- atac_samples[,c(7,6)]
m4 <- m3[,c(1,2,3,4,5,6,7,8,9,10,11,12,19,20,21,22,23,13,14,15,16,17,18)]

# Create a heatmap PDF
pdf(file = "colour-adjusted_A-485_ATAC-seq_Heatmap.pdf", width = 8, height= 10)
b1 <- pheatmap(m4, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE, 
                annotation_row=p300_cluster_info, annotation_col=annot.df, 
                color = myColorAdjusted, breaks=myBreaksAdjusted, 
                main=plot_title, width=10, height=7, border_color=NA, gaps_row=cluster_gaps)
dev.off()

# Check if row names of m4 match those of cluster_info
all(rownames(m4) == rownames(cluster_info[rownames(m4), ]))

# Reorder m4 to match cluster info
m4 <- m4[rownames(cluster_info),]

# Combine m4 with cluster info
m4_CI <- cbind(m4, cluster_info)

# Melt the combined data frame for plotting
m4_CI_melt <- reshape2::melt(m4_CI)

# Create a boxplot for cluster comparisons
cluster_boxplot <- ggplot(m4_CI_melt, aes(x="", y=value, group=variable, fill=variable)) +
  geom_boxplot(notch=FALSE, outlier.shape=NA) +
  facet_grid(. ~ k.means.cluster) + theme_bw() +
  scale_fill_manual(values=c("#785EF0", "#785EF0","#785EF0", "#785EF0","#785EF0", "#785EF0","#785EF0", "#785EF0","#785EF0","#1C2371","#1C2371","#1C2371","#00FF1B","#00FF1B","#00FF1B", "#E4FF24", "#E4FF24","#DC267F","#DC267F","#DC267F","#FFB000","#FFB000","#FFB000")) +
  labs(title=NULL, x ="Drug treatment", y = "Median-centered ATAC-seq signal") +
  scale_y_continuous(limits = c(-2, 2))

# Save the boxplot
ggsave("A-485_ATAC-seq.boxplot.pdf", cluster_boxplot, width=8, height=11)