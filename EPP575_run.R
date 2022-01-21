# Set-Up for Analysis
# Change the destination directory to match your own
setwd("~/Downloads/EPP575_testrun/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")
library(DESeq2)

# Get File Names
fileName <- list.files(path = ".", pattern = "*.counts.txt")
fileName

# Get Sample Names
sampleName <- sapply(strsplit(fileName, "[.]"), "[[", 1)
sampleName

# Get Mutant Level
mut <- factor(c("max2", "max2", "control", "control"))
mut

# Get Replicate Level
rep <- factor(c(1, 2, 1, 2))
rep

# Build the Sample Table
MaxSampleTable <- data.frame(sampleName, fileName, mut, rep)
MaxSampleTable

# DESeq Analysis Mutant vs Control
# Build DESeq Dataset
DESeqData_Max <- DESeqDataSetFromHTSeqCount(sampleTable = MaxSampleTable, directory = ".", design = ~ rep + mut)
DESeqData_Max
head(counts(DESeqData_Max))

# DESeq2 Analysis
DESeqData_Max.Res <- DESeq(DESeqData_Max)

# PCA Plots - Visualizes how related samples are to one another
rld <- rlog(DESeqData_Max, blind = FALSE)
plotPCA(rld, intgroup = c("rep"))

# Mutant Comparisons
MaxRes <- results(DESeqData_Max.Res, alpha = 0.05, contrast = c("mut", "max2", "control"))
summary(MaxRes)

# MA Plot - Visualizes log fold-change vs mean expression between treatments
plotMA(MaxRes, ylim=c(-6,6))

# Get Significant Genes
MaxRes_sig <- as.data.frame(MaxRes[ which(MaxRes$padj < 0.05),])
MaxRes_sig <- MaxRes_sig[order(MaxRes_sig$padj),]
head(MaxRes_sig, n=10)

# Optional 1 - Filter out genes with no expression
dim(DESeqData_Max)
DESeqData_Max <- DESeqData_Max[ rowSums(counts(DESeqData_Max)) > 1, ]
dim(DESeqData_Max)

# Then go back to rerun the DESeq2 Analysis step

# Optional 2 - Limit LFC
MaxRes_FC <- as.data.frame(MaxRes[ which(abs(MaxRes$log2FoldChange) > 1),])
# You can further subset this to include only significant genes with absolute log fold change >1
MaxRes_FC_sig <- as.data.frame(MaxRes_FC[ which(MaxRes_FC$padj < 0.05),])

# DESeq Analysis 2: Replicates 1 vs 2
# Build DESeq Dataset
DESeqData_Rep <- DESeqDataSetFromHTSeqCount(sampleTable = MaxSampleTable, directory = ".", design = ~ mut + rep)
DESeqData_Rep
head(counts(DESeqData_Rep))

# DESeq2 Analysis
DESeqData_Rep.Res <- DESeq(DESeqData_Rep)

# Replicate Comparisons
RepRes <- results(DESeqData_Rep.Res, alpha = 0.05, contrast = c("rep", 1, 2))
summary(RepRes)

# PCA Plot - What is the more important factor?
rld2 <- rlog(DESeqData_Rep, blind = FALSE)
plotPCA(rld2, intgroup = c("rep"))
