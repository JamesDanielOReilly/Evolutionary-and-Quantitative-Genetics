###############################################
###############################################
##### Analyzing RNA-seq data with DESeq2 ######
###############################################
###############################################

# script partly based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html and 'Beginner's guide to using the DESeq2 package', Love et al. 2014
# data from:

# Swaegers, J., Spanier, K. I., & Stoks, R. (2020). Genetic compensation rather than 
# genetic assimilation drives the evolution of plasticity in response to mild warming 
# across latitudes in a damselfly. Molecular Ecology. 2020 doi: 10.1111/mec.15676



# In this script we explore the expression variation between the samples, 
# infer differentially expressed genes (DEGs) from several contrasts
# and visualise the DEG analysis

# ! Attention: in the paper sample replicates were collapsed for analysis purposes further on.
# Here we keep them as seperate samples. Differences in results hence occur.

# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("DESeq2")
BiocManager::install('EnhancedVolcano')
install.packages('pheatmap')

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# read count matrix and metadata
setwd("/home/james/Documents/leuven/second-year/Evolutionary-and-Quantitative-Genetics/part-two/sixth-assignment")
count_data <- read.table("count_data_exercise.txt", row.names = 1, header = T)
str(count_data)
design <- read.table("design.txt", row.names = 1, header = T)
design$treat <- factor(design$treat)
design$lat <- factor(design$lat)
group2<-factor(design$group2)

# run DESeq2

# this pipleline consists of the estimation of size factors 
# (which control for differences in the library size of the sequencing experiments),
# the estimation of dispersion for each gene (a measure of spread or variability in the data), 
# and fitting a generalized linear model.

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = design,
                              design = ~ group2)
dds <- DESeq(dds)

# Principal component analysis

# First apply a variance stabilizing transformation on dds
vsd=vst(dds)

# PCA
pcaDataDEG <- plotPCA(vsd,intgroup=c("treat","lat"),returnData=TRUE)
percentVar <- round(100 * attr(pcaDataDEG, "percentVar"))

ggplot(pcaDataDEG, aes(x = PC1, y = PC2, color = treat, shape = lat)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()


# Detection of DEGs per contrast

# N24 vs N20
res_N24_vsN20 <- results(dds, contrast = c("group2", "N24", "N20"))
sig_N24_vsN20 <- res_N24_vsN20[ which(res_N24_vsN20$padj < 0.05 ), ] 

# S24 vs S20
res_S24_vsS20 <- results(dds, contrast = c("group2", "S24", "S20"))
sig_S24_vsS20 <- res_S24_vsS20[ which(res_S24_vsS20$padj < 0.05 ), ] 

# S24 vs N24
res_S24_vsN24 <- results(dds, contrast = c("group2", "S24", "N24"))
sig_S24_vsN24 <- res_S24_vsN24[ which(res_S24_vsN24$padj < 0.05 ), ]

# S20 vs N20
res_S20_vsN20 <- results(dds, contrast = c("group2", "S20", "N20"))
sig_S20_vsN20 <- res_S20_vsN20[ which(res_S20_vsN20$padj < 0.05 ), ]

# explore the result table:
# baseMean is the average of the normalized count values, dividing by size factors, 
# taken over all samples. 
# The column log2FoldChange is the effect size estimate. It tells us how much the gene's expression
# seems to have changed due to the temperature change or latitude. This value is reported on
# a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene's expression
# is increased by a multiplicative factor of 2^1.5, which is 2.82. 
# upregulation: log.fc > 0
# downregulation: log.fc < 0

# Make a MA-plot

# In DESeq2, the function plotMA shows the log2 fold changes attributable to a given 
# variable over the mean of normalized counts for all the samples in the DESeqDataSet. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.

plot.new()
plotMA(res_N24_vsN20, ylim=c(-2,2)) # do the same for the other contrasts

# Make a volcano plot
BiocManager::install('EnhancedVolcano')

plot.new()
EnhancedVolcano(res_N24_vsN20,
                lab = rownames(res_N24_vsN20),
                x = 'log2FoldChange',
                y = 'padj',
                labSize = 1.5)

# do the same for the other contrasts


# Make a heatmap with all DEGs

install.packages("pheatmap")

# Combine all unique differentially expressed genes (DEGs) for all contrasts 
# (with unique = each DEG transcript ID occuring only once)
DEGs_N24_vsN20=rownames(sig_N24_vsN20)
DEGs_S24_vsS20=rownames(sig_S24_vsS20)
DEGs_S24_vsN24=rownames(sig_S24_vsN24)
DEGs_S20_vsN20=rownames(sig_S20_vsN20)
allDEGs<-c(DEGs_N24_vsN20,DEGs_S24_vsS20,DEGs_S24_vsN24,DEGs_S20_vsN20)
allDEGsUnique=unique(allDEGs)
select=allDEGsUnique

nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("treat","lat","pop")])

pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Venn diagrams
# make a Venn diagram where you combine the DEGs for both temperature contrasts
# and one where you combine the DEGs for both latitude contrasts

# extract the names of the DEGs like this (do also for other contrasts)
write.table(DEGs_N24_vsN20,file="DEGs_N24_vsN20.txt",row.names=FALSE,col.names=FALSE, sep="\t", quote = FALSE)

https://bioinfogp.cnb.csic.es/tools/venny/



# It can also be useful to examine the counts of reads for a single gene across the groups.'
# A simple function for making this plot is plotCounts, which normalizes counts by sequencing depth and adds a pseudocount of 1/2 to allow for log scale plotting
# The counts are grouped by the variables in intgroup, where more than one variable can be specified. 
# Here we specify the gene which had the smallest p value from the results table created above. 
# You can select the gene to plot by rowname or by numeric index.

# example
d <- plotCounts(dds, gene=which.min(res_N24_vsN20$padj), intgroup="group2", 
                returnData=TRUE)
ggplot(d, aes(x=group2, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


# FUNCTIONAL EXPLORATION (OPTIONAL)

# What do the most significant differtial expressed genes within each of the four contrast code for? 
# And those with the largest log fold change? Blast them using the ncbi translated nucleotide database.

# download the de novo assembly (GSE158138_assembly.fasta.gz) here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158138
# and search for the gene names, copy the associated sequence and blast like this:
# 1. https://blast.ncbi.nlm.nih.gov
# 2. Blastx
# 3. Select Non-redundant protein sequences database
# 4. Paste the transcript in 
