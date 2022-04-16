#This script performs differential gene expression analysis and pathaway analysis
#The following packages should be installed to perform differential gene expression analysis with DESeq2 as well as pathway enrichment with kegg.
#Uncomment to install.

#install.package("magrittr")
#install.packages("NMF")
#install.packages("ggdendro")
#install.packages("gridExtra")
#install.packages("gtable")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("apeglm")
#BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Hs.eg.db"))
#BiocManager::install('EnhancedVolcano', force=T)


#Loading the library magrittr
library(magrittr)

#get the table of read counts by indicating the path to the file
read.counts <- read.table("feature_counts.txt", header = TRUE, row.names = 1)
head(read.counts)

#one of the requirement of of the assay() slot is that the row.names corresponds to the gene ID and the col.name to the sample names
rownames(read.counts)
colnames(read.counts)

#exclude all column that do not contain read counts
readcounts <- read.counts[ , c(6:21)]
head(readcounts)
colnames(readcounts)
#Renaming the columns in readcounts to sample ID only
colnames(readcounts)= c('SRR6059552','SRR6059554','SRR6059560','SRR6059562','SRR6059564','SRR6059566','SRR6059568',
                        'SRR6059570','SRR6059572','SRR6059574','SRR6059576','SRR6059578','SRR6059580','SRR6059582',
                        'SRR6059600','SRR6059602')
head(readcounts)        
#checking the data
str(readcounts)
head(readcounts)

#Making the sample_info dataframe with metadata where row.names should match the indivdual samplenames
sample_info <- data.frame(condition = gsub("_[0-9]+", "", names(readcounts)), row.names = names(readcounts))
head(sample_info)
#Adding more columns names for our sample_info dataframe
sample_info$'ethinicity'= c('African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans')
sample_info$'condition'= c('normal_prostate','prostate_tumor','prostate_tumor','normal_prostate','prostate_tumor','normal_prostate',
                           'normal_prostate','prostate_tumor','prostate_tumor','normal_prostate','prostate_tumor','normal_prostate',
                           'prostate_tumor','normal_prostate','prostate_tumor','normal_prostate')

#converting the values in the added columns to characters
sample_info$condition <- as.factor(sample_info$condition)
sample_info$ethinicity <- as.factor(sample_info$ethinicity)
head(sample_info)

#Checking if the name of the columns in the count matrix (readcounts) in the same order as that in name of the sample data (sample_info)
all(colnames(readcounts) == rownames(sample_info))
all(colnames(readcounts) %in% rownames(sample_info))

#Loading the DESeq2 Library
library(DESeq2)

#Generating the DESeqDataSet for DGE analysis
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = sample_info,
                                   design = ~ condition)

#checking results using accessors described above
head(colData(DESeq.ds))
head(assay(DESeq.ds, "counts")) 
head(rowData(DESeq.ds))

#test what counts() returns
str(counts(DESeq.ds))

#removing genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))> 0, ]

#investigate different library sizes
colSums(counts(DESeq.ds))

#calculating the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

#check to see that it now contains the sizeFactors
colData(DESeq.ds)

#use counts() to immediately retrieve the normalized read count
counts.sf.normalized <- counts(DESeq.ds, normalized = TRUE)

#transform size factor normalized readcounts to log2 scale using pseudo count of 1
log.norm.counts <- log2(counts.sf.normalized + 1)

#plotting to visualize non-transformed and log2 transformed read counts
par(mfrow=c(2,1))

#plotting a box plot for non-transformed read counts
boxplot(counts.sf.normalized, notch = TRUE,
        main = "untransformed read counts", ylab = "read counts",
        col = "blue", points='orange')

#plotting a box plot for log2-transformed read counts
boxplot(log.norm.counts, notch = TRUE,
        main = "log2-transformed read counts",
        ylab = "log2 read counts")

#Heatmap of the count matrix to explore a count matrix
library(pheatmap)
select <- order(rowMeans(counts(DESeq.ds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(DESeq.ds)[,c("condition","ethinicity")])
ntd <- normTransform(DESeq.ds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#visualizing the dataset in a PCA
library(ggplot2)

#performing Variance stabilizing transformation
#performing Regularized log transformation
DESeq.vst = vst(DESeq.ds, blind = F, fitType='mean')
DESeq.rlog = rlog(DESeq.ds, blind = F)

par(mfrow=c(1,2))
p <- plotPCA(DESeq.vst)
p <- p + theme_bw() + ggtitle("vst transformed counts")
print(p)

q <- plotPCA(DESeq.rlog)
q <- p + theme_bw() + ggtitle("rlog transformed counts")
print(q)

#Running differential expression analysis using DESeq2
#setting control condition as reference
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, ref="normal_prostate")

#Running Differential expression analysis using the DESeq function
DGE <- DESeq(DESeq.ds)
#Viewing all comparisons 
resultsNames(DGE)

#Using the result() function to extract the base means across samples 
DGEresults <- results(DGE,contrast = c("condition","prostate_tumor","normal_prostate"), independentFiltering = TRUE, alpha = 0.1)
summary(DGEresults)
# get gene expression table 
head(DGEresults)

#Ordering gene expression table by adjusted p value (Benjamini-Hochberg FDR method) 
DGEresults[order(DGEresults$padj),]  

#converting the results to a dataframe
head(DGEresults)
table(DGEresults$padj < 0.05)
rownames(subset(DGEresults, padj < 0.05))

#Exporting differential gene expression analysis table to CSV file,
write.csv(as.data.frame(DGEresults[order(DGEresults$padj),] ), file="condition_normal_prostate_vs_prostate_tumor.csv")

#Get summary of differential gene expression with adjusted p value cut-off at 0.05,
summary(results(DGE, alpha=0.05))

#Representing results in plots
#Histogram showing p values
hist(DGEresults$pvalue,
     col = "#0080ff", border = "white", xlab = "P value", ylab = "Frequency",
     main = "Frequencies of p-values")

#Shrinkage estimation of log2 fold changes (LFCs)
#The shrinkage of effect size (LFC) helps to remove the low count genes (by shrinking towards zero).
#The low or highly variable read count genes can give large estimates of LFCs which may not represent true difference in changes in 
#gene expression between two conditions.
library(apeglm)

resLFC <- lfcShrink(DGE, coef="condition_prostate_tumor_vs_normal_prostate", type="apeglm")
head(resLFC)

#Visualize the shrinkage estimation of LFCs with MA plot and compare it without shrinkage of LFCs
par(mfrow = c(1, 2))
plotMA(resLFC, main="Shrinkage of LFCs")
plotMA(DGEresults, main="No shrinkage of LFCs")

# Load libraries to view MA plot using gglot2
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(DGEresults)

# Examining this data frame
head(deseq2ResDF)

# Setting a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# Plotting the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + 
  scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + theme_bw()

# adding some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + geom_point(size=1) + 
  scale_y_continuous(limits=c(-3, 3), oob=squish) + scale_x_log10() + 
  geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") +
  labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') +
  theme_bw() + geom_density_2d(colour="black", size=2)

#plotting a Heatmap 
library(NMF)
DGEresults.sorted <- DGEresults[order(DGEresults$padj), ]
DGEgenes <- rownames(subset(DGEresults.sorted, padj < 0.05))

library("RColorBrewer")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]
aheatmap(hm.mat_DGEgenes, Rowv = NA, Colv = NA, color="Blues")
aheatmap(hm.mat_DGEgenes,Rowv = "correlation", Colv = "correlation",
         distfun = "homosapien", hclustfun = "average", scale = "row",color="YlOrRd")

#Transforming the count data using the variance stablilizing transform
deseq2VST <- vst(DESeq.ds, blind = F)

# Converting the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keeping only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
sigGenes

# Converting the VST counts to long format for ggplot2
library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwriting our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Making a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + 
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


#Converting the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Computing a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Clustering the based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construcing a dendogram for samples
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factoring samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Constructing the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") +
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combining the dendrogram and the heatmap
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))

# Loading in libraries necessary for modifying plots
library(gtable)
library(grid)

# Modifying the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Converting both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Checking the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Adding in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Making sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arranging the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Drawing the plot
grid.draw(finalGrob)
sampleData = sample_info

# Re-ordering the sample data to match the clustering we did
sampleData$Run <- factor(row.names(sampleData), levels=clusterSample$labels[clusterSample$order])

# Constructing a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52")
sampleClinical <- ggplot(sampleData, aes(x=Run, y=1, fill=condition)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + 
scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

# Converting the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Making sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arranging the output in the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

#  making a dendrogram with ggplot
geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + 
scale_x_continuous(expand=c(0, 0)) + theme_dendro()

#Step 2: Re-arranging the heatmap cells

# re-factoring genes for ggplot2
deseq2VST$Gene <- factor(deseq2VST$Gene, levels=clusterGene$labels[clusterGene$order])

# recreating the heatmap with this new factoring
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + 
theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())

#Step 3: converting to everything to grobs 
# convert the heatmap to a grob
heatmapGrob <- ggplotGrob(heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)))

# converting the dendrogram to a grob
# flipping the axis above so the x-axis is now displayed as what we would think of as the y-axis
geneDendrogramGrob <- ggplotGrob(geneDendrogram + scale_x_discrete(expand=c(0, 0)))

#Step 4: aligning the gene dendrograms to match the heatmap 
# checking that the both the heatmap and gene dendrogram have the same number of vertical elements
length(heatmapGrob$heights) == length(geneDendrogramGrob$heights)

# make sure every height between the two grobs is the same
maxHeight <- unit.pmax(geneDendrogramGrob$heights, heatmapGrob$heights)
geneDendrogramGrob$heights <- as.list(maxHeight)
heatmapGrob$heights <- as.list(maxHeight)

# Step 4b: we have a new heatmap so we need to re-align the horizontal elements 
# repeating the steps in the tutorial
# checking the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths
sampleClinicalGrob$widths

# adding in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# making sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Step 5: creating a blank panel 
# we can use grid graphics for this
blankPanel <- grid.rect(gp=gpar(col="white"))

# Step 6: Arranging the final result 
# arranging all the plots together
finalGrob_v2 <- arrangeGrob(blankPanel, sampleDendrogramGrob, blankPanel, sampleClinicalGrob, geneDendrogramGrob, heatmapGrob, ncol=2, nrow=3, widths=c(1,5), heights=c(2,.8,6))

# drawing the final result
grid.draw(finalGrob_v2)

## plotting a volcano plot
library(EnhancedVolcano)
EnhancedVolcano(DGEresults,
                lab = rownames(DGEresults), x = 'log2FoldChange', y = 'pvalue')

##Performing pathway analysis
library(gage)
tumor_v_normal_DE <- results(DGE,
                             contrast=c("condition", "prostate_tumor", "normal_prostate"))
# setting up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# setting up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# loading in libraries to annotate data
library(AnnotationDbi)
library(org.Hs.eg.db)

#annotating the deseq2 results with additional gene identifiers
tumor_v_normal_DE$symbol <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="SYMBOL", keytype="SYMBOL", multiVals="first")
tumor_v_normal_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="ENTREZID", keytype="SYMBOL", multiVals="first")
tumor_v_normal_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(tumor_v_normal_DE), column="GENENAME", keytype="SYMBOL", multiVals="first")

# grabbing the log fold changes for everything
tumor_v_normal_DE.fc <- tumor_v_normal_DE$log2FoldChange
names(tumor_v_normal_DE.fc) <- tumor_v_normal_DE$entrez

# Running enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(tumor_v_normal_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(tumor_v_normal_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(tumor_v_normal_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(tumor_v_normal_DE.fc, gsets = go.cc.gs)

# converting the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# converting the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# getting the list of genes associated with pathways
a <- function(sigPathways, geneSetDB){
  
  # grabbing the names of all significant pathways
  sigPathwaysID <- rownames(sigPathways)
  
  # subsetting the geneSet list to only those in the significant pathways
  sigPathwaysGenes <- geneSetDB[which(names(geneSetDB) %in% sigPathwaysID)]
  numberSigPathways <- length(sigPathwaysGenes)
  sigPathwaysGenes <- unlist(sigPathwaysGenes)
  
  # counting the number of times a gene occurs in these significant pathways
  sigPathwaysTable <- plyr::count(sigPathwaysGenes)
  
  # annotating these final results with the gene symbol and some extra information
  sigPathwaysTable$symbol <- mapIds(org.Hs.eg.db, keys=as.character(sigPathwaysTable$x), column="SYMBOL", keytype="ENTREZID", multiVals="first")
  sigPathwaysTable$proportion <- sigPathwaysTable$freq/numberSigPathways
  
  # sorting and returning the results
  sigPathwaysTable <- sigPathwaysTable[order(-sigPathwaysTable$proportion),]
  
  return(sigPathwaysTable)
}

# running the function defined above
geneTable <- a(na.omit(fc.go.bp.p.up[fc.go.bp.p.up$q.val <= .05,]), go.bp.gs)
geneTable
write.csv(as.data.frame(geneTable ), file="genes_associated_pathway.csv")

# Installing the pathview from bioconductor
library(pathview)
####################################################################################################################
# Viewing the hsa05215 (prostate cancer) and hsa59351 pathway from the pathway analysis 
fc.kegg.sigmet.p.up[grepl("hsa05215", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]
fc.kegg.sigmet.p.down[grepl("hsa59351", rownames(fc.kegg.sigmet.p.down), fixed=TRUE),]

# Overlaying the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa05215",kegg.native=T)      
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa05215")

#####################################################################################################################
# Viewing the hsa00140 pathway (Steroid hormone biosynthesis) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa00140", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlaying the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa00140", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa00140")

#####################################################################################################################
# Viewing the hsa04010 pathway (	MAPK signaling pathway) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04010", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlayying the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04010", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04010")

#####################################################################################################################
# View the hsa04060 pathway (	Cytokine-cytokine receptor interaction) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04060", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04060", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04060")

####################################################################################################################
# View the hsa04110 pathway (Cell cycle) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04110", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04110", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04110")

####################################################################################################################
# View the hsa04110 pathway (Drug metabolism cytochrome) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa00982", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa00982", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa00982")

####################################################################################################################
# View the hsa04110 pathway (Drug metabolism cytochrome) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04610", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04610", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04610")

####################################################################################################################
# View the hsa04110 pathway (TGF-beta signaling pathway) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04350", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04350", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04350")

####################################################################################################################
# View the hsa04110 pathway (Drug metabolism cytochrome) from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa04390", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04390", kegg.native=FALSE)
pathview(gene.data=tumor_v_normal_DE.fc, species="hsa", pathway.id="hsa04390")
