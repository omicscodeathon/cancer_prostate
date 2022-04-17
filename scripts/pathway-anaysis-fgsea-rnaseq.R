
#importing read counts from omics
BiocManager::install('edgeR')
BiocManager::install(c("gage","GO.db","AnnotationDbi","org.Hs.eg.db"))
BiocManager::install('reactome.db')
install.packages(c("magrittr","dplyr"))
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
library('edgeR')
library('data.table')
library(reactome.db)
library('org.Hs.eg.db')
library(data.table)
library(fgsea)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(biomaRt)

countdata <- read.table('feature_counts',header = TRUE, skip = 1, row.names = 1)
countdata <- countdata[,-c(1:5)]
colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("..", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("outputs.", "", colnames(countdata), fixed = T)
colnames(countdata) <- gsub("_trimmed_sorted", "", colnames(countdata), fixed = T)
colnames(countdata)

keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]

countdata
sample_info <- data.frame(condition = gsub("_[0-9]+", "", names(countdata)), row.names = names(countdata) )
sample_info$'ethinicity'= c('African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans','African_Americans','African_Americans','African_Americans','African_Americans',
                            'African_Americans','African','African','African','African','African','African','African','African',
                            'African','African','African','African','African','African','African','African','African','African',
                            'African','African','African','African','African','African')
sample_info$'condition'= c('normal_prostate','prostate_tumor','prostate_tumor','normal_prostate','prostate_tumor','normal_prostate',
                           'normal_prostate','prostate_tumor','prostate_tumor','normal_prostate','prostate_tumor','normal_prostate',
                           'prostate_tumor','normal_prostate','prostate_tumor','normal_prostate','prostate_tumor','prostate_tumor',
                           'prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor',
                           'prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor',
                           'prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor',
                           'prostate_tumor','prostate_tumor','prostate_tumor','prostate_tumor')

sample_info$condition <- as.factor(sample_info$condition)
sample_info$ethinicity <- as.factor(sample_info$ethinicity)
design <- as.formula(~condition+ethinicity)
#modelMatrix <- model.matrix(design, data = sample_info)
design

ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata,colData = sample_info,design = design)
vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup=c("condition", "ethinicity"))
ddsObj <- estimateSizeFactors(ddsObj.raw)
ddsObj <- estimateDispersions(ddsObj)
ddsObj <- nbinomWaldTest(ddsObj)
# Run DESeq
ddsObj <- DESeq(ddsObj.raw)
resultsNames(ddsObj)

res <- results(ddsObj, alpha=0.05)

resNvP <- res

listids <- as.list(rownames(resNvP)) 
#Coverting HGNC gene names to ensembl IDS

mart <- useMart('ENSEMBL_MART_ENSEMBL',host="https://www.ensembl.org")
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id'),
  uniqueRows = TRUE)
table(listids %in% annotLookup$hgnc_symbol)
for (i in rownames(resNvP)){
  nm <-  annotLookup$ensembl_gene_id[annotLookup$hgnc_symbol == i]
  if (!is_empty(nm)){
    rownames(resNvP)[rownames(resNvP) == i] <- nm
  }
}

#baCK TO THE PIPELINE

# view the available databases
listMarts(host="https://www.ensembl.org")
## set up connection to ensembl database
ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="https://www.ensembl.org")
listDatasets(ensembl) %>% 
  filter(str_detect(description, "Human"))

ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
listFilters(ensembl) %>% 
  filter(str_detect(name, "ensembl"))


#Filtering
ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(resNvP)[1:1000]
listAttributes(mart) %>% 
  head(20)
attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name')

# run the query
annot <- getBM(attributes=attributeNames, 
               filters = ourFilterType, 
               values = filterValues, 
               mart = ensembl)
head(annot)
annot %>%  
  add_count(ensembl_gene_id) %>%  
  filter(n>1)

#Querying for the entire table
filterValues <- rownames(resNvP)
attributeNames <- c('ensembl_gene_id', 'c', 'external_gene_name')
annot <- getBM(attributes=attributeNames, 
               filters = ourFilterType, 
               values = filterValues, 
               mart = mart)
colnames(annot) 
colnames(annot)[colnames(annot)=="ensembl_gene_id"] <- "GeneID"

annotNvP <- as.data.frame(resNvP) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)

#skipping the Shrink
ddsShrink <- lfcShrink(ddsObj,coef="condition_prostate_tumor_vs_normal_prostate")
shrinkNvP <- as.data.frame(ddsShrink) %>%
  rownames_to_column("GeneID") %>% 
  left_join(annot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)
#Shrink yielded NA values because the rownames have not been changed to EntrezId
gseaDat <- filter(annotNvP, !is.na(entrezgene_id))
rankData <- gseaDat$logFC
names(rankData) <- gseaDat$entrezgene_id
head(rankData)





library(reactome.db)
BiocManager::install('reactome.db')
library('org.Hs.eg.db')
#Loaing the pathway db
human.universe <- keys(org.Hs.eg.db, "ENTREZID")


# Selecting reactome gene sets for 
pathways <- na.omit(select(reactome.db, keys=human.universe, c("PATHID"),
                           keytype = 'ENTREZID'))
pathways <- split(pathways$ENTREZID, pathways$PATHID)

pathway2name <- as.data.table(na.omit(select(reactome.db, names(pathways),
                                             c("PATHNAME"), 'PATHID')))
# Remove organism prefix
pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
pathway2name <- setNames(pathway2name$PATHNAME, pathway2name$PATHID)
pathway2name <- iconv(pathway2name, "latin1", "ASCII", sub="")

pathway.lines <- sapply(names(pathways), function(p) {
  link <- p
  name <- paste0(p, "_", pathway2name[p])
  name <- gsub("[ ()/]+", "_", name)
  sprintf("%s\t%s\t%s", name, link, paste0(pathways[[p]], collapse="\t"))
})

#Writing the pathway files.
write(pathway.lines, file="pathway.lines")

pathwayLines <- strsplit(readLines("pathway.lines"), "\t")
inbuiltPathways <- lapply(pathwayLines, tail, -2)
names(inbuiltPathways) <- sapply(examplePathways, head, 1)


#Loading the pathway files. 
#using the downloaded gmt files from humans
pathwayLines <- strsplit(readLines("Human_AllPathways_April_01_2022_entrezgene.gmt"), "\t")
newPathways <- lapply(pathwayLines, tail, -2)
names(newPathways) <- sapply(pathwayLines, head, 1)

fgseaRes <- fgsea(pathways = newPathways, 
                  stats    = rankData,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes)
table(is.na(fgseaRes$pathway))
head(fgseaRes[order(pval), ])
fgseaRes <- fgsea(pathways = newPathways, 
                  stats    = rankData,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ],20)
