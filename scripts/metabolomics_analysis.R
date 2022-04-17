# A Deep-learning enabled, personalized pathway-based R package for diagnosis and prognosis predictions using metabolomics data
# Uncomment to install
# install.packages("lilikoi")
# install.packages("pathview")

# Loading the installed package
library(lilikoi)

## Load the data with the Loaddata function:
dt <- lilikoi.Loaddata(file="prostateEdited.csv")
Metadata <- dt$Metadata
Metadata$Label <- as.factor(Metadata$Label)
class(Metadata$Label)
dataSet <- dt$dataSet

## Check the preprocess functions:
lilikoi.preproc_norm(inputdata=Metadata, method="standard")
lilikoi.preproc_knn(inputdata=Metadata, method="knn")


## Transform the metabolite names to the HMDB ids using Lilikoi MetaTOpathway function
convertResults=lilikoi.MetaTOpathway('name')
Metabolite_pathway_table = convertResults$table
head(Metabolite_pathway_table)

## Transform metabolites into pathway using Pathtracer algorithm
PDSmatrix= lilikoi.PDSfun(Metabolite_pathway_table)

## Using PDSfun, we generate a new PDSmatrix( a normalized score in the range [0,1] that measures the degree of dysregulation of 
## a pathway relative to the norm (controls))  which has pathways as columns instead of metabolites.
head(t(PDSmatrix))
dim(t(PDSmatrix))

## Select the most significant pathway related to phenotype.
selected_Pathways_Weka= lilikoi.featuresSelection(PDSmatrix,threshold= -0.5,method="info")

## Metabolites-pathway regression
## Before running this function, users need to download Cytoscape from http://www.cytoscape.org/download.php.
## Keep Cytoscape running, then you can connect to R package RCy3.
## Load RCy3 library 
library(RCy3)
lilikoi.meta_path(PDSmatrix = PDSmatrix, selected_Pathways_Weka = selected_Pathways_Weka, Metabolite_pathway_table = Metabolite_pathway_table, pathway = "Metabolic Pathways")

## Machine learning using the metabolite features as the inputs for deep learning–based classification, along with other 
## popular methods: LDA, SVM, RF, RPART, LOG, and GBM (Methods). The objective is to distinguish the controls from samples

dt <- lilikoi.Loaddata(file="prostateEdited.csv")
Metadata <- dt$Metadata
Metadata$Label <- as.factor(Metadata$Label)
class(Metadata$Label)
dataSet <- dt$dataSet

Metadata$Label <- make.names(c(Metadata$Label))
library(h2o)
lilikoi.machine_learning(MLmatrix = Metadata, measurementLabels = Metadata$Label,
significantPathways = 0,
trainportion = 0.8, cvnum = 10, dlround=50,nrun=10, Rpart=FALSE,
LDA=FALSE,SVM=FALSE,RF=FALSE,GBM=TRUE,PAM=TRUE,LOG=TRUE,DL=TRUE)

## KEGG (pathview) plot
## Reload the dat
prostate_data <- read.csv('prostateEdited.csv', check.names=F, row.names=1)
sampleInfo <- prostate_data$Label
names(sampleInfo) <- row.names(prostate_data)

metamat <- t(t(prostate_data[-1]))
metamat <- log2(metamat)
grouporder <- c('Normal','Cancer')

# If encounter error: .onLoad failed in loadNamespace() for 'org.Hs.eg.db', run this line first:
# options(connectionObserver = NULL)

library(pathview)
lilikoi.KEGGplot(metamat = metamat, sampleinfo = sampleInfo, grouporder = grouporder,
                 pathid = '05215', specie = 'hsa',
                 filesuffix = 'PR000570', 
                 Metabolite_pathway_table = Metabolite_pathway_table)

#Thus, from the KEGGplot function, a graph called "hsa05215" .GSE16873.png" has been saved at your working direction

