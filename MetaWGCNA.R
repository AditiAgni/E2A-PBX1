setwd(" ") #or set session directory

library(tidyverse)
library(dbplyr)
library(WGCNA)
library(rio)
library(sva)
library(ggplot2)
library(ggfortify)
library(limma)


#Importing datasets saved after pre-processing
GSE26281= readRDS("MD_gse26281.rds")
GSE26366= readRDS("MD_gse26366.rds")
GSE79533= readRDS("MD_gse79533.rds")



############# CREATING METAFILE #################



#merging multiple data frames
datasetlist = list(GSE26281, GSE26366, GSE79533)
Metafile <- datasetlist %>% reduce(inner_join, by='gene.symbol')
dim(Metafile) #14215   391

#making model matrix with batches of pheno-data for batch correction
Pmeta= import("Pmeta.xlsx")
modphenoBC = model.matrix(~1, data=Pmeta)
PhenoBatch = Pmeta$Dataset_ID

#converting gene symbols to row names and forming uniform matrix of metafile
Metafile= column_to_rownames(Metafile, var="gene.symbol")
Metafile.numeric = mutate_all(Metafile, as.numeric)

#batch correction 
MetaBC = ComBat(dat= Metafile.numeric, batch=PhenoBatch, mod = modphenoBC)
dim(MetaBC) #14215   390


# PCA plot before batch correction
pca_before<-prcomp(t(Metafile.numeric))
autoplot(pca_before, data=Pmeta, colour= "Dataset_ID" ,main="Before batch correction", scale. = TRUE)+ theme_bw()

# PCA plot after batch effects removal
pca_after<-prcomp(t(MetaBC))
autoplot(pca_after, data=Pmeta, colour= "Dataset_ID", main="After batch correction", scale. = TRUE)+ theme_bw()


#Differential gene expression with E2A-PBX1 as case and other translocations as control
f.source=factor(Pmeta$mutation, levels = c("control", "case"))
design_metafile <- model.matrix(~ 0+factor(f.source))
colnames(design_metafile) <-c("control", "case")
fit_metafile=lmFit(MetaBC, design_metafile)

contrast.matrix <- makeContrasts(case-control, levels = design_metafile)
fit2 <- contrasts.fit(fit_metafile, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_metafile<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)

#adding differential expression column in table
DGE_metafile$DE = "NO"
DGE_metafile$DE[DGE_metafile$logFC>1 & DGE_metafile$adj.P.Val<10^-25] = "UP"
DGE_metafile$DE[DGE_metafile$logFC< -1 & DGE_metafile$adj.P.Val<10^-25] = "DOWN"

#adding top 10 DEG names in separate column in table
DGE_metafile=rownames_to_column(DGE_metafile, var = "symbol")
DGE_metafile$GeneSymbol <-ifelse(DGE_metafile$symbol %in% head(DGE_metafile[order(DGE_metafile$adj.P.Val), 
                                                                            "symbol"], 50), DGE_metafile$symbol, NA)

#ggplot for 30 top expressed genes
ggplot(data = DGE_metafile, aes(x = logFC, y = -log10(adj.P.Val), col = DE, label=GeneSymbol)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(10^-25), col = "black", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  scale_color_manual(values = c("blue", "grey", "red"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_text(check_overlap = TRUE, hjust=0.5,vjust=-1)

saveRDS(MetaBC, file = "Metadata.rds")



########### DATA PRE-PROCESSING for WGCNA ################



# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
datExpr= as.data.frame(t(MetaBC))

# removing genes with missing values
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #TRUE

# Updating fusion gene names
Pmeta$Types=gsub("E2A-PBX1", "TCF3::PBX1", Pmeta$Types)
Pmeta$Types=gsub("BCR-ABL1", "BCR::ABL1", Pmeta$Types)
Pmeta$Types=gsub("TEL-AML1", "TEL::AML1", Pmeta$Types)

# binarize categorical variables for creating levels in pheno-data
f.source=factor(Pmeta$Types, levels = c("TEL::AML1", "TCF3::PBX1", "MLL", "Hyperdiploid", "BCR::ABL1"))
pdatabi <- binarizeCategoricalColumns(f.source, includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      dropFirstLevelVsAll = FALSE,
                                      includePrefix = FALSE,
                                      prefixSep = ".", nameForAll = "")
pdata <- cbind(Pmeta, pdatabi) #to combine GSM IDs



# outlier detection; sample clustering
sizeGrWindow(22,9) #width = 22, height = 9
par(cex = 1)
par(mar = c(0,4,2,0))
sampleTree1 = hclust(dist(datExpr), method = "average")
plot(sampleTree1, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 120, col = "red")

# Determine cluster under the line to remove outliers
clust = cutreeStatic(sampleTree1, cutHeight = 120, minSize = 10)
table(clust) #samples we want to keep- 709 #removed 4 outliers
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2, main = "Sample clustering without outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# removing outliers from pheno-data
outliers= c("GSM2097329", "GSM645439", "GSM647339", "GSM2097413", "GSM647305", "GSM645388")
datTraits=pdata %>%
  column_to_rownames(var = "Sample_IDs")%>%
  filter(!row.names(.)%in% outliers)
datTraits= datTraits[,-c(1:4)]

# checking if samples and their order match in expression data and pheno-data
all(rownames(datTraits) %in% rownames(datExpr))
all(rownames(datTraits) == rownames(datExpr))
#should return true for both


####### Re-cluster samples with clinical traits
sampleTree3 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree3, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")



############### SIGNED NETWORK CONSTRUCTION #################



# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(12,9)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off 
abline(h = 0.8, col="red") 


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



############# MODULE CONSTRUCTION & TRAIT RELATIONSHIPS ###############


# Constructing signed network in single block using sft power 12 
cor <- WGCNA::cor
bwnet = blockwiseModules(datExpr, maxBlockSize = 15000,
                         power = 12, TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3, networkType = "signed")


# get number of genes for each module and plot dendrogram
table(bwnet$colors)
sizeGrWindow(12, 9)
plotDendroAndColors(bwnet$dendrograms[[1]], bwnet$colors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


## Relating modules with external traits

moduleColors = bwnet$colors
MEs = bwnet$MEs
# Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


## Display the correlation values within a heatmap plot
sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.45,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



################## Intra-modular analysis ###################



# Define variable weight containing the TCF3::PBX1 column of datTrait
TCF3_PBX1 = as.data.frame(datTraits$`TCF3::PBX1`)
names(TCF3_PBX1) = "TCF3::PBX1"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

# calculating module membership and p-value by using Pearson correlation between expression data and module eigengens
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# calculating gene significance and p-value
geneTraitSignificance = as.data.frame(cor(datExpr, TCF3_PBX1, use = "p")) #pearson method
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(TCF3_PBX1), sep="")
names(GSPvalue) = paste("p.GS.", names(TCF3_PBX1), sep="")


# analysing module associated with TCF3::PBX1
# positively correlated
module = "midnightblue"
column = match(module, modNames)
moduleGenes = moduleColors==module
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for TCF3::PBX1",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                   col = module, abline = TRUE)
midnightblue=names(datExpr)[moduleColors=="midnightblue"]
                

## Ordering genes according to significance
# modulewise 
Genesymbols= midnightblue
module = "midnightblue"
geneInfo = data.frame(geneSymbol = Genesymbols,
                       moduleColor = module,
                       geneTraitSignificance=abs(geneTraitSignificance[moduleGenes, 1]),
                       GSPvalue=abs(GSPvalue[moduleGenes, 1]))
write.csv(geneInfo, file = "MidnightBlue genes.csv")


## Module Eigengenes for QC

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, TCF3_PBX1))

# Plot the dendrogram
sizeGrWindow(15,10)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 0.7)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(6,6,2,1),
                      plotDendrograms = FALSE)
