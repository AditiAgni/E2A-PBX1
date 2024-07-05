setwd(" ")  #or session-choose directory

#cel files download and unpacking
gse26281<-list.celfiles(" ", pattern="CEL")


#reading raw data
RD_gse26281<-ReadAffy(verbose = TRUE, filenames = gse26281) #Containing all samples (n=154)
boxplot(RD_gse26281) #visualization of raw data before normalization


#phenodata import and read
file.exists("GSE26281_updata.xlsx")
UPD_gse26281= import("GSE26281_updata.xlsx") #phenodata having only common samples (n=62)
View(UPD_gse26281)


#filtering samples from expression data; retaining only common samples
keepSamples <- UPD_gse26281$Sample_IDs
rownames(RD_gse26281@phenoData@data) <- sub("_.*", "", rownames(RD_gse26281@phenoData@data))
rownames(RD_gse26281@protocolData@data) <- sub("_.*", "", rownames(RD_gse26281@protocolData@data))
FtRD_gse26281 <- RD_gse26281[, colnames(RD_gse26281) %in% keepSamples]


#RMA normalization
ND_gse26281 = rma(FtRD_gse26281)
boxplot(ND_gse26281@assayData$exprs) 


#getting expression value or probe IDs for genes
ED_gse26281<-as.data.frame(exprs(ND_gse26281))
dim(ED_gse26281) #22283 62


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse26281=rownames(ED_gse26281)
PIDs_gse26281= getBM(attributes = c("affy_hg_u133a_2", "hgnc_symbol"),
                     filters = "affy_hg_u133a_2",
                     values = PIDs.gse26281,
                     mart = Martfunction)
dim(PIDs_gse26281) #23652     2


#removing unassigned probe IDs
PIDs_gse26281.df=PIDs_gse26281[!(PIDs_gse26281$hgnc_symbol==""),]
dim(PIDs_gse26281.df) #22287     2

#bringing row names to column names
ED_gse26281.df=rownames_to_column (ED_gse26281, "PIDs")
colnames(PIDs_gse26281.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse26281= merge(PIDs_gse26281.df, ED_gse26281.df, by= "PIDs")
dim(GS_gse26281) #22454 64


#limma avereps- to avg out duplicate symbols
AR_gse26281=as.data.frame((limma::avereps(GS_gse26281, GS_gse26281$hgnc_symbol)))
dim(AR_gse26281) #14215  64


#final dataset - having unique gene symbols as row names
AR_gse26281=column_to_rownames(AR_gse26281, "hgnc_symbol")
dim(AR_gse26281) #14215  63
FD_gse26281=AR_gse26281[-c(1)] #removing PIDs column
dim(FD_gse26281) #14215  62
FD_gse26281.matrix=as.matrix(FD_gse26281) #creating matrix file of final dataset
class(FD_gse26281.matrix)= "numeric"


#Creating common gene symbol column for merging data later
MD_gse26281 <- rownames_to_column (FD_gse26281, "gene.symbol")
dim(FD_gse26281) #14215  62


#Differential gene expression with E2A-PBX1 as case and other translocations as controls
f.source=factor(UPD_gse26281$mutation, levels = c("control", "case"))
design_26281 <- model.matrix(~ 0+factor(f.source))
colnames(design_26281) <-c("control", "case")
fit_26281=lmFit(FD_gse26281.matrix, design_26281)

contrast.matrix <- makeContrasts(case-control, levels = design_26281)
fit2 <- contrasts.fit(fit_26281, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_gse26281<-topTable(fit2, coef=1, adjust="BH", confint=0.95)
DGE2_gse26281<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)


#saving final final for meta-analysis
saveRDS(MD_gse26281, file="MD_gse26281")
