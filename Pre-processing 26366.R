setwd(" ")  #or session-choose directory

#cel files listing
gse26366<-list.celfiles(" ", pattern="CEL")


#reading raw data
RD_gse26366<-ReadAffy(verbose = TRUE, filenames = gse26366) #contain all samples (n=229)
boxplot(RD_gse26366)


#phenodata import and read
file.exists("GSE26366_updata.xlsx")
UPD_gse26366= import("GSE26366_updata.xlsx") #phenodata having only common samples (n=109)


#filtering samples from expression data; retaining only common samples
keepSamples <- UPD_gse26366$Sample_ID
rownames(RD_gse26366@phenoData@data) <- sub("_.*", "", rownames(RD_gse26366@phenoData@data))
rownames(RD_gse26366@protocolData@data) <- sub("_.*", "", rownames(RD_gse26366@protocolData@data))
FtRD_gse26366 <- RD_gse26366[, colnames(RD_gse26366) %in% keepSamples]


#RMA normalization
ND_gse26366 = rma(FtRD_gse26366)
boxplot(ND_gse26366@assayData$exprs) 


#getting expression value or probe IDs for genes
ED_gse26366<-as.data.frame(exprs(ND_gse26366))
dim(ED_gse26366) #22283 109


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse26366=rownames(ED_gse26366)
PIDs_gse26366= getBM(attributes = c("affy_hg_u133a_2", "hgnc_symbol"),
                     filters = "affy_hg_u133a_2",
                     values = PIDs.gse26366,
                     mart = Martfunction)
dim(PIDs_gse26366) #23652     2
PIDs_gse26366.df=PIDs_gse26366[!(PIDs_gse26366$hgnc_symbol==""),] #removing unassigned probe IDs
dim(PIDs_gse26366.df) #22454     2


#bringing row names to column
ED_gse26366.df=rownames_to_column (ED_gse26366, "PIDs")
colnames(PIDs_gse26366.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse26366= merge(PIDs_gse26366.df, ED_gse26366.df, by= "PIDs")
dim(GS_gse26366) #22454 111


#limma avereps- to avg out duplicate symbols
AR_gse26366=as.data.frame((limma::avereps(GS_gse26366, GS_gse26366$hgnc_symbol)))
dim(AR_gse26366) #14215 111


#final dataset - having unique gene symbols as row names
AR_gse26366=column_to_rownames(AR_gse26366, "hgnc_symbol")
dim(AR_gse26366) #14215 110
FD_gse26366=AR_gse26366[-c(1)] #removing PIDs column
dim(FD_gse26366) #14215 109
FD_gse26366.matrix=as.matrix(FD_gse26366) #creating matrix file of final dataset
class(FD_gse26366.matrix)= "numeric"


##Creating a common gene symbol column for merging data
MD_gse26366 <- rownames_to_column (FD_gse26366, "gene.symbol")
dim(MD_gse26366) #14215 110


#Differential gene expression with E2A-PBX1 as case and other translocations as controls
f.source=factor(PD_gse26366$mutation, levels = c("control", "case"))
design_26366 <- model.matrix(~ 0+factor(f.source))
colnames(design_26366) <-c("control", "case")
fit_26366=lmFit(FD_gse26366.matrix, design_26366)

contrast.matrix <- makeContrasts(case-control, levels = design_26366)
fit2 <- contrasts.fit(fit_26366, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_gse26366<-topTable(fit2, coef=1, adjust="BH", confint=0.95)
DGE2_gse26366<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)


#saving final final for meta-analysis
saveRDS(MD_gse26366, file="MD_gse26366")
