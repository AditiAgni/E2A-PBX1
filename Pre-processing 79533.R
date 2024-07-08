setwd(" ")  #or session-choose directory


#cel files listing
gse79533<-list.celfiles("GSE79533/", pattern="CEL")


#reading raw data
RD_gse79533<-ReadAffy(verbose = TRUE, filenames = gse79533) #containing all samples (n=229)
boxplot(RD_gse79533)


#phenodata import and read
file.exists("GSE79533_updata.xlsx")
UPD_gse79533= import("GSE79533_updata.xlsx") #phenodata having only common samples (n=229)


#filtering samples from expression data; retaining only common samples
keepSamples <- UPD_gse79533$Sample_IDs
rownames(RD_gse79533@phenoData@data) <- sub("_.*", "", rownames(RD_gse79533@phenoData@data))
rownames(RD_gse79533@protocolData@data) <- sub("_.*", "", rownames(RD_gse79533@protocolData@data))
FtRD_gse79533 <- RD_gse79533[, colnames(RD_gse79533) %in% keepSamples]


#RMA normalization
ND_gse79533 = rma(FtRD_gse79533)
boxplot(ND_gse79533@assayData$exprs) #viewing normalized data


#getting expression value or probe IDs for genes
ED_gse79533<-as.data.frame(exprs(ND_gse79533))
dim(ED_gse79533) #54675 219


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse79533=rownames(ED_gse79533)
PIDs_gse79533= getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
                     filters = "affy_hg_u133_plus_2",
                     values = PIDs.gse79533,
                     mart = Martfunction)
PIDs_gse79533.df=PIDs_gse79533[!(PIDs_gse79533$hgnc_symbol==""),] #removing unassigned probe IDs
dim(PIDs_gse79533.df)  #44539     2


#bringing row names to column
ED_gse79533.df=rownames_to_column (ED_gse79533, "PIDs")
colnames(PIDs_gse79533.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse79533= merge(PIDs_gse79533.df, ED_gse79533.df, by= "PIDs")
dim(GS_gse79533) #44539   221


#limma avereps- to avg out duplicate symbols
AR_gse79533=as.data.frame((limma::avereps(GS_gse79533, GS_gse79533$hgnc_symbol)))
dim(AR_gse79533) #22442   221
AR_gse79533=column_to_rownames(AR_gse79533, "hgnc_symbol") #final dataset - having unique gene symbols as row names
dim(AR_gse79533) #22442   220


#removing PIDs column
FD_gse79533=AR_gse79533[-c(1)]
dim(FD_gse79533) #22442   219
FD_gse79533.matrix=as.matrix(FD_gse79533) #converting data frame to matrix for numeric values
class(FD_gse79533.matrix)= "numeric"
dim(FD_gse79533.matrix) #22442   219


#Creating a common gene symbol column for merging data
MD_gse79533 <- rownames_to_column (FD_gse79533, "gene.symbol")
dim(MD_gse79533) #22442   220



#Differential gene expression with E2A-PBX1 as case and other translocations as control
f.source=factor(UPD_gse79533$mutation, levels = c("control", "case"))
design_79533 <- model.matrix(~ 0+factor(f.source))
colnames(design_79533) <-c("control", "case")
fit_79533=lmFit(FD_gse79533.matrix, design_79533)

contrast.matrix <- makeContrasts(case-control, levels = design_79533)
fit2 <- contrasts.fit(fit_79533, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_gse79533<-topTable(fit2, coef=1, adjust="BH", confint=0.95)
DGE2_gse79533<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)


#saving final dataset for meta-analysis
saveRDS(MD_gse79533, file="MD_gse79533.rds")
