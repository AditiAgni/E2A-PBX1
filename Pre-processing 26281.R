#set wd- session-choose directory-"C:/Users/agnih/OneDrive/Desktop/dissertation/wd/GSE26281"

#cel files unpacking
gse26281<-list.celfiles("GSE26281/", pattern="CEL")


#reading raw data
RD_gse26281<-ReadAffy(verbose = TRUE, filenames = gse26281)

#visualization of raw data before normalization
#hist(RD_gse26281)
#RAWEXP<-exprs(RD_gse26281)
#ma.plot(rowMeans(log2(RAWEXP)), log2(RAWEXP[,1])-log2(RAWEXP[,2]), cex=1)
#boxplot(RD_gse26281)


#phenodata import and read
file.exists("GSE26281_pdata.txt")
PD_gse26281<-read.table("GSE26281_pdata.txt", header = TRUE)
#alternate- read.delim("GSE26281_pdata.txt", header = TRUE, fill = TRUE)
#options(max.print = 380)
#alternate through rio- PD_gse26281<-import("GSE26281_pdata.txt")
#head(Phdata)
#for excel file- PD_gse26281<-import("GSE26281_series_matrix.xlsx")
#head(rio_xlsx)

#viewing phenodata = View(Phdata), for dimensions= dim(file name)


#RMA normalization
ND_gse26281 = rma(RD_gse26281)

#viewing normalized data= View(ND_gse26281) 
#boxplot(ND_gse26281@assayData$exprs)


#getting expression value or probe IDs for genes
ED_gse26281<-as.data.frame(exprs(ND_gse26281))


#merging ED and Phdata for allotment of probe IDs to gene IDs
colnames(ED_gse26281) = PD_gse26281$Sample_IDs
#to view- colnames(ED_gse26281)
#dim(ED_gse26281)- 22283 154


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse26281=rownames(ED_gse26281)
PIDs_gse26281= getBM(attributes = c("affy_hg_u133a_2", "hgnc_symbol"),
                     filters = "affy_hg_u133a_2",
                     values = PIDs.gse26281,
                     mart = Martfunction)
# alternate- if(interactive()){
#  mart <- useEnsembl(biomart = "ensembl",
#                    dataset = "hsapiens_gene_ensembl")
# 
#  getBM(attributes = c("affy_hg_u133a", "hgnc_symbol", "chromosome_name", "band"),
#        filters    = "affy_hg_u133a",
#        values     = "PIDs.gse26281",
#        mart       = mart)
#}
#dim(PIDs_gse26281) [1] 23467     2
#to view lists- listMarts(),  listDatasets(mart= ensembl), listAttributes(mart = ensembl)


#removing unassigned probe IDs
PIDs_gse26281.df=PIDs_gse26281[!(PIDs_gse26281$hgnc_symbol==""),]
#dim(PIDs_gse26281.df)- 22287     2

#bringing rownames to colnames
ED_gse26281.df=rownames_to_column (ED_gse26281, "PIDs")
colnames(PIDs_gse26281.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse26281= merge(PIDs_gse26281.df, ED_gse26281.df, by= "PIDs")
#dim(GS_gse26281)- 22279 85


#limma avereps- to avg out duplicate symbols
AR_gse26281=as.data.frame((limma::avereps(GS_gse26281, GS_gse26281$hgnc_symbol)))
#dim 14117  85

#final dataset - having umique gene symbols as rownames
AR_gse26281=column_to_rownames(AR_gse26281, "hgnc_symbol")
#dim(AR_gse26281)- 14117  84

#removing PIDs column
FD_gse26281=AR_gse26281[-c(1)]
#dim(FD_gse26281)- 14117  83

#creating matrix file of final dataset
FD_gse26281.matrix=as.matrix(FD_gse26281)
class(FD_gse26281.matrix)= "numeric"

#to create common gene symbol column for merging data
MD_gse26281 <- rownames_to_column (FD_gse26281, "gene.symbol")
#dim(FD_gse26281)- 14117  84


#DGE
f.source=factor(PD_gse26281$mutation, levels = c("control", "case"))
design_26281 <- model.matrix(~ 0+factor(f.source))
colnames(design_26281) <-c("control", "case")
fit_26281=lmFit(FD_gse26281.matrix, design_26281)

contrast.matrix <- makeContrasts(case-control, levels = design_26281)
fit2 <- contrasts.fit(fit_26281, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_gse26281<-topTable(fit2, coef=1, adjust="BH", confint=0.95)
DGE2_gse26281<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)

#saving the dge file
write.csv(DGE_gse26281, file = "1DGE_26281.csv")
write.csv(DGE2_gse26281, file = "1DGE2_26281.csv")

#saving final final for meta-analysis
saveRDS(MD_gse26281, file="1MD_gse26281")
