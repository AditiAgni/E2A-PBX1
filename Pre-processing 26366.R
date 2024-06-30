#set wd- session-choose directory-"C:/Users/agnih/OneDrive/Desktop/dissertation/wd/GSE26366"

#cel files download and unpacking
gse26366<-list.celfiles("GSE26366/", pattern="CEL")


#reading raw data
RD_gse26366<-ReadAffy(verbose = TRUE, filenames = gse26366)

#visualization of raw data before normalization
#hist(RD_gse26366)
#RAWEXP<-exprs(RD_gse26366)
#ma.plot(rowMeans(log2(RAWEXP)), log2(RAWEXP[,1])-log2(RAWEXP[,2]), cex=1)
#boxplot(RD_gse26366)


#phenodata import and read
file.exists("GSE26366_pdata.txt")
PD_gse26366<-read.table("GSE26366_pdata.txt", header = TRUE)
#alternate- read.delim("GSE26366_pdata.txt", header = TRUE, fill = TRUE)
#options(max.print = 380)
#alternate through rio- PD_gse26366<-import("GSE26366_pdata.txt")
#head(Phdata)
#for excel file- PD_gse26366<-import("GSE26366_series_matrix.xlsx")
#head(rio_xlsx)

#viewing phenodata = View(Phdata), for dimensions= dim(file name)


#RMA normalization
ND_gse26366 = rma(RD_gse26366)

#viewing normalized data= View(ND_gse26366) 
#boxplot(ND_gse26366@assayData$exprs)


#getting expression value or probe IDs for genes
ED_gse26366<-as.data.frame(exprs(ND_gse26366))


#merging ED and Phdata for allotment of probe IDs to gene IDs
colnames(ED_gse26366) = PD_gse26366$Sample_ID
#to view- colnames(ED_gse26366)
#dim(ED_gse26366)- 22283 206


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse26366=rownames(ED_gse26366)
PIDs_gse26366= getBM(attributes = c("affy_hg_u133a_2", "hgnc_symbol"),
                     filters = "affy_hg_u133a_2",
                     values = PIDs.gse26366,
                     mart = Martfunction)
# alternate- if(interactive()){
#  mart <- useEnsembl(biomart = "ensembl",
#                    dataset = "hsapiens_gene_ensembl")
# 
#  getBM(attributes = c("affy_hg_u133a", "hgnc_symbol", "chromosome_name", "band"),
#        filters    = "affy_hg_u133a",
#        values     = "PIDs.gse26366",
#        mart       = mart)
#}
#dim(PIDs_gse26366) [1] 23467     2
#to view lists- listMarts(),  listDatasets(mart= ensembl), listAttributes(mart = ensembl)


#removing unassigned probe IDs
PIDs_gse26366.df=PIDs_gse26366[!(PIDs_gse26366$hgnc_symbol==""),]
#dim(PIDs_gse26366.df) [1] 22287     2

#bringing rownames to colnames
ED_gse26366.df=rownames_to_column (ED_gse26366, "PIDs")
colnames(PIDs_gse26366.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse26366= merge(PIDs_gse26366.df, ED_gse26366.df, by= "PIDs")
#dim(GS_gse26366)- 22287 208


#limma avereps- to avg out duplicate symbols
AR_gse26366=as.data.frame((limma::avereps(GS_gse26366, GS_gse26366$hgnc_symbol)))
#dim(AR_gse26366)- 14118 207


#final dataset - having umique gene symbols as rownames
AR_gse26366=column_to_rownames(AR_gse26366, "hgnc_symbol")
#dim(AR_gse26366)- 14118 207

#removing PIDs column
FD_gse26366=AR_gse26366[-c(1)]
#dim(FD_gse26366)- 14118 206

#creating matrix file of final dataset
FD_gse26366.matrix=as.matrix(FD_gse26366)
class(FD_gse26366.matrix)= "numeric"

##to create common gene symbol column for merging data
MD_gse26366 <- rownames_to_column (FD_gse26366, "gene.symbol")
#dim(MD_gse26366)- 14118 207


#DGE
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

#saving the dge file
write.csv(DGE_gse26366, file = "1DGE_26366.csv")
write.csv(DGE2_gse26366, file = "1DGE2_26366.csv")

#saving final final for meta-analysis
saveRDS(MD_gse26366, file="1MD_gse26366")