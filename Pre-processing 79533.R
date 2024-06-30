#set wd- session-choose directory-dekstop-wd

#cel files download and unpacking
gse79533<-list.celfiles("GSE79533/", pattern="CEL")


#reading raw data
RD_gse79533<-ReadAffy(verbose = TRUE, filenames = gse79533)

#visualization of raw data before normalization
#hist(RD_gse79533)
#RAWEXP<-exprs(RD_gse79533)
#ma.plot(rowMeans(log2(RAWEXP)), log2(RAWEXP[,1])-log2(RAWEXP[,2]), cex=1)
boxplot(RD_gse79533)


#phenodata import and read
file.exists("GSE79533_pdata.txt")
PD_gse79533<-read.table("GSE79533_pdata.txt", header = TRUE, sep = "\t")
#alternate- read.delim("GSE79533_pdata.txt", header = TRUE, fill = TRUE)
#options(max.print = 380)
#alternate through rio- PD_gse79533<-import("GSE79533_pdata.txt")
#head(Phdata)
#for excel file- PD_gse79533<-import("GSE79533_series_matrix.xlsx")
#head(rio_xlsx)

#viewing phenodata = View(Phdata), for dimensions= dim(file name)


#RMA normalization
ND_gse79533 = rma(RD_gse79533)

#viewing normalized data= View(ND_gse79533) 
boxplot(ND_gse79533@assayData$exprs)


#getting expression value or probe IDs for genes
ED_gse79533<-as.data.frame(exprs(ND_gse79533))


#merging ED and Phdata for allotment of probe IDs to gene IDs
colnames(ED_gse79533) = PD_gse79533$Sample_IDs
#to view- colnames(ED_gse79533)
#dim(ED_gse79533)- 22283  83


#gene annotation through biomart hgnc
Martfunction=useMart("ENSEMBL_MART_ENSEMBL")
Martfunction=useDataset("hsapiens_gene_ensembl",Martfunction)
PIDs.gse79533=rownames(ED_gse79533)
PIDs_gse79533= getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
                     filters = "affy_hg_u133_plus_2",
                     values = PIDs.gse79533,
                     mart = Martfunction)
# alternate- if(interactive()){
#  mart <- useEnsembl(biomart = "ensembl",
#                    dataset = "hsapiens_gene_ensembl")
# 
#  getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "chromosome_name", "band"),
#        filters    = "affy_hg_u133a",
#        values     = "PIDs.gse79533",
#        mart       = mart)
#}
#dim(PIDs_gse79533) [1] 48315     2
#to view lists- listMarts(),  listDatasets(mart= ensembl), listAttributes(mart = ensembl)


#removing unassigned probe IDs
PIDs_gse79533.df=PIDs_gse79533[!(PIDs_gse79533$hgnc_symbol==""),]
#dim(PIDs_gse79533.df) [1] 22279     2

#bringing rownames to colnames
ED_gse79533.df=rownames_to_column (ED_gse79533, "PIDs")
colnames(PIDs_gse79533.df)[1]<-"PIDs"


#merging gene symbols with ED
GS_gse79533= merge(PIDs_gse79533.df, ED_gse79533.df, by= "PIDs")
#dim(GS_gse79533)- 44137 231


#limma avereps- to avg out duplicate symbols
AR_gse79533=as.data.frame((limma::avereps(GS_gse79533, GS_gse79533$hgnc_symbol)))
#dim(AR_gse79533)- 22269 231


#final dataset - having unique gene symbols as rownames
AR_gse79533=column_to_rownames(AR_gse79533, "hgnc_symbol")
#dim(AR_gse79533)- 22269 230

#removing PIDs column
FD_gse79533=AR_gse79533[-c(1)]
#dim(FD_gse79533)- 22269 229

#converting data frame to matrix for numeric values
FD_gse79533.matrix=as.matrix(FD_gse79533)
class(FD_gse79533.matrix)= "numeric"
#dim(FD_gse79533.matrix) 22269  229

#write.csv(FD_gse79533, file = "FD_gse79533.csv") - to export files in xcel.
##to create common gene symbol column for merging data
MD_gse79533 <- rownames_to_column (FD_gse79533, "gene.symbol")
#dim(MD_gse79533)- 22269 230


#DGE
f.source=factor(PD_gse79533$mutation, levels = c("control", "case"))
design_79533 <- model.matrix(~ 0+factor(f.source))
colnames(design_79533) <-c("control", "case")
fit_79533=lmFit(FD_gse79533.matrix, design_79533)

contrast.matrix <- makeContrasts(case-control, levels = design_79533)
fit2 <- contrasts.fit(fit_79533, contrast.matrix)
fit2 <- eBayes(fit2)

#making top table of DEGs (top 10 and all)
DGE_gse79533<-topTable(fit2, coef=1, adjust="BH", confint=0.95)
DGE2_gse79533<-topTable(fit2, n=Inf, adjust="BH", confint=0.95)

#plotting DEGs
volcanoplot(fit2)

#adding differential expression column in table
DGE2_gse79533$DE = "NO"
DGE2_gse79533$DE[DGE2_gse79533$logFC>0.5 & DGE2_gse79533$adj.P.Val<0.05] = "UP"
DGE2_gse79533$DE[DGE2_gse79533$logFC< -0.5 & DGE2_gse79533$adj.P.Val<0.05] = "DOWN"

#adding top 10 DEG names in separate column in table
DGE2_gse79533=rownames_to_column(DGE2_gse79533, var = "symbol")
DGE2_gse79533$GeneSymbol <-ifelse(DGE2_gse79533$symbol %in% head(DGE2_gse79533[order(DGE2_gse79533$adj.P.Val), 
                                                                               "symbol"], 10), DGE2_gse79533$symbol, NA)
#ggplot for 30 top expressed genes
ggplot(data = DGE2_gse79533, aes(x = logFC, y = -log10(P.Value), col = DE, label=GeneSymbol)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = "dashed") +
scale_color_manual(values = c("blue", "grey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_text_repel(max.overlaps = Inf) #for gene names


#test
#volcano plot of the DEgs
#highlighting the top 25 DEGs
volcanoplot(fit2,coef=1,highlight=25, names=rownames(FD_gse79533.matrix))
dev.off()

#factorising the samples and the condition columns
sc2$Sample<- as.factor(PD_gse79533$Sample_IDs)
sc2$condition<- as.factor(PD_gse79533$mutation)


#plotting a heatmap with  a defined colour scheme
heatmap3(data.matrix.genes2, labCol =paste( sc2$Sample, sc2$condition, sep = " - " ), col=colorRampPalette(c('navy','white','firebrick2'))(1024))

heatmap3(FD_gse79533.matrix, labCol = paste(sc2$Sample, sc2$condition, sep = " - "))
heatmap3(DGE_gse79533, Rowv = DGE2_gse79533$symbol)

#saving metafiles
saveRDS(MD_gse12995, file="MD_gse12995")
saveRDS(MD_gse13425, file="MD_gse13425")
saveRDS(MD_gse79533, file="1MD_gse79533")
saveRDS(MD_gse26281, file="MD_gse26281")
saveRDS(MD_gse26366, file="MD_gse26366")

#saving dge file
write.csv(DGE_gse79533, file = "1DGE_79533.csv")
write.csv(DGE2_gse79533, file = "1DGE2_79533.csv")
