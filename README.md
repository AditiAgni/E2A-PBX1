# E2A-PBX1 Project
The project aimed to identify the genes dysregulated in E2A-PBX1 translocation. The expression microarray datasets from GEO (GSE26281, GSE26366, and GSE79533) on ALL were used to perform WGCNA to identify genes associated with E2A-PBX1. Only common translocations among all three datasets were retained to ensure homogeneity. T-ALL samples and samples with common sample IDs (from GSE26281 and GSE26366) were removed. The updated phenodata is provided along with the scripts. All the datasets were from Affymetrix platform. The raw data and the series matrix file were downloaded and unzipped manually. Raw data for each dataset was filtered according to updated pheno-data and then normalized using RMA function in RStudio. Differential gene expression was performed to visualize the perturbed genes. The datasets were merged to form metafile that was batch-corrected using "comBat" function. This batch-corrected file was utilized to construct a signed network for WGCNA. The modules significantly correlated to E2A-PBX1 were analyzed for correlation within the module and the gene list was extracted according to the order of gene significance.

# Analysis steps
Raw data was downloaded from GEO and the cel files were extracted manually. The series matrix data was edited in excel for each dataset, and the updated phenodata is provided with the scripts. A common library script was created and ran first to load all the packages required for the analysis. After these steps, pre-processing of each dataset was performed individually. The final matrix file was saved and imported again to create metafile and perform WGCNA. The phenodata for metafile was created by merging phenodata from individual datasets and is provided along with the MetaWGCNA script.

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26281
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26366
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79533



# Details of the uploaded files:

GSE26281_updata: Updated phenodata for GSE26281 dataset.

GSE26366_updata: Updated phenodata for GSE26366 dataset.

GSE79533_updata: Updated phenodata for GSE79533 dataset.

Pmeta: excel file for metafile phenodata. contain samples combined from GSE26281, GSE26366, GSE79533 in the respective order. Alongside Batch and Dataset ID information for batch correction.

libraries: Common library of required packages for the complete analysis.

Pre-processing 26281: R script for pre-processing of GSE26281 dataset. 

Pre-processing 26366: R script for pre-processing of GSE26366 dataset.

Pre-processing 79533: R script for pre-processing of GSE79533 dataset.

MetaWGCNA: R script for creation of metafile and weighted gene co-expression network analysis.

Contact: aditiagnihotri39@gmail.com



