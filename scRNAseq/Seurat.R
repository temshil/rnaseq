setwd("C:/Users/shili/Desktop/git/scRNAseq")
options(stringsAsFactors=F)

library(Seurat)
library(limma)

# Load the PBMC dataset
pbmc.data <- Read10X("zhong2020")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "zhong2020", min.cells = 30)

#QC
#Calculating percent of mitochondrial RNAs
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# We filter cells that have unique feature counts over 7000 or less than 800
# We filter cells that have >15% mitochondrial counts

pbmc <- subset(pbmc, subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 15)

#Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Finding most variable genes for PCA
# pbmc <- FindVariableFeatures(pbmc, selection.method = "mvp",
#                              mean.cutoff = c(0.0125,8), dispersion.cutoff = c(2,Inf))

pbmc <- FindVariableFeatures(pbmc)

#Scaling the data
memory.limit(50000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


#Determining dimensionality
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc, dims = 1:20)

#Clustering cells.
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.15) #1 --> 28

#t-SNE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20, dim.embed = 2)

DimPlot(pbmc, reduction = "tsne")

FeaturePlot(pbmc, features = c("ASCL1", "NEUROD2", "GAD1", "PROX1", "MEIS2", "LHX6", 
                               "NR2F2","AQP4","OLIG2","MBP","PTPRC", "SPARC"),
            reduction = "tsne")

#Subletting excitatory neurons
exn <- subset(pbmc, idents = c(1,2))

#Sub-clustering cells
exn <- FindNeighbors(exn, dims = 1:20)
exn <- FindClusters(exn, resolution = 0.5)

#t-SNE
exn <- RunTSNE(exn, reduction = "pca", dims = 1:20, dim.embed = 2)

DimPlot(exn)


FeaturePlot(exn, features = c("MEIS2", "MPPED1", "MGAT4C","DPF3", 
                              "PROX1", "SCGN","SEMA5A", "PID1",
                              "SULF2","NRIP3", "NEUROD6"), 
            reduction = "tsne")


#Selecting all the genes expressed by excitatory neurons
exn_exp <- AverageExpression(exn)

#colnames(exn_exp$RNA) <- c("CA2","CA1","CA3", "DG", "4", "5", "6")

colnames(exn_exp$RNA) <- c("CA2","CA1","DG", "CA3", "4", "5", "6")

#Selecting genes of each cluster
CA1_h <- rownames(subset(exn_exp$RNA, exn_exp$RNA[,"CA1"]>0))
CA2_h <- rownames(subset(exn_exp$RNA, exn_exp$RNA[,"CA2"]>0))
CA3_h <- rownames(subset(exn_exp$RNA, exn_exp$RNA[,"CA3"]>0))
DG_h <- rownames(subset(exn_exp$RNA, exn_exp$RNA[,"DG"]>0))

INH_h <- AverageExpression(subset(pbmc, idents = c(1)))
INH_h <- rownames(subset(INH_h$RNA, INH_h$RNA[,"all"]>0))

CA1_hp <- nrow(data.frame(CA1_h))/17787*100
CA2_hp <- nrow(data.frame(CA2_h))/17787*100
CA3_hp <- nrow(data.frame(CA3_h))/17787*100
DG_hp <- nrow(data.frame(DG_h))/17787*100
INH_hp <- nrow(data.frame(INH_h))/17787*100

convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  mousex <- unique(genesV2[, 2])
  
  return(mousex)
}

#Converting human genes to mouse genes
CA1_m <- data.frame("genes"=convertHumanGeneList(CA1_h))
CA2_m <- data.frame("genes"=convertHumanGeneList(CA2_h))
CA3_m <- data.frame("genes"=convertHumanGeneList(CA3_h))
DG_m <- data.frame("genes"=convertHumanGeneList(DG_h))
INH_m <- data.frame("genes"=convertHumanGeneList(INH_h))

CA1_mp <- nrow(CA1_m)/nrow(data.frame(CA1_h))*100
CA2_mp <- nrow(CA1_m)/nrow(data.frame(CA2_h))*100
CA3_mp <- nrow(CA1_m)/nrow(data.frame(CA3_h))*100
DG_mp <- nrow(CA1_m)/nrow(data.frame(DG_h))*100
INH_mp <- nrow(CA1_m)/nrow(data.frame(INH_h))*100


Up <- read.table("Up.txt", col.names = "Genes")
Down <- read.table("Down.txt", col.names = "Genes")


CA1_perc_up <- sum(as.numeric(CA1_m$genes %in% Up$Genes))/nrow(Up)*100
CA1_perc_down <- sum(as.numeric(CA1_m$genes %in% Down$Genes))/nrow(Down)*100

CA2_perc_up <- sum(as.numeric(CA2_m$genes %in% Up$Genes))/nrow(Up)*100
CA2_perc_down <- sum(as.numeric(CA2_m$genes %in% Down$Genes))/nrow(Down)*100

CA3_perc_up <- sum(as.numeric(CA3_m$genes %in% Up$Genes))/nrow(Up)*100
CA3_perc_down <- sum(as.numeric(CA3_m$genes %in% Down$Genes))/nrow(Down)*100

DG_perc_up <- sum(as.numeric(DG_m$genes %in% Up$Genes))/nrow(Up)*100
DG_perc_down <- sum(as.numeric(DG_m$genes %in% Down$Genes))/nrow(Down)*100

INH_perc_up <- sum(as.numeric(INH_m$genes %in% Up$Genes))/nrow(Up)*100
INH_perc_down <- sum(as.numeric(INH_m$genes %in% Down$Genes))/nrow(Down)*100



CA1_mar <- data.frame("genes"=rownames(FindMarkers(exn, ident.1 = 1, min.pct = 0.01)))
CA2_mar <- data.frame("genes"=rownames(FindMarkers(exn, ident.1 = 0, min.pct = 0.01)))
CA3_mar <- data.frame("genes"=rownames(FindMarkers(exn, ident.1 = 3, min.pct = 0.01)))
DG_mar <- data.frame("genes"=rownames(FindMarkers(exn, ident.1 = 2, min.pct = 0.01)))
INH_mar <- data.frame("genes"=rownames(FindMarkers(pbmc, ident.1 = 0, min.pct = 0.01)))


AllMar <- unique.data.frame(rbind(CA1_mar,CA2_mar,CA3_mar, DG_mar, INH_mar))
AllMar$CA1 <- as.numeric(AllMar$genes %in% CA1_mar$genes)
AllMar$CA2 <- as.numeric(AllMar$genes %in% CA2_mar$genes)
AllMar$CA3 <- as.numeric(AllMar$genes %in% CA3_mar$genes)
AllMar$DG <- as.numeric(AllMar$genes %in% DG_mar$genes)
AllMar$INH <- as.numeric(AllMar$genes %in% INH_mar$genes)
vennDiagram(AllMar[,c(2,3,4,5,6)])


CA1_marm <- data.frame("genes"=convertHumanGeneList(CA1_mar$genes))
CA2_marm <- data.frame("genes"=convertHumanGeneList(CA2_mar$genes))
CA3_marm <- data.frame("genes"=convertHumanGeneList(CA3_mar$genes))
DG_marm <- data.frame("genes"=convertHumanGeneList(DG_mar$genes))
INH_marm <- data.frame("genes"=convertHumanGeneList(INH_mar$genes))


CA1_perc_up_mar <- sum(as.numeric(CA1_marm$genes %in% Up$Genes))/nrow(Up)*100
CA1_perc_down_mar <- sum(as.numeric(CA1_marm$genes %in% Down$Genes))/nrow(Down)*100

CA2_perc_up_mar <- sum(as.numeric(CA2_marm$genes %in% Up$Genes))/nrow(Up)*100
CA2_perc_down_mar <- sum(as.numeric(CA2_marm$genes %in% Down$Genes))/nrow(Down)*100

CA3_perc_up_mar <- sum(as.numeric(CA3_marm$genes %in% Up$Genes))/nrow(Up)*100
CA3_perc_down_mar <- sum(as.numeric(CA3_marm$genes %in% Down$Genes))/nrow(Down)*100

DG_perc_up_mar <- sum(as.numeric(DG_marm$genes %in% Up$Genes))/nrow(Up)*100
DG_perc_down_mar <- sum(as.numeric(DG_marm$genes %in% Down$Genes))/nrow(Down)*100

INH_perc_up_mar <- sum(as.numeric(INH_marm$genes %in% Up$Genes))/nrow(Up)*100
INH_perc_down_mar <- sum(as.numeric(INH_marm$genes %in% Down$Genes))/nrow(Down)*100


#Clustering cells.
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 1) #1 --> 28

#t-SNE
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:20, dim.embed = 2)

DimPlot(pbmc, reduction = "tsne", label = TRUE)

FeaturePlot(pbmc, features = c("MEIS2", "MPPED1", "MGAT4C","DPF3", 
                              "PROX1", "SCGN","SEMA5A", "PID1",
                              "SULF2","NRIP3", 'NEUROD6'), 
            reduction = "tsne")



ac <- data.frame("clus"=Idents(pbmc))
ca1 <- subset(ac, ac$"clus"==3)
ca1 <- data.frame("cells" = rownames(ca1))

ca2 <- subset(ac, ac$"clus"==2)
ca2 <- data.frame("cells" = rownames(ca2))

ca3 <- subset(ac, ac$"clus"==4)
ca3 <- data.frame("cells" = rownames(ca3))

dg <- subset(ac, ac$"clus"==6)
dg <- data.frame("cells" = rownames(dg))

acex <- data.frame("clus"=Idents(exn))
ca1ex <- subset(acex, acex$"clus"==1)
ca1ex <- data.frame("cells" = rownames(ca1ex))

ca2ex <- subset(acex, acex$"clus"==0)
ca2ex <- data.frame("cells" = rownames(ca2ex))

ca3ex <- subset(acex, acex$"clus"==2)
ca3ex <- data.frame("cells" = rownames(ca3ex))

dgex <- subset(acex, acex$"clus"==3)
dgex <- data.frame("cells" = rownames(dgex))

nrow(ca1)
nrow(ca1ex)
ca1p <- sum(as.numeric(ca1$cells %in% ca1ex$cells))/nrow(ca1)*100

nrow(ca2)
nrow(ca2ex)
ca2p <- sum(as.numeric(ca2$cells %in% ca2ex$cells))/nrow(ca2)*100

nrow(ca3)
nrow(ca3ex)
ca3p <- sum(as.numeric(ca3$cells %in% ca3ex$cells))/nrow(ca3)*100

nrow(dg)
nrow(dgex)
dgp <- sum(as.numeric(dg$cells %in% dgex$cells))/nrow(dg)*100

ca1mar <- data.frame("names"=rownames(FindMarkers(pbmc, ident.1 = 3, min.pct = 0.01)))
ca2mar <- data.frame("names"=rownames(FindMarkers(pbmc, ident.1 = 2, min.pct = 0.01)))
ca3mar <- data.frame("names"=rownames(FindMarkers(pbmc, ident.1 = 4, min.pct = 0.01)))
dgmar <- data.frame("names"=rownames(FindMarkers(pbmc, ident.1 = 6, min.pct = 0.01)))
inhmar <- data.frame("names"=rownames(FindMarkers(pbmc, ident.1 = 0, min.pct = 0.01)))

nrow(CA1_mar)
nrow(ca1mar)
ca1mar_p <- sum(as.numeric(CA1_mar$genes %in% ca1mar$names))/nrow(ca1mar)*100

nrow(CA2_mar)
nrow(ca2mar)
ca2mar_p <- sum(as.numeric(CA2_mar$genes %in% ca2mar$names))/nrow(ca2mar)*100

nrow(CA3_mar)
nrow(ca3mar)
ca3mar_p <- sum(as.numeric(CA3_mar$genes %in% ca3mar$names))/nrow(ca3mar)*100

nrow(DG_mar)
nrow(dgmar)
dgmar_p <- sum(as.numeric(DG_mar$genes %in% dgmar$names))/nrow(dgmar)*100

allmar <- unique.data.frame(rbind(ca1mar,ca2mar,ca3mar, dgmar, inhmar))
allmar$ca1 <- as.numeric(allmar$names %in% ca1mar$names)
allmar$ca2 <- as.numeric(allmar$names %in% ca2mar$names)
allmar$ca3 <- as.numeric(allmar$names %in% ca3mar$names)
allmar$dg <- as.numeric(allmar$names %in% dgmar$names)
allmar$inh <- as.numeric(allmar$names %in% inhmar$names)
vennDiagram(allmar[,c(2,3,4,5,6)])

ca1mar_m <- data.frame("names"=convertHumanGeneList(ca1mar$names))
ca2mar_m <- data.frame("names"=convertHumanGeneList(ca2mar$names))
ca3mar_m <- data.frame("names"=convertHumanGeneList(ca3mar$names))
dgmar_m <- data.frame("names"=convertHumanGeneList(dgmar$names))
inhmar_m <- data.frame("names"=convertHumanGeneList(inhmar$names))

ca1_perc_up <- sum(as.numeric(Up$Genes %in% ca1mar_m$names))/nrow(Up)*100
ca1_perc_down <- sum(as.numeric(Down$Genes %in% ca1mar_m$names))/nrow(Down)*100

ca2_perc_up <- sum(as.numeric(Up$Genes %in% ca2mar_m$names))/nrow(Up)*100
ca2_perc_down <- sum(as.numeric(Down$Genes %in% ca2mar_m$names))/nrow(Down)*100

ca3_perc_up <- sum(as.numeric(Up$Genes %in% ca3mar_m$names))/nrow(Up)*100
ca3_perc_down <- sum(as.numeric(Down$Genes %in% ca3mar_m$names))/nrow(Down)*100

dg_perc_up <- sum(as.numeric(Up$Genes %in% dgmar_m$names))/nrow(Up)*100
dg_perc_down <- sum(as.numeric(Down$Genes %in% dgmar_m$names))/nrow(Down)*100

inh_perc_up <- sum(as.numeric(Up$Genes %in% inhmar_m$names))/nrow(Up)*100
inh_perc_down <- sum(as.numeric(Down$Genes %in% inhmar_m$names))/nrow(Down)*100
