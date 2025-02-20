#This code was created with the help of Jenny Drnevich

setwd("C:/Users/shili/Desktop/MOV10")
options(stringsAsFactors=F)

library(tximport)
library(rtracklayer)
library(Glimma)
library(dplyr)
library(edgeR)
library(RUVSeq)
library(gplots)
library(WGCNA)
library(org.Mm.eg.db)
library(KEGGREST)

#Importing salmon data
samples <- read.csv(file.path("samples.csv"))
samples$condition <- factor(samples$condition)
files <- file.path("data",samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

tx.all <- tximport(files, type="salmon", txOut=TRUE)

#Create in pre-made transcript - gene matching file
gtf0 <- import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.primary_assembly.annotation.gtf.gz")
gtf0 <- gtf0[gtf0$type=="exon"]
gtf1 <- gtf0[!duplicated(gtf0$transcript_id)]
gtf1b <- gtf1[!is.na(gtf1$transcript_id)]
head(gtf1b)

#use metadata columns transcript_id, gene_id, and gene_name
all_tx_info <- mcols(gtf1b)[,c("transcript_id", "gene_id", "gene_name")]

#Get Entrez Gene IDs
con <- gzcon(url("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.metadata.EntrezGene.gz"))
txt <- readLines(con)
Txid2Entrezid <- read.delim(textConnection(txt), header = FALSE)

names(Txid2Entrezid) <- c("transcript_id", "ENTREZID")

#joining together and keep only one entrezid
all_tx_info <- as.data.frame(all_tx_info) %>% left_join(Txid2Entrezid) %>% group_by(transcript_id) %>%
  summarise(gene_id = unique(gene_id), ENTREZID = as.character(min(as.numeric(ENTREZID))))

gene_info <- mcols(gtf1b)[, c("gene_id", "gene_name")][!duplicated(mcols(gtf1b)$gene_id),] %>% 
  as.data.frame()%>%
  left_join(all_tx_info[!duplicated(all_tx_info$gene_id),c("gene_id", "ENTREZID")])

#get gene-level counts
gene.all <- summarizeToGene(tx.all, all_tx_info[,1:2], countsFromAbundance = "lengthScaledTPM")

#make geneID and genename table

gene_info <- left_join(data.frame(gene_id = rownames(gene.all$counts)), gene_info)

rownames(gene_info) <- gene_info$gene_id


#Make DGEList

d <- DGEList(gene.all$counts, samples = samples, genes = gene_info)


#Do TMM norm ----

d <- calcNormFactors(d)


logCPM <- cpm(d, log=T)

glMDSPlot(logCPM, top = 5000, labels = samples$sample_id, groups = samples[,c("condition")],
          html = "MDSclustering_preFiltering", folder = "/stats/interactive_plots")


####Filtering####
min_cpm <- 0.25
min_samp <- 3

i.filter <- rowSums(logCPM > log2(min_cpm)) >= min_samp


#cut down to filtered gene list

d.filt <- d[i.filter, , keep.lib.sizes=F]


#re-do TMM factors
d.filt <- calcNormFactors(d.filt)

#do final normalization----
logCPM.filt <- cpm(d.filt, log = T, prior.count =3)

#Do glMDSPlot

temp <- glMDSPlot(logCPM.filt, top = 5000, labels = samples$sample_id, groups = samples[,c("condition")],
                  html = "MDSclustering_postFiltering", folder = "stats/interactive_plots/", launch = FALSE)
plot.values <- as.data.frame(temp$points)
colnames(plot.values) <- paste0("Dim", 1:ncol(plot.values))
plot.values$var.exp <- round(temp$eig*100/sum(temp$eig),1)
plot.values$Treatment <- factor(d.filt$samples$condition)
plot.values$Label <- factor(d.filt$samples$sample_id)

####Statistical analysis----

e <- new("EList", list(E =logCPM.filt, genes=d.filt$genes, targets=samples))

#make design matrix
Group <- factor(samples$condition)
design0 <- model.matrix(~0+Group)
colnames(design0) <- c("Del", "WT")

#fit model using edgeR, without sv ----
fit0 <- lmFit(e,design0)

cont.matrix0 <- makeContrasts(DelvsWT.0 = Del-WT,
                              levels = design0)
fit1 <- contrasts.fit(fit0, cont.matrix0) %>% eBayes(trend = TRUE)

# Get oneway anova to select most non-significant genes for RUV

res.temp <- topTable(fit1, coef = 1, n = Inf, sort.by = "none")

diff <- makeGroups(Group)

#find negative control genes, use least 5000 significant genes
genes <- res.temp %>% arrange(P.Value) %>% tail(n = 5000) %>% pull(gene_id)

#calculate surrogate variables
wobj <- RUVs(logCPM.filt, genes, k = 2, diff, isLog = TRUE)

no.w <- removeBatchEffect(logCPM.filt, covariates= wobj$W)

glMDSPlot(no.w, top = 5000, labels = d.filt$samples$sample_id,
          groups = samples[,c("condition","rep")],
          html = "MDSclustering_RUV", folder = "/stats/interactive_plots/")

temp <- glMDSPlot(no.w, top = 5000, labels = d.filt$samples$sample_id, groups = samples[,c("condition")],
                  html = "MDSclustering_RUV", folder = "stats/interactive_plots/", launch = F)
plot.values2 <- as.data.frame(temp$points)
colnames(plot.values2) <- paste0("Dim", 1:ncol(plot.values2))
plot.values2$var.exp <- round(temp$eig*100/sum(temp$eig),1)
plot.values2$Treatment <- factor(d.filt$samples$condition)
plot.values2$Label <- factor(d.filt$samples$sample_id)


#add RUV to design matrix

design.sv <- cbind(design0, wobj$W)

fit2 <- lmFit(e, design.sv)

cont.matrix <- makeContrasts(Del_vs_WT = Del-WT,
                             levels = design.sv)

fit3 <- eBayes(contrasts.fit(fit2, cont.matrix), trend = TRUE)

#Pull out results using global FDR correction

source("summarizeFit.R")
results.out <- summarizeFit(fit3, method = "global")
#get oneway anova test to select genes for heatmap

anova.out <- topTable(fit3, coef = 1, num = Inf, sort.by = "none")

#add to results
results.out <- cbind(anova.out[,c(1:3,6:8)], results.out[,-c(1:4)])
names(results.out)[4:6] <- paste0(c("Fstat.","rawP.", "FDR."), "ANOVA")

#heatmap ----

sig.index.int <- anova.out$adj.P.Val < 0.05

#put WT, DEL
new.order <- c(1:3,4:6)

col.pan <- colorpanel(100, "blue", "white", "red")
heat.h0 <- t(scale(t(no.w[sig.index.int,new.order])))
x.cluster.h0 <- hclust(dist(heat.h0))
rowCols.h0 <- WGCNA::cutreeStaticColor(x.cluster.h0, cutHeight = 5, minSize = 1 )
heatVal.h0 <- data.frame(gene_id=rownames(heat.h0), 
                         heat.h0,
                         ModuleColor=rowCols.h0)
heatVal.h0 <- heatVal.h0[rev(x.cluster.h0$order),]
heatVal.h0$TopToBottom <- 1:nrow(heatVal.h0)
results.out  <-  results.out %>% left_join(heatVal.h0 %>% dplyr::select(gene_id, ModuleColor, TopToBottom))
#Get GO and updated KEGG terms

source("getGO.R")
source("getPath.R")

#get unique egids without NAs

egids <- results.out$ENTREZID[!is.na(results.out$ENTREZID)]
egids <- egids[!duplicated(egids)]

tempGO <- getGO(org.Mm.eg.db, keys = egids,keytype = "ENTREZID")

results.out <- results.out %>% left_join(tempGO%>%tibble::rownames_to_column("ENTREZID"))

tempKEGG <- getPath(keys = egids, updatePath = TRUE, species = "mmu")

results.out  <-results.out %>% left_join(data.frame(ENTREZID = rownames(tempKEGG), KEGGpathway = tempKEGG))

#re-order
results.out <- results.out %>% dplyr::select(gene_id:ENTREZID, everything())

#Add in columns showing which genes are in the 513 gene list or in the
#BrainClip list

#read in brain clip list ----

clip930 <- read.csv("BrainClip_680genes.txt")

clip930$ENTREZID <- as.character(clip930$ENTREZID)
clip930 <- clip930[!duplicated(clip930$ENTREZID),]

#add in column to results with TRUE/FALSE for each list

results.out$Clip_930 <- results.out$ENTREZID %in% clip930


#Do venn diagram - add extra lane on codes----

venn.codes <- cbind(decideTests(fit3, method = "global") ,Clip_930 = as.numeric(results.out$Clip_930))

#have to do a bit of manipulation to check cross-over of down-regulated genes. Make new object and
#change the Mov10 and Clip values to -1s

venn.codes.down <- venn.codes
venn.codes.down[,2] <- venn.codes.down[,2] * -1

#Do heatmaps of these genes----

heat.h2 <- t(scale(t(no.w[results.out$Clip_930,new.order])))
x.cluster.h2 <- hclust(dist(heat.h2))
rowCols.h2 <- WGCNA::cutreeStaticColor(x.cluster.h2, cutHeight = 5, minSize = 1 )

heatVal.h2 <- data.frame(gene_id=rownames(heat.h2), 
                         heat.h2,
                         ModCol_Clip930=rowCols.h2)
heatVal.h2 <- heatVal.h2[rev(x.cluster.h2$order),]
heatVal.h2$TopToBot_Clip930<- 1:nrow(heatVal.h2)

results.out  <-  results.out %>% left_join(heatVal.h2 %>% dplyr::select(gene_id, ModCol_Clip930, TopToBot_Clip930))

#write out results
write.table(results.out, 
            file = "stats/results_Mouse_RUV_2021-02-26.csv",
            row.names = FALSE, sep = ",")

cpm.out <- cbind(results.out[,1:5], logCPM.filt[,new.order])
RUV.out <- cbind(results.out[,1:5], no.w[,new.order])

write.table(cpm.out, file = "stats/logCPMvalues_2021-02-26.csv",
            row.names = FALSE, sep = ",")
write.table(RUV.out, file = "stats/logCPMValues_RUVcorrected_2021-02-26.csv",
            row.names = FALSE, sep = ",")
#Do interactive results plots 

res.codes <- decideTests(fit3, method = "global")
md_anno <- fit3$genes
idx <- new.order #change order if necessary
rawPs <-  results.out %>% dplyr::select(contains("rawP"))
rawPs <- rawPs[,-1]  #remove anova
globalFDRs <- results.out %>% dplyr::select(contains("FDR"))
#globalFDRs <- globalFDRs[,-1]

#Group2 <- factor(Group, levels = c("Del","WT"))
Group2 <- factor(c("Del","Del","Del","WT","WT","WT"))
md_anno$rawP <- rawPs
md_anno$global_FDR <- globalFDRs[,1]
glMDPlot(fit3, coef = 1, counts = no.w[,idx],
         groups = Group2[idx],
         samples = samples$sample_id[idx],
         status = res.codes[,1], anno = md_anno,
         main = colnames(res.codes)[1], #display.columns = c("gene_id", "rawP", "global_FDR"),
         sample.cols = idx, folder = "/stats/interactive_results-RUV/",
         html = colnames(res.codes)[1], launch = FALSE)

heatmap.2(heat.h0, col=col.pan, Rowv=as.dendrogram(x.cluster.h0),Colv=FALSE, scale="none",
          trace="none", dendrogram="row", cexRow=0.9, cexCol=1.1, labRow = "",
          keysize = 0.2,lwid = c(1.4,4), lhei = c(0.8,4), key.par=list(mar = c(4,2,2,1)),
          margins = c(7,2),density.info = "none", key.xlab = "SD from mean",
          colsep =c(3,6),
          rowsep,
          sepcolor="black",
          sepwidth=c(0.05,0.05),
          main = paste0("\n\n",nrow(heat.h0), " genes\n", "ANOVA\nFDR<0.05"),
          RowSideColors = rowCols.h0)

temp1 <- decideTests(fit3, method = "global")
#temp1[,1:5] <- trend_global_0.1[,1:5]
write.csv(summary(temp1),"stats/DEG.csv")

vennDiagram(venn.codes[,c(1,2)], include = "up", main = "UP regulated genes", cex=c(1.2,1,0.7))
vennDiagram(venn.codes.down[,c(1,2)], include = "down", main = "DOWN regulated genes", cex=c(1.2,1,0.7))

save.image("Analysis.RData")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")