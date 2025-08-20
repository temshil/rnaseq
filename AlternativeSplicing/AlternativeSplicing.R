setwd("/Users/shilikbay/Desktop")
options(stringsAsFactors=F)

library(IsoformSwitchAnalyzeR)
library("biomaRt")

# if (!requireNamespace("devtools", quietly = TRUE)){
#   install.packages("devtools")
# }
# devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)

salmonQuant <- importIsoformExpression(parentDir = "data")

myDesign <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = gsub('.{1}$', '', colnames(salmonQuant$abundance)[-1])
)

aSwitchList <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "data/gencode.vM26.annotation.gtf.gz",
  isoformNtFasta       = "data/gentrome.fa.gz",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE
)

summary(aSwitchList)

exampleSwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist   = aSwitchList,
  pathToOutput = 'data/output',
  outputSequences      = TRUE, # change to TRUE when analyzing your own data 
  prepareForWebServers = FALSE  # change to TRUE if you will use webservers for external sequence analysis
)

extractSwitchSummary(exampleSwitchList)

mySwitchList <- preFilter( exampleSwitchList )

mySwitchList <- isoformSwitchTestDEXSeq(  mySwitchList, reduceToSwitchingGenes=FALSE )

mySwitchList <- analyzeAlternativeSplicing( mySwitchList )

pdf(file = 'data/output/SplicingSummary.pdf', onefile = FALSE, height=6, width = 22)
extractSplicingSummary(
  mySwitchList,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)
dev.off()

pdf(file = 'data/output/SplicingEnrichment.pdf', onefile = FALSE, height=6, width = 9)
splicingEnrichment <- extractSplicingEnrichment(
  mySwitchList,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)
dev.off()

IsoformSwitching <- extractTopSwitches(mySwitchList, n=91)
isoforms <- data.frame(IsoformSwitching$gene_name)
colnames(isoforms) <- "gene_name"

clip <- read.csv("data/DEG_CLIP.txt")
deg <- read.csv("data/DEG.txt")

all_clip <-  unique.data.frame(rbind(isoforms, clip))
all_deg <-  unique.data.frame(rbind(isoforms, deg))

all_clip$clip <- as.numeric(all_clip$gene_name %in% clip$gene_name)
all_clip$isoforms <- as.numeric(all_clip$gene_name %in% isoforms$gene_name)

all_deg$deg <- as.numeric(all_deg$gene_name %in% deg$gene_name)
all_deg$isoforms <- as.numeric(all_deg$gene_name %in% isoforms$gene_name)

pdf(file = 'data/output/Venn.pdf', onefile = FALSE, height=10, width = 5)
layout(matrix(1:2))
vennDiagram(all_clip[,c(2,3)]) 
vennDiagram(all_deg[,c(2,3)])
dev.off()

clip_all <- read.csv("clip_all.txt")

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mapping <- getBM(
  attributes = c('entrezgene_id', 'mgi_symbol'), 
  filters = "entrezgene_id",
  values = clip_all,
  mart = mouse
)

clip_all <- data.frame("gene_name"=mapping[,2])
all_genes <- unique.data.frame(rbind(clip_all,isoforms))
all_genes$clip_all <- as.numeric(all_genes$gene_name %in% clip_all$gene_name)
all_genes$isoforms <- as.numeric(all_genes$gene_name %in% isoforms$gene_name)
all_genes$clip_deg <- as.numeric(all_genes$gene_name %in% clip$gene_name)

pdf(file = 'Venn_2.pdf', onefile = FALSE, height=5, width = 5)
vennDiagram(all_genes[,c(2,3,4)]) 
dev.off()

subset(all_clip,all_clip$clip>0 & all_clip$isoforms>0)

pdf(file = 'data/output/Acly.pdf', onefile = FALSE, height=10, width = 14)
switchPlot(
  mySwitchList,
  gene='Acly',
  condition1 = 'Del',
  condition2 = 'WT',
  localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
)
dev.off()

pdf(file = 'data/output/Rbm33.pdf', onefile = FALSE, height=10, width = 14)
switchPlot(
  mySwitchList,
  gene='Rbm33',
  condition1 = 'Del',
  condition2 = 'WT',
  localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
)
dev.off()

pdf(file = 'data/output/Adamts2.pdf', onefile = FALSE, height=10, width = 14)
switchPlot(
  mySwitchList,
  gene='Adamts2',
  condition1 = 'Del',
  condition2 = 'WT',
  localTheme = theme_bw(base_size = 13) # making text sightly larger for vignette
)
dev.off()

save.image("data/AlternativeSplicing.RData")
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")