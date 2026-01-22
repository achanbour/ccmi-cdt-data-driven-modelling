# We use this script to install the required packages and download the data-sets required for going through the OSCA books.

## installing relevant packages

# this is for using bioconductor packages
install.packages("BiocManager")

# Eachbook is a biconductor package the following commands install all the relevant packages we need for OSCA book

BiocManager::install("OSCA.intro")
BiocManager::install("OSCA.basics")
BiocManager::install("OSCA.advanced")

## Loading the relevant data-sets in SingleCellExperiments objects

# Loading sce.zeisel data from workflow 2
library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))

library(org.Mm.eg.db)
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db, 
                                      keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")

stats <- perCellQCMetrics(sce.zeisel, subsets=list(
  Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)



# pbmc data loading, filtering, dimensional reduction from workflow 3 

library(DropletTestFiles)
raw.path <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.tar.gz")
out.path <- file.path(tempdir(), "pbmc4k")
untar(raw.path, exdir=out.path)

library(DropletUtils)
fname <- file.path(out.path, "raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

unfiltered <- sce.pbmc

stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]
summary(high.mito)


library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

summary(sizeFactors(sce.pbmc))

set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)


set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")


# 416b data loading from workflow 1

library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                   keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                    keytype="GENEID", column="SEQNAME")

library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
                                           rowData(sce.416b)$SYMBOL)

unfiltered <- sce.416b

mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)


dec.416b <- modelGeneVar(sce.416b)
dec.416b[order(dec.416b$bio, decreasing=TRUE),]

# Load data sce.nest data from workflow 10, this is only for trajectory inference, if you get there uncomment and run

# library(scRNAseq)
# sce.nest <- NestorowaHSCData()
# 
# library(AnnotationHub)
# ens.mm.v97 <- AnnotationHub()[["AH73905"]]
# anno <- select(ens.mm.v97, keys=rownames(sce.nest), 
#                keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
# rowData(sce.nest) <- anno[match(rownames(sce.nest), anno$GENEID),]
# 
# library(scater)
# stats <- perCellQCMetrics(sce.nest)
# qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent")
# sce.nest <- sce.nest[,!qc$discard]
# 
# library(scran)
# set.seed(101000110)
# clusters <- quickCluster(sce.nest)
# sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
# sce.nest <- logNormCounts(sce.nest)
# 
# set.seed(00010101)
# dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
# top.nest <- getTopHVGs(dec.nest, prop=0.1)
# 
# set.seed(101010011)
# sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
# sce.nest <- runTSNE(sce.nest, dimred="PCA")
# 
# snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
# colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

