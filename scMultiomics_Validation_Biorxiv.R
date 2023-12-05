#### Analysis of scMultiomics in BCP-ALL
options(bedtools.path = "/Users/sdigiove/miniforge3/envs/py36/bin/")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(bedtoolsr)
library(cicero)
library(leiden)
library(dplyr)
library(igraph)
library(tidyverse)
# plan("multicore", workers = 4)
# options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
# memory.limit(16000)
# Annotations
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr' , seqlevels(annotation))
genome(annotation) <- "hg19"

### Create Common Peaks list
peaks.E1 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E1/atac_peaks.bed",
  col.names = c('chr','start','end')
)
peaks.E2 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E2/atac_peaks.bed",
  col.names = c('chr','start','end')
)
peaks.R1 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R1/atac_peaks.bed",
  col.names = c('chr','start','end')
)
peaks.R2 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R2/atac_peaks.bed",
  col.names = c('chr','start','end')
)
# convert to genomic ranges
gr.E1 <- makeGRangesFromDataFrame(peaks.E1)
gr.E2 <- makeGRangesFromDataFrame(peaks.E2)
gr.R1 <- makeGRangesFromDataFrame(peaks.R1)
gr.R2 <- makeGRangesFromDataFrame(peaks.R2)
# Combine peaks
combined.peaks <- reduce(x = c(gr.E1,gr.E2,gr.R1,gr.R2))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

### sample E1
counts.E1 <- Read10X_h5("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E1/filtered_feature_bc_matrix.h5")
fragpath.E1 <- "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E1/atac_fragments.tsv.gz"
barcodes.E1 <- read.delim("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E1/barcodes.tsv" , col.names = c('filtered.cells'))
md.E1 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E1/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]
md.E1 <- md.E1[barcodes.E1$filtered.cells,]

frags.E1 <- CreateFragmentObject(
  path = fragpath.E1,
  cells = rownames(md.E1)
)
E1.counts.atac <- FeatureMatrix(
  fragments = frags.E1,
  features = combined.peaks,
  cells = rownames(md.E1)
)

E1.obj <- CreateSeuratObject(
  counts = counts.E1$`Gene Expression`[,rownames(md.E1)],
  assay = 'RNA',
  meta.data = md.E1
)
E1.obj[['ATAC']] <- CreateChromatinAssay(
  counts = E1.counts.atac,
  sep = c(':','-'),
  fragments = fragpath.E1,
  annotation = annotation
)

### sample E2
counts.E2 <- Read10X_h5("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E2/filtered_feature_bc_matrix.h5")
fragpath.E2 <- "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E2/atac_fragments.tsv.gz"
barcodes.E2 <- read.delim("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E2/barcodes.tsv" , col.names = c('filtered.cells'))
md.E2 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_E2/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]
md.E2 <- md.E2[barcodes.E2$filtered.cells,]

frags.E2 <- CreateFragmentObject(
  path = fragpath.E2,
  cells = rownames(md.E2)
)
E2.counts.atac <- FeatureMatrix(
  fragments = frags.E2,
  features = combined.peaks,
  cells = rownames(md.E2)
)

E2.obj <- CreateSeuratObject(
  counts = counts.E2$`Gene Expression`[,rownames(md.E2)],
  assay = 'RNA',
  meta.data = md.E2
)
E2.obj[['ATAC']] <- CreateChromatinAssay(
  counts = E2.counts.atac,
  sep = c(':','-'),
  fragments = fragpath.E2,
  annotation = annotation
)

### sample R1
counts.R1 <- Read10X_h5("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R1/filtered_feature_bc_matrix.h5")
fragpath.R1 <- "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R1/atac_fragments.tsv.gz"
barcodes.R1 <- read.delim("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R1/barcodes.tsv" , col.names = c('filtered.cells'))
md.R1 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R1/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]
md.R1 <- md.R1[barcodes.R1$filtered.cells,]

frags.R1 <- CreateFragmentObject(
  path = fragpath.R1,
  cells = rownames(md.R1)
)
R1.counts.atac <- FeatureMatrix(
  fragments = frags.R1,
  features = combined.peaks,
  cells = rownames(md.R1)
)

R1.obj <- CreateSeuratObject(
  counts = counts.R1$`Gene Expression`[,rownames(md.R1)],
  assay = 'RNA',
  meta.data = md.R1
)
R1.obj[['ATAC']] <- CreateChromatinAssay(
  counts = R1.counts.atac,
  sep = c(':','-'),
  fragments = fragpath.R1,
  annotation = annotation
)

### sample R2
counts.R2 <- Read10X_h5("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R2/filtered_feature_bc_matrix.h5")
fragpath.R2 <- "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R2/atac_fragments.tsv.gz"
barcodes.R2 <- read.delim("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R2/barcodes.tsv" , col.names = c('filtered.cells'))
md.R2 <- read.table(
  file = "/Users/sdigiove/Documents/Work/Projects/BCP_ALL/data/scMultiomics/sample_R2/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]
md.R2 <- md.R2[barcodes.R2$filtered.cells,]

frags.R2 <- CreateFragmentObject(
  path = fragpath.R2,
  cells = rownames(md.R2)
)
R2.counts.atac <- FeatureMatrix(
  fragments = frags.R2,
  features = combined.peaks,
  cells = rownames(md.R2)
)

R2.obj <- CreateSeuratObject(
  counts = counts.R2$`Gene Expression`[,rownames(md.R2)],
  assay = 'RNA',
  meta.data = md.R2
)
R2.obj[['ATAC']] <- CreateChromatinAssay(
  counts = R2.counts.atac,
  sep = c(':','-'),
  fragments = fragpath.R2,
  annotation = annotation
)


### Merge Objects
E1.obj$dataset <- 'sample_E1'
E2.obj$dataset <- 'sample_E2'
R1.obj$dataset <- 'sample_R1'
R2.obj$dataset <- 'sample_R2'

E1.obj$patient <- 'Patient1'
E2.obj$patient <- 'Patient2'
R1.obj$patient <- 'Patient1'
R2.obj$patient <- 'Patient2'


E1.obj$condition <- 'Malignant'
E2.obj$condition <- 'Malignant'
R1.obj$condition <- 'Relapse'
R2.obj$condition <- 'Relapse'


# Save single Objects merged.
# saveRDS(E1.obj , "/Volumes/SSD1T/singleCellBCPALL/sampleE1Obj.rds")
# saveRDS(E2.obj , "/Volumes/SSD1T/singleCellBCPALL/sampleE2Obj.rds")
# saveRDS(R1.obj , "/Volumes/SSD1T/singleCellBCPALL/sampleR1Obj.rds")
# saveRDS(R2.obj , "/Volumes/SSD1T/singleCellBCPALL/sampleR2Obj.rds")

# E1.obj <- readRDS("/Volumes/SSD1T/singleCellBCPALL/sampleE1Obj.rds")
# E2.obj <- readRDS("/Volumes/SSD1T/singleCellBCPALL/sampleE2Obj.rds")
# R1.obj <- readRDS("/Volumes/SSD1T/singleCellBCPALL/sampleR1Obj.rds")
# R2.obj <- readRDS("/Volumes/SSD1T/singleCellBCPALL/sampleR2Obj.rds")


# merge all datasets, adding cell ID to make sure cell names are unique
combined.samples.obj <- merge(
  x = E1.obj,
  y = list(E2.obj, R1.obj, R2.obj),
  add.cell.ids = c('E1' , 'E2' , 'R1' , 'R2')
)

### Quality Control
# Cell filtering
DefaultAssay(combined.samples.obj) <- "ATAC"
combined.samples.obj <- NucleosomeSignal(object = combined.samples.obj)
combined.samples.obj <- TSSEnrichment(combined.samples.obj, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
combined.samples.obj$pct_reads_in_peaks <- combined.samples.obj$atac_peak_region_fragments / combined.samples.obj$atac_fragments * 100
combined.samples.obj$blacklist_fraction <- FractionCountsInRegion(
  object = combined.samples.obj, 
  assay = 'ATAC',
  regions = blacklist_hg19
)

DefaultAssay(combined.samples.obj) <- "RNA"
combined.samples.obj[["percent.mt"]] <- PercentageFeatureSet(combined.samples.obj, pattern = "^MT-")

# saveRDS(combined.samples.obj , "/Volumes/SSD1T/singleCellBCPALL/AllSamplesObj.rds")
combined.samples.obj <- readRDS("/Volumes/SSD1T/singleCellBCPALL/AllSamplesObj.rds")
### Plot Quality Control
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/qcATAC_divided.png" , height = 5 , width = 9 , units = 'in' , res = 300)
# png("/Users/sdigiove/Downloads/qcATAC_divided.png" , height = 5 , width = 9 , units = 'in' , res = 300)
VlnPlot(
  object = combined.samples.obj,
  # features = c('pct_reads_in_peaks', "atac_peak_region_fragments",
  #              'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  features = c('pct_reads_in_peaks' , "nFeature_RNA", "percent.mt"),
  pt.size = 0,
  ncol = 3,
  flip = TRUE,
  group.by = 'dataset'
)
# dev.off()

# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/qcRNA_divided.png" , height = 5 , width = 9 , units = 'in' , res = 300)
VlnPlot(
  object = combined.samples.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3,
  group.by = 'dataset'
)
# dev.off()

combined.samples.obj.filt <- subset(
  x = combined.samples.obj,
  subset = pct_reads_in_peaks > 50 &
    atac_peak_region_fragments < 200000 &
    TSS.enrichment <20 &
    TSS.enrichment > 1 &
    nucleosome_signal < 2 &
    nucleosome_signal > 0.2 &
    nFeature_RNA < 10000 &
    nCount_RNA < 100000 &
    percent.mt < 20
)

# saveRDS(combined.samples.obj.filt , "/Volumes/SSD1T/singleCellBCPALL/AllSamplesObj_filtered.rds")
# combined.samples.obj.filt <- readRDS("/Volumes/SSD1T/singleCellBCPALL/AllSamplesObj_filtered.rds")

# combined.samples.obj.filt <- subset(combined.samples.obj.filt,
#                                     subset = dataset == 'sample_E1' | dataset == 'sample_R1')

# RNA normalization
DefaultAssay(combined.samples.obj.filt) <- "RNA"
combined.samples.obj.filt <- SCTransform(combined.samples.obj.filt)
all.genes <- rownames(combined.samples.obj.filt)
combined.samples.obj.filt <- ScaleData(combined.samples.obj.filt, features = all.genes)
combined.samples.obj.filt <- RunPCA(combined.samples.obj.filt)

# ATAC normalization
DefaultAssay(combined.samples.obj.filt) <- "ATAC"
combined.samples.obj.filt <- FindTopFeatures(combined.samples.obj.filt, min.cutoff = 5)
combined.samples.obj.filt <- RunTFIDF(combined.samples.obj.filt)
combined.samples.obj.filt <- RunSVD(combined.samples.obj.filt)

# combined.samples.obj.filt <- RunHarmony(combined.samples.obj.filt, lambda =c(1,1)  , c("patient" , "dataset") , ncores=4) #lambda = c(1,1,1)
# RunHarmony(combined.samples.obj.filt, lambda = c(1,1) , c("patient" , 'dataset') , plot_convergence = TRUE)
DefaultAssay(combined.samples.obj.filt) <- "ATAC"
combined.samples.obj.filt <- RunHarmony(
  object = combined.samples.obj.filt,
  lambda = c(1,1),
  group.by.vars =  c("dataset" , "patient") ,
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE,
  reduction.save = "harmony_atac"
)
combined.samples.obj.filt <- RunHarmony(
  object = combined.samples.obj.filt,
  lambda = c(1,1),
  group.by.vars =  c("dataset" , "patient") ,
  reduction = 'pca',
  assay.use = 'RNA',
  project.dim = FALSE,
  reduction.save = "harmony_rna"
)

# DimPlot(combined.samples.obj.filt, reduction = 'harmony' , label = TRUE, repel = TRUE , group.by = 'dataset' )

### Visualization UMAP
# build a joint neighbor graph using both assays
combined.samples.obj.filt <- FindMultiModalNeighbors(
  object = combined.samples.obj.filt,
  reduction.list = list("harmony_rna", "harmony_atac"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
combined.samples.obj.filt <- RunUMAP(
  object = combined.samples.obj.filt,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "umap",
  #reduction = 'harmony',
  verbose = TRUE
)

table(combined.samples.obj.filt@meta.data$dataset) |> as.data.frame() |> 
  ggplot(aes(x = Var1 , y = Freq , fill = Var1))+
  geom_bar(stat = 'identity' , width = 0.5)+
  theme_classic()+
  ylab("Number of Cells")+
  xlab("")+
  theme(legend.position="none")+
  theme(text = element_text(family="", size=20),
        plot.title = element_text(hjust = 0.5, face="bold", size=20),
        axis.title.y = element_text(size=15, margin=margin(r=25)),
        axis.title.x = element_text(size=10, margin=margin(t=25)),
        axis.text = element_text(size=15, color="black"))
# ggsave("/Users/sdigiove/Downloads/NumberOfCellsPerSample.png" , height = 7 , width = 5 , units = 'in' , dpi = 300)
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/UMAP_multiomics.png" , height = 5 , width = 9 , units = 'in' , res = 300)
DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE, reduction = "umap")
# dev.off()
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/UMAP_multiomics_divided.png" , height = 5 , width = 9 , units = 'in' , res = 300)
DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE, reduction = "umap" , group.by = 'dataset')
# dev.off()
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/UMAP_multiomics_divided.png" , height = 5 , width = 9 , units = 'in' , res = 300)
DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE, reduction = "umap" , group.by = 'Phase')
# dev.off()

# FeaturePlot(combined.samples.obj.filt, features = "pct_reads_in_peaks" , reduction = 'umap')

### UMAP of RNA
# combined.samples.obj.filt <- FindNeighbors(combined.samples.obj.filt, dims = 1:50 , reduction = 'pca' , assay = 'RNA')
# combined.samples.obj.filt <- RunUMAP(combined.samples.obj.filt, dims = 1:50 , reduction = 'pca' , assay = 'RNA')
# DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE, reduction = "umap" , group.by = 'dataset')
# 
# combined.samples.obj.filt <- FindNeighbors(combined.samples.obj.filt, dims = 2:40 , reduction = 'lsi' , assay = 'ATAC')
# combined.samples.obj.filt <- RunUMAP(combined.samples.obj.filt, dims = 2:40 , reduction = 'lsi' , assay = 'ATAC')
# DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE, reduction = "umap" , group.by = 'dataset')


# Clustering
DefaultAssay(combined.samples.obj.filt) <- 'RNA'
combined.samples.obj.filt <- FindClusters(object = combined.samples.obj.filt,graph.name = 'wknn' , algorithm = 4 , random.seed = 100,
                                          resolution = c(0.2 , 0.4, 0.6, 0.8, 1.0,1.2, 1.4))

### Evaluate Modularity for each resolution

wknn.graph <- graph_from_adjacency_matrix(combined.samples.obj.filt@graphs$wknn , mode = 'undirected',diag=F)

mod <- c()
clust.num <- c()
for (clust.resolution in c("wknn_res.0.2","wknn_res.0.4","wknn_res.0.6","wknn_res.0.8","wknn_res.1","wknn_res.1.2","wknn_res.1.4")) {
  mod.value <- modularity( wknn.graph , combined.samples.obj.filt@meta.data[[clust.resolution]])
  mod <- c(mod , mod.value)
  clust.num <- c(clust.num , length(levels(combined.samples.obj.filt@meta.data[[clust.resolution]])))
}

mod.plot <- data.frame(list('resolution.seurat' = c(0.2 , 0.4, 0.6, 0.8, 1.0,1.2, 1.4) , 'modularity.eval' = mod)) |> 
  ggplot(aes(x = resolution.seurat , y = modularity.eval ))+
  geom_bar(stat = 'identity') +
  theme_classic()

clust.num.plot <- data.frame(list('resolution.seurat' = c(0.2 , 0.4, 0.6, 0.8, 1.0,1.2, 1.4) , 'clust.num' = clust.num)) |> 
  ggplot(aes(x = resolution.seurat , y = clust.num ))+
  geom_bar(stat = 'identity') +
  theme_classic()

mod.plot + clust.num.plot

### Plot the cluster stability at different resolutino
combined <- combined.samples.obj.filt@meta.data[,c('wknn_res.0.4','wknn_res.0.6','wknn_res.0.8','wknn_res.1','wknn_res.1.4')]
library(clustree)
set.seed(1111)
clustree(combined, prefix="wknn_res.", edge_arrow=FALSE)


#### Find Clusters at the selected resolution
reticulate::use_python("/Users/sdigiove/miniforge3/envs/leidenalg/bin/python")
combined.samples.obj.filt <- FindClusters(combined.samples.obj.filt, resolution = 0.8, verbose = FALSE , graph.name = 'wknn' , algorithm = 4 , random.seed = 100)
# png("/Users/sdigiove/Downloads/PCA_COrrectionPhase.png" , height = 5 , width = 9 , res = 300, units = 'in')
DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE)
# dev.off()

# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/UMAP_clusters.png" , height = 5 , width = 9 , res = 300, units = 'in')
DimPlot(combined.samples.obj.filt, label = TRUE, repel = TRUE , group.by = 'seurat_clusters')
# dev.off()


### Plot number of cells per clusters
patientsPerClust <- combined.samples.obj.filt@meta.data |> dplyr::select(dataset ,wknn_res.0.8) |> group_by(wknn_res.0.8,dataset) |> summarise(count=n())

ggplot2::ggplot(data = patientsPerClust , aes(x = wknn_res.0.8 , y = count , fill = dataset) ) +
  geom_bar(position="fill", stat="identity" , color = "black") +
  theme_classic()
ggsave("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/PerClusterCompositionPatients_0_8.png" , height = 5 , width = 9 , units = 'in' , dpi = 300 )


### singleR Annotation
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
ref <- BlueprintEncodeData()

library(SingleR)
perd.combined <- SingleR(test =combined.samples.obj.filt@assays$SCT@data , ref = hpca.se, assay.type.test=1,
                         labels = hpca.se$label.fine , clusters = combined.samples.obj.filt$seurat_clusters)
# perd.combined <- SingleR(test =combined.samples.obj.filt@assays$SCT@counts , ref = ref, assay.type.test=1,
#                          labels = ref$label.fine , clusters = combined.samples.obj.filt$seurat_clusters)
png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/CellANnotationHeatmap.png" , height = 5 , width = 9 , res = 300, units = 'in')
plotScoreHeatmap(perd.combined , clusters = rownames(perd.combined) ,   max.labels = 20)
dev.off()

combined.samples.obj.filt$SingleRClust <- ""
for (clustNum in rownames(perd.combined)){
  combined.samples.obj.filt@meta.data[combined.samples.obj.filt@meta.data$seurat_clusters == clustNum,]$SingleRClust <- perd.combined[clustNum , "pruned.labels"]
}
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/UMAP_patients_annotations.png" , height = 5,width = 9,units = 'in',res = 300)
DimPlot(combined.samples.obj.filt , group.by = "SingleRClust" )
# dev.off()

# Find normal cells according the marker used in the publication https://www.nature.com/articles/s41556-021-00814-7 
DefaultAssay(combined.samples.obj.filt) <- 'SCT'
png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/CD20Plot.png" , height = 7 , width = 4 , units = 'in' , res = 300)
VlnPlot(combined.samples.obj.filt , features = 'MS4A1' , group.by = 'SingleRClust' , pt.size= 0 )+
  theme(legend.position="none")+
  geom_jitter( height = 0.1 , size = 0.5)
dev.off()


### Evaluate cell cycle status
DefaultAssay(combined.samples.obj.filt) <- 'RNA'
combined.samples.obj.filt <- NormalizeData(combined.samples.obj.filt)
combined.samples.obj.filt <- FindVariableFeatures(combined.samples.obj.filt, selection.method = "vst")
combined.samples.obj.filt <- ScaleData(combined.samples.obj.filt, features = rownames(combined.samples.obj.filt))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
combined.samples.obj.filt <- CellCycleScoring(combined.samples.obj.filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/UMAP_Phase.png" , height = 5,width = 9,units = 'in',res = 300)
DimPlot(combined.samples.obj.filt , group.by = 'Phase' )
# dev.off()

PhasePerClust <- combined.samples.obj.filt@meta.data |> dplyr::select(seurat_clusters ,Phase) |> group_by(seurat_clusters,Phase) |> summarise(count=n())

ggplot2::ggplot(data = PhasePerClust , aes(x = seurat_clusters , y = count , fill = Phase) ) +
  geom_bar(position="fill", stat="identity" , color = "black") +
  theme_classic()
# ggsave("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/PerClusterCompositionPhase_0_8.png" , height = 5 , width = 9 , units = 'in' , dpi = 300 )

### Remove clusters with high number of cells ins G1/M or S cell cycle
combined.samples.obj.filt.nocycle <- subset(combined.samples.obj.filt,
                                            subset = seurat_clusters != 4 &
                                              seurat_clusters != 12)


#### Evaluation on clusters from paper heatmap
combined.sites.bt <- rownames(combined.samples.obj.filt.nocycle@assays$ATAC@data) |> as.data.frame() |> 
  separate("rownames(combined.samples.obj.filt@assays$ATAC@data)",c('chrom','start','end') , sep = '-')
cluster1 <- read.delim("/Users/sdigiove/Downloads/cluster1_11k_hg38.bed" , col.names = c('chr','start','end') )
combined.sites.bt.cluster1 <- bedtoolsr::bt.intersect(combined.sites.bt, cluster1 , wa = T)
cluster2 <- read.delim("/Users/sdigiove/Downloads/cluster2_11k_hg38.bed" , col.names = c('chr','start','end') )
combined.sites.bt.cluster2 <- bedtoolsr::bt.intersect(combined.sites.bt, cluster2 , wa = T)
cluster3 <- read.delim("/Users/sdigiove/Downloads/cluster3_11k_hg38.bed" , col.names = c('chr','start','end') ) 
combined.sites.bt.cluster3 <- bedtoolsr::bt.intersect(combined.sites.bt, cluster3 , wa = T)
cluster4 <- read.delim("/Users/sdigiove/Downloads/cluster4_11k_hg38.bed" , col.names = c('chr','start','end') )
combined.sites.bt.cluster4 <- bedtoolsr::bt.intersect(combined.sites.bt, cluster4 , wa = T)

#Create list 
clusters.list <- list('cluster1'= unique(unite(combined.sites.bt.cluster1,'sites',c(V1,V2,V3),sep='-')$sites), 'cluster2'= unique(unite(combined.sites.bt.cluster2,'sites',c(V1,V2,V3),sep='-')$sites) , 'cluster3'= unique(unite(combined.sites.bt.cluster3,'sites',c(V1,V2,V3),sep='-')$sites) , 'cluster4'= unique(unite(combined.sites.bt.cluster4,'sites',c(V1,V2,V3),sep='-')$sites) )

library(BSgenome.Hsapiens.UCSC.hg38.masked)
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(combined.samples.obj.filt.nocycle))) %in% main.chroms)
combined.samples.obj.filt.nocycle[["ATAC"]] <- subset(combined.samples.obj.filt.nocycle[["ATAC"]], features = rownames(combined.samples.obj.filt.nocycle[["ATAC"]])[keep.peaks])
DefaultAssay(combined.samples.obj.filt.nocycle) <- 'ATAC'
combined.samples.obj.filt.nocycle <- AddChromatinModule(
  combined.samples.obj.filt.nocycle,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'ATAC',
  features = clusters.list
)

combined.samples.obj.filt.nocycle@meta.data |> as.data.frame() |> dplyr::select(SingleRClust , cluster1 , cluster2,cluster3,cluster4) |> gather(key = clusters_name , value = ChromatinModules , -SingleRClust) |> 
  ggplot(aes(x = SingleRClust , y = ChromatinModules , fill = SingleRClust )) + #p <- 
  # geom_boxplot(scale = "width" , width = 0.5 , trim = FALSE , color = 'black')+
  geom_violin(width = 0.7 , trim = FALSE , color = 'black')+
  theme_classic()+
  geom_jitter(width = 0.1 , size = 0.3 , alpha=0.1)+
  xlab("")+
  ylab("")+
  theme(legend.position="none") +
  facet_wrap(~clusters_name )
ggsave("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/11kSitesEvaluation/VlnPlotChromatinModule.png",height = 5,width = 9,dpi = 300)

### Find Transcription factors for each clusters of sites
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
combined.samples.obj.filt.nocycle <- AddMotifs(combined.samples.obj.filt.nocycle, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

clusters.list <- list('cluster1'= unique(unite(combined.sites.bt.cluster1,'sites',c(V1,V2,V3),sep='-')$sites), 'cluster2'= unique(unite(combined.sites.bt.cluster2,'sites',c(V1,V2,V3),sep='-')$sites) , 'cluster3'= unique(unite(combined.sites.bt.cluster3,'sites',c(V1,V2,V3),sep='-')$sites) , 'cluster4'= unique(unite(combined.sites.bt.cluster4,'sites',c(V1,V2,V3),sep='-')$sites) )

DefaultAssay(combined.samples.obj.filt.nocycle) <- 'ATAC'
i <- 0 
for (comm in names(clusters.list)) {
  enriched.motifs <- FindMotifs(
    object = combined.samples.obj.filt.nocycle,
    features = clusters.list[[comm]],
    assay = 'ATAC'
  )
  enriched.motifs$community <- comm
  
  if (i == 0)
    final.enriched.motifs = enriched.motifs
  else
    final.enriched.motifs <- rbind(final.enriched.motifs , enriched.motifs)
  
  i <- 1
}

enriched.motifs.filtered <- final.enriched.motifs[final.enriched.motifs$p.adjust <= 0.01,]

enriched.motifs.filtered <- enriched.motifs.filtered |> mutate(ScoringValue= fold.enrichment*-log10(p.adjust))
heatmap.df.TFS <- enriched.motifs.filtered |> dplyr::select(motif.name , community , ScoringValue) |> spread(key = community , value = ScoringValue) |> column_to_rownames('motif.name')

heatmap.df.TFS[is.na(heatmap.df.TFS)] <- 0
heatmap.df.TFS <- t(scale(t(heatmap.df.TFS)))
# col_fun = circlize::colorRamp2(c(-1,0,1), c("blue","white","red"))
r.cluster=hclust(stats::dist(heatmap.df.TFS, method = "euclidean"), method="ward.D2")
c.cluster=hclust(stats::dist(t(heatmap.df.TFS), method = "euclidean"), method="ward.D2")
heatmap.df.TFS[heatmap.df.TFS==0] <- NA

HM <- ComplexHeatmap::Heatmap(
  name = "zscore(normalizedIntesity)",
  as.matrix(heatmap.df.TFS),
  # col=col_fun,
  show_column_names = T,
  show_row_names = T,
  cluster_rows = r.cluster,
  cluster_columns =c.cluster,
  row_split = 4,
  na_col = "black",
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  row_title = " ",
  column_title = " ",
  column_names_gp = grid::gpar(fontsize = 10),
  row_names_gp = grid::gpar(fontsize = 5),
  row_dend_reorder = TRUE,
  show_row_dend = TRUE
)

# png("/Users/sdigiove/Documents/Work/Projects/CardiacDevelopment/data/CellOracle/cebpe/ModuleScore_CEBPA.png" ,  height = 5 , width = 9 , units = 'in' , res = 300)
HM
# dev.off()


### Cicero Analysis to validate the paper findings
library(cicero)
library(SeuratWrappers)
DefaultAssay(combined.samples.obj.filt.nocycle) <- 'ATAC'
combined.samples.obj.filt.cds <- as.cell_data_set(x = combined.samples.obj.filt.nocycle)

set.seed(2017)
combined.samples.obj.filt.cds <- detect_genes(combined.samples.obj.filt.cds)
combined.samples.obj.filt.cds <- estimate_size_factors(combined.samples.obj.filt.cds)
combined.samples.obj.filt.cds <- preprocess_cds(combined.samples.obj.filt.cds, method = "LSI")
combined.samples.obj.filt.cds <- reduce_dimension(combined.samples.obj.filt.cds, reduction_method = 'UMAP', 
                                                  preprocess_method = "LSI")
plot_cells(combined.samples.obj.filt.cds)

combined.samples.obj.filt.cicero <- make_cicero_cds(combined.samples.obj.filt.cds, reduced_coordinates = reducedDims(combined.samples.obj.filt.cds)$UMAP , k = 100 )

data("human.hg19.genome")

## Usually use the whole mouse.mm9.genome ##
## Usually run with sample_num = 100 ##
conns <- run_cicero(combined.samples.obj.filt.cicero, human.hg19.genome) 
# saveRDS(conns , "/Volumes/SSD1T/singleCellBCPALL/CiceroConns.rds")
conns <- readRDS("/Volumes/SSD1T/singleCellBCPALL/CiceroConns.rds")
head(conns)

# chr4:182,807,812-182,886,094
# chr4:182809294 -182809457 ## liftover coord
### Find the DCTD enhancer
dctd.enhancer <- data.frame('chr'=c('chr4') , 'start'=c(182809294) , 'end'=c(182809457))
dctd.combined.coord <- bedtoolsr::bt.intersect(combined.sites.bt, dctd.enhancer , wa = T)
# chr4:182,916,259-182,917,628
dctd.promoter <- data.frame('chr'=c('chr4') , 'start'=c(182916259) , 'end'=c(182917628))
dctd.promoter.combined.coord <- bedtoolsr::bt.intersect(combined.sites.bt, dctd.promoter , wa = T)
dctd.promoter.combined.coord

DefaultAssay(combined.samples.obj.filt) <- 'ATAC'
FeaturePlot(combined.samples.obj.filt , features = c("chr4-182808859-182809757") )

ranges.show <- StringToGRanges("chr4-182808859-182809757")
DimPlot(combined.samples.obj.filt.nocycle , reduction = 'umap', group.by = 'wknn_res.0.8' , label = T)
DimPlot(combined.samples.obj.filt.nocycle , reduction = 'umap', group.by = 'dataset' , label = T)
# Patientse2 <- subset(combined.samples.obj.filt,
#                      subset = dataset == "sample_E2")

cov_plot <- CoveragePlot(
  object = combined.samples.obj.filt.nocycle,
  region = "chr4-182808840-182917840",
  # annotation = TRUE,
  peaks = FALSE,
  ranges = ranges.show,
  group.by = 'seurat_clusters' ,
  # split.by = 'dataset',
  assay = 'ATAC' #wknn_res.0.8'
  # idents = c("Pro-B_cell_CD34+", "B_cell")
)
# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/Harmony/cov_dctd_enhancer_clusters_all.png" , height = 5,res=300,width = 9 , units = 'in')
cov_plot
# dev.off()

tile_plot <- TilePlot(
  object = combined.samples.obj.filt,
  region = "chr4-182808840-182917840",
  # group.by = 'SingleRClust',
  tile.cells = 100,
  idents = c("Pro-B_cell_CD34+", "B_cell")
)
tile_plot

# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/CombinedPlotDCTD.png" , height = 5,res=300,width = 9,units = 'in')
CombineTracks(
  plotlist = list(cov_plot, tile_plot),
  heights = c(10, 10),
  widths = c(10, 1)
)
# dev.off()


conns[conns$Peak1 == "chr4-182808859-182809757"  & conns$Peak2 == "chr4-182917013-182917827", ]
plot_connections(conns, "chr4", 182808840, 182917840,
                 # gene_model = gene_anno, 
                 coaccess_cutoff = .25, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

plot_connections(conns, "chr4", 182808840, 182917840,
                 viewpoint = "chr2_9811600_9811800",
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_track = cons_4c,
                 comparison_connection_width = .5,
                 include_axis_track = F,
                 collapseTranscripts = "longest") 

#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.
# To do this we need the first exon of each transcript
gene_anno <- rtracklayer::readGFF("/Users/sdigiove/Downloads/gencode.v38.annotation.gtf")

# rename some columns to match requirements
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("seqid", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"
names(gene_annotation_sub)[1] <- "chromosome"
combined.samples.obj.filt.cds <- annotate_cds_by_site(combined.samples.obj.filt.cds, gene_annotation_sub)


#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(combined.samples.obj.filt.cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- pData(combined.samples.obj.filt.cds)$num_genes_expressed
names(num_genes) <- row.names(pData(combined.samples.obj.filt.cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
cicero_gene_activities['DCTD',]


dctd.ga <- as.data.frame(cicero_gene_activities['DCTD',])
plotting.coord.original.df <- as.data.frame(combined.samples.obj.filt@reductions$umap@cell.embeddings)
dctd.ga <- dctd.ga[rownames(plotting.coord.original.df) , ]
plotting.coord.original.df$dctd_ga <- dctd.ga


ggplot(data = plotting.coord.original.df , aes(x = UMAP_1 , y = UMAP_2 , color = dctd_ga))+
  geom_point(size = 1)


### FindLoops for Bcells vs CD34+
Bcells <- subset(combined.samples.obj.filt.nocycle,
                 subset = SingleRClust == 'B_cell:immature')
ProBCell.cd34 <- subset(combined.samples.obj.filt,
                        subset = SingleRClust == 'Pro-B_cell_CD34+')


Bcells.cds <- as.cell_data_set(x = Bcells)
Bcells.cicero <- make_cicero_cds(Bcells.cds, reduced_coordinates = reducedDims(Bcells.cds)$UMAP , k = 10 )

ProBCell.cd34.cds <- as.cell_data_set(x = ProBCell.cd34)
ProBCell.cd34.cicero <- make_cicero_cds(ProBCell.cd34.cds, reduced_coordinates = reducedDims(ProBCell.cd34.cds)$UMAP , k = 100 )

Bcells.conns <- run_cicero(Bcells.cicero, human.hg19.genome)
ProBCell.cd34.conns <- run_cicero(ProBCell.cd34.cicero, human.hg19.genome) 

# saveRDS(Bcells.conns , "/Volumes/SSD1T/singleCellBCPALL/BCellCons.rds")
# saveRDS(ProBCell.cd34.conns , "/Volumes/SSD1T/singleCellBCPALL/ProBCellcd34.rds")

# Bcells.conns <- readRDS("/Volumes/SSD1T/singleCellBCPALL/BCellCons.rds")
# ProBCell.cd34.conns <- readRDS("/Volumes/SSD1T/singleCellBCPALL/ProBCellcd34.rds")

#Find Loop In BCells
Bcells.conns[Bcells.conns$Peak1 == "chr4-182808859-182809757"  & Bcells.conns$Peak2 == "chr4-182917013-182917827", ]

names(gene_anno)[1] = 'chromosome'

# png("/Users/sdigiove/Documents/Work/Projects/BCP_ALL/pictures/PlotConnection.png", height = 5 , width = 9 , units = 'in' , res = 300)
plot_connections(Bcells.conns, "chr4", 182808840, 182917840,
                 gene_model = gene_anno,
                 coaccess_cutoff = .20,
                 connection_width = .5,
                 comparison_track = ProBCell.cd34.conns,
                 comparison_connection_width = .5,
                 comparison_coaccess_cutoff = 0.20,
                 collapseTranscripts = "longest")
# dev.off()
