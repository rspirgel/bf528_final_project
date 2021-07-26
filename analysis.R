BiocManager::install("Seurat")
BiocManager::install("tximeta")
BiocManager::install("tximport")
BiocManager::install("fishpond")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("SingleCellExperiment")
BiocManager::install("org.Hs.eg.db")
install.packages("umap")
reticulate::py_install(packages = 'umap-learn')

library(Seurat)
library(tximport)
library(SingleCellExperiment)
library(org.Hs.eg.db)
library(umap)
library(dplyr)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(ggplot2)


dir <- '/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant'
files <- file.path(dir, "alevin", "quants_mat.gz")
txi <- tximport(files, type="alevin")
cts <- txi$counts
rownames(cts) <- substr(rownames(cts),1,15) #remove decimal from Ensembl ids
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(cts),1,15), keytype = "GENEID", columns = c("SYMBOL"))
cts <- cts[rownames(cts) %in% geneIDs$GENEID,]
rownames(cts) <- geneIDs$SYMBOL #replace ensembl ids with gene names


#import to Seurat
pbmc <- CreateSeuratObject(cts) #does this filter the cells and genes according to instructions?

#labeling the genes that are from mitochondria
#mt.genes <- rownames(sce)[as.logical(seqnames(sce) == "chrM")]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features=mt.genes)
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
VlnPlot(pbmc, features=feats, ncol=3)

pbmc_filtered <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 100 & percent.mt < 30)
VlnPlot(pbmc_filtered, features=feats, ncol=3)


#normalize data
pbmc_filtered <- NormalizeData(pbmc_filtered)

pbmc_filtered <- FindVariableFeatures(
  pbmc_filtered,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc_filtered), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#scaling
all.genes <- rownames(pbmc_filtered)
pbmc_filtered <- ScaleData(pbmc_filtered, features = all.genes)

#create PCA
pbmc_filtered <- RunPCA(pbmc_filtered, features = VariableFeatures(object = pbmc_filtered))
VizDimLoadings(pbmc_filtered, dims = 1:5, reduction = "pca")
DimPlot(pbmc_filtered, reduction = "pca")
print(pbmc_filtered[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(pbmc_filtered, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc_filtered, dims = 1:5, cells = 500, balanced = TRUE)


#determine dimensionality by UMAP
pbmc_filtered <- JackStraw(pbmc_filtered, num.replicate = 100)
pbmc_filtered <- ScoreJackStraw(pbmc_filtered, dims = 1:20)
JackStrawPlot(pbmc_filtered, dims = 1:20)
ElbowPlot(pbmc_filtered)

pbmc_filtered <- FindNeighbors(pbmc_filtered, dims = 1:10)
pbmc_filtered <- FindClusters(pbmc_filtered, resolution = 0.5)
head(Idents(pbmc_filtered), 10)

pbmc_filtered <- RunUMAP(pbmc_filtered, dims = 1:10)
DimPlot(pbmc_filtered, reduction = "umap")

pbmc_filtered <- RunTSNE(pbmc_filtered, dims = 1:10)
tsne <- DimPlot(pbmc_filtered, reduction = "tsne") 
tsne
cell_count <- data.frame(table(Idents(pbmc_filtered)))
ggplot(data=cell_count, aes(x=Var1, y=Freq)) +geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  geom_bar(stat="identity",  fill="steelblue") + xlab("Cluster")


#visualize 
deg <- FindAllMarkers(pbmc_filtered, max.p_val_adj = 0.05,logfc.threshold = 0.8) 


gene_list <- c("GCG","INS","SST","PPY","GHRL","KRT19","CPA1","PDGFRB",
               "VWF", "PECAM1", "CD34", "CD163", "CD68", "IGG", "CD3", 
               "CD8", "TPSAB1", "KIT", "CPA3", 'RGS5', 'PDGFRA', 'SOX10', 
               'SDS', 'TRAC', "IRX1", "IRX2", "ESR1", "MAFA", "NKX6-1", "ETV1", 'IAPP')

cluster_names <- c("Gamma", 1, 2, "Delta", "Ductal", "Alpha", 6, 
                   "Acinar", "Beta", 9, "Stellate", "Macrophage", "Vascular")
cluster_nums <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
names(cluster_names) <- levels(pbmc_filtered)
pbmc_filtered <-RenameIdents(pbmc_filtered,cluster_names)

tsne <- DimPlot(pbmc_filtered, reduction = "tsne") 


#violin plot of gene expression
VlnPlot(pbmc_filtered, features = 'INS')

for (x in gene_list){
print(x)
print(VlnPlot(pbmc_filtered, features = x))
}

#create heatmaps for each cluster
deg_by_cluster<-split(deg, deg$cluster)

for (item in deg_by_cluster){
  print(names(item))
  top5 <- deg %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  print(DoHeatmap(pbmc_filtered, features = top5$gene) 
        + NoLegend())
}

#create heatmap for top5 genes and top3 genes
top5 <- deg %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
htmp<-DoHeatmap(pbmc_filtered, features = top5$gene) + NoLegend()
htmp
FeaturePlot(pbmc_filtered, features = new_gene_list)

cluster2.markers <- FindMarkers(pbmc_filtered, ident.1 = 2, min.pct = 0.25)

top3 <- deg %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)



