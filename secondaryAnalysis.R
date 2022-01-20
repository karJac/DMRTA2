library(DoubletFinder)
library(Matrix)
library(Seurat)
library(readxl)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(parallel)
library(R.devices)
library(pbmcapply)

setwd("/home/rstudio/all/DMRTA2/neurodevelopment/filtered/BTonly")

subfolders <- list.dirs()[-1]
names(subfolders) <- subfolders
subfolders

seuratObjects <- pbmclapply(subfolders,function(i){
     setwd(i)
     barcodes <- read_tsv("barcodes.tsv",col_names=FALSE)
     genes <- read_tsv("genes.tsv",col_names=FALSE)
     mtx <- readMM("matrix.mtx")
     colnames(mtx) <- as.vector(barcodes$X1)
     rownames(mtx) <- as.vector(genes$X2)
     return(mtx)
     setwd("..")
},mc.cores=length(subfolders))

# data for 4 patients was in two mtx tabels
BT338 <- cbind(seuratObjects[[3]],seuratObjects[[4]])
BT363 <- cbind(seuratObjects[[6]],seuratObjects[[7]])
BT364 <- cbind(seuratObjects[[8]],seuratObjects[[9]])
BT397 <- cbind(seuratObjects[[13]],seuratObjects[[14]])

seuratObjects <- seuratObjects[-c(3,4,6,7,8,9,13,14)]
seuratObjects <- append(seuratObjects, c(BT338,BT363,BT364,BT397))
namesSerObj <- lapply(names(seuratObjects),function(x){
        return(substring(x,3,7))
})
names(seuratObjects) <- c(namesSerObj[1:10],"BT338","BT363","BT364","BT307")
names(seuratObjects)

patientsList <- mclapply(1:length(seuratObjects),function(i){
     return(CreateSeuratObject(seuratObjects[[i]], project = names(seuratObjects)[[i]]))
},mc.cores=length(seuratObjects))

getwd()


patientsListQC <- pbmclapply(patientsList,function(patient){
     
     patient <- PercentageFeatureSet(patient, "^MT-", col.name = "percent_mito")
     patient <- PercentageFeatureSet(patient, "^RP[SL]", col.name = "percent_ribo")
     patient <- PercentageFeatureSet(patient, "^HB[^(P)]", col.name = "percent_hb")
     patient <- PercentageFeatureSet(patient, "PECAM1|PF4", col.name = "percent_plat")
     
     feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
     
     nCountFeaturePlot <- FeatureScatter(patient, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
     
     #filter the extreme outliers to get properl violin plot
     selected_mito <- WhichCells(patient, expression = percent_mito < 15)
     test <- subset(patient, cells = selected_mito)
     
     QCplot <- VlnPlot(test, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
          NoLegend()
     
    
     
     exp <- patient@assays$RNA@counts
     
     #mouse_mito <- read_excel("Mouse.MitoCarta3.0.xls") # lepsze przyrownanie niż "^mt-" ale nie skonczylem go pisac
     #mouseMito <- mouse_mito[,3]
     
     
     selected_c <- WhichCells(patient, expression = nFeature_RNA > 200)
     #dodaj usuwanie komórek z mniejsza niz 3k? UMI
     selected_f <- rownames(patient)[Matrix::rowSums(patient) > 20]
     
     data.filt <- subset(patient, features = selected_f, cells = selected_c)
     dim(data.filt)
     
     selected_mito <- WhichCells(data.filt, expression = percent_mito < 8) #normaly ~5,5
     selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 5) # cuz it correlates with high levels of mito genes

     data.filt <- subset(data.filt, cells = selected_mito)
     data.filt <- subset(data.filt, cells = selected_ribo)
     
     dim(data.filt)
     
     table(data.filt$orig.ident)
     
     feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
     
     data.filt <- PercentageFeatureSet(data.filt, "^MT-", col.name = "percent_mito") #idk if these two lines are neccesery
     data.filt <- PercentageFeatureSet(data.filt, "^RP[SL]", col.name = "percent_ribo")
     
     QCfilteredPlot <- VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + 
          NoLegend()
     
     # Compute the relative expression of each gene per cell Use sparse matrix
     # operations, if your dataset is large, doing matrix devisions the regular way
     # will take a very long time.
     par(mar = c(4, 8, 2, 1))
     C <- data.filt@assays$RNA@counts
     C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
     most_expressed <- order(apply(C, 1, median), decreasing = T)[30:1]
     tmp <- data.frame(t(C[most_expressed, ]))
     tmp2 <- pivot_longer(tmp, cols = colnames(tmp))
     tmp2$name <- factor(tmp2$name, levels = colnames(tmp), ordered = TRUE)
     highestExprsPlot <- ggplot(tmp2, aes(x=value,y=name, fill = name)) + geom_boxplot()
     
     patID <- as.character(unique(patient@meta.data$orig.ident))
     
     R.devices::suppressGraphics({
          png(paste("QCplot_",patID,".png",sep=''),height = 1200, width = 1200)
          print(plot_grid(QCplot,nCountFeaturePlot,QCfilteredPlot,highestExprsPlot, nrow=2,
                          labels = c(paste("QC for patient: ",patID) ,"nCountFeature", "QC after filtering", "most expressed"),
                          vjust = 0.85, greedy = FALSE)) # improve annotations of graphs because they are overlapping with figures
          dev.off()
     })
     return(data.filt)
}, mc.cores = length(patientsList))



patientsListFilt <- pbmclapply(patientsListQC,function(data.filt){
     
     # Filter MALAT1
     data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
     # Filter Mitocondrial
     #data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
     # Filter Ribosomal
     data.filt <- data.filt[!grepl("^RP[SL][[:digit:]]", rownames(data.filt)), ]
     # Filter Ribosomal rRNA
     data.filt <- data.filt[!grepl("rRNA",  ignore.case = TRUE, rownames(data.filt)), ]
     
     return(data.filt)
}, mc.cores = length(patientsListQC))

m.s.genes <- cc.genes.updated.2019$s.genes
m.g2m.genes <- cc.genes.updated.2019$g2m.genes

patientsListNorm <- mclapply(patientsListFilt,function(data.filt){
     
     # Before running CellCycleScoring the data need to be normalized and
     # logtransformed.
     data.filt = NormalizeData(data.filt)
     
     glio <- CellCycleScoring(object = glio, g2m.features = m.g2m.genes, 
                                    s.features = m.s.genes)
     cellCycleGraph <- VlnPlot(glio, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
                                 pt.size = 0.1)
     
     patID <- as.character(unique(data.filt@meta.data$orig.ident))
     
     # R.devices::suppressGraphics({
     #         png(paste("CellCycleplot_",patID,".png",sep=''))
     #         print(cellCycleGraph)
     #         dev.off()
     # })
     # 
     data.filt = FindVariableFeatures(data.filt, verbose = F)
     data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), 
                           verbose = F)
     data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
     data.filt = RunUMAP(data.filt, dims = 1:20, verbose = F)
     
     return(data.filt)
     
}, mc.cores=length(patientsListFilt))




patientsAfterDF2 <- pbmclapply(patientsListNorm,function(patient){
     
     
     sweep.res <- paramSweep_v3(patient) 
     sweep.stats <- summarizeSweep(sweep.res,GT = FALSE) 
     bcmvn <- find.pK(sweep.stats)
     
     pK=as.numeric(as.character(bcmvn$pK))
     BCmetric=bcmvn$BCmetric
     pK_choose = pK[which(BCmetric %in% max(BCmetric))]
     
     par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
     #plot(x = pK, y = BCmetric, pch = 16,type="b",  #visualtion of BCmetrics (but function finds the maximum automaticly)
     #     col = "blue",lty=1)
     #abline(v=pK_choose,lwd=2,col='red',lty=2)
     #title("The BCmvn distributions")
     #text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
     
     # define the expected number of doublet cellscells.
     nExp <- round(ncol(patient) * 0.039)  # we expect 3.9% doublets (10X doublets prop)
     patient <- doubletFinder_v3(patient, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10)
     n=dim(patient@meta.data)[[2]]-1
     m=dim(patient@meta.data)[[2]]
     colnames(patient@meta.data)[[n]] <- "pANN"
     colnames(patient@meta.data)[[m]] <- "DoubletFinder"
     return(patient)
     
}, mc.cores=length(patientsListNorm)) 


# Merge datasets into one single seurat object

alldata <- merge(patientsAfterDF2[[1]], c(patientsAfterDF2[2:length(patientsAfterDF2)]))
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata, verbose = F)
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"), assay = "RNA")
alldata <- RunPCA(alldata, verbose = F, npcs = 50)
#alldata <- JackStraw(alldata,reduction = "pca", dims = 30, prop.freq = 0.1) #find the amount of PC u want to use in UMAP
alldata <- RunUMAP(alldata, dims = 1:30, verbose = F)

DF.name = colnames(alldata@meta.data)[grepl("DoubletFinder", colnames(alldata@meta.data))]

png("doublets.png",width = 1200, height = 500)
cowplot::plot_grid(ncol = 2, DimPlot(alldata, group.by = "orig.ident") + NoAxes(), 
                   DimPlot(alldata, group.by = "DoubletFinder") + NoAxes())
dev.off()

png("doubletsViolin.png",width = 1200, height = 500)
VlnPlot(alldata, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

alldata = alldata[, alldata@meta.data[, DF.name] == "Singlet"]




########################### Dimensionality reduction #############################


top20 <- head(VariableFeatures(alldata), 20)

png("volcanoPlot.png")
LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
dev.off()

plot_grid(ncol = 3, DimPlot(alldata, reduction = "pca", group.by = "orig.ident", 
                            dims = 1:2), DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4), 
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 5:6))
VizDimLoadings(alldata, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)

png("PCAcontributions.png")
ElbowPlot(alldata, reduction = "pca", ndims = 30) # amonut of variance explained by each PC
dev.off()

#### Check if your data needs integration
png("origIdent.png",height = 800, width = 1800)
DimPlot(alldata, reduction = "umap", group.by = "orig.ident", dims=c(8,2)) + 
     ggplot2::ggtitle(label = "UMAP_on_PCA") +
     facet_grid(~orig.ident)
dev.off()


saveRDS(alldata,"tmp.rds")





#################### if the batch effect exist (cells are clustering by sample) then we are running sample integration #####
#                                optional:
#########################  SAMPLE INTEGRATION #########################
# 
# library(venn)
# 
# 
# alldata.list <- SplitObject(alldata, split.by = "orig.ident")
# 
# for (i in 1:length(alldata.list)) {
#         alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
#         alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst",
#                                                   nfeatures = 2000, verbose = FALSE)
# }
# 
# which(alldata.list[[1]]@assays$RNA@var.features == "ENSG00000142700.11")
# 
# hvgs_per_dataset <- lapply(alldata.list, function(x) {
#         x@assays$RNA@var.features
# })
# venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1,
#            cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)
# 
# alldata.features <- SelectIntegrationFeatures(object.list = alldata.list)
# alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,
#                                           reduction = "cca", anchor.features = alldata.features)
# 
# allFeatures <- lapply(alldata.list, row.names) %>% Reduce(intersect, .)
# 
# alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA", features.to.integrate = allFeatures,
#                              verbose = T)
# 
# alldata.int@active.assay
# 
# # Run Dimensionality reduction on integrated space
# alldata.int <- ScaleData(alldata.int, verbose = FALSE)
# alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)
# alldata.int <- RunUMAP(alldata.int, dims = 1:30)
# 
# #check the effect of integration on UMAP and/or clustering
# DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident") +
#         ggplot2::ggtitle(label = "UMAP_on_PCA") +
#         facet_grid(~orig.ident)
# 
# saveRDS(alldata.int,"afterIntegration.rds")
# alldata <- alldata.int

############################ Clustering ############################

library(pheatmap)
library(enrichR)
library(rafalib)
library(clustree)
library(leiden)


alldata@active.assay


alldata <- FindNeighbors(alldata, dims = 1:20, k.param = 60, prune.SNN = 1/15)
names(alldata@graphs)

pheatmap(alldata@graphs$RNA_nn[1:200, 1:200], col = c("white", "black"), border_color = "grey90", 
         legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)

# Clustering with Leiden (algorithm 4)
for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
     alldata <- FindClusters(alldata, graph.name = names(alldata@graphs)[[2]], resolution = res, algorithm = 4) #4 = Leiden ALGORITHM
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.

plot_grid(ncol = 3, DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.0.5") + 
               ggtitle("leibens_0.5"), DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.1") + 
               ggtitle("leibens_1"), DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.2") + 
               ggtitle("leibens_2"))

png("clustree.png")
clustree(alldata@meta.data, prefix = "RNA_snn_res.")
dev.off()

# K-means clustering
# for (k in c(5, 7, 10, 12, 15, 17, 20)) {
#         alldata@meta.data[, paste0("kmeans_", k)] <- kmeans(x = alldata@reductions[["pca"]]@cell.embeddings, 
#                                                             centers = k, nstart = 100)$cluster
# }
# 
# plot_grid(ncol = 3, DimPlot(alldata, reduction = "umap", group.by = "kmeans_5") + 
#                   ggtitle("kmeans_5"), DimPlot(alldata, reduction = "umap", group.by = "kmeans_10") + 
#                   ggtitle("kmeans_10"), DimPlot(alldata, reduction = "umap", group.by = "kmeans_15") + 
#                   ggtitle("kmeans_15"))
# 
# clustree(alldata@meta.data, prefix = "kmeans_")

saveRDS(alldata,"afterClustering.rds")




############################ Differential gene expression ############################

# Set the identity as leibens with resolution

library(tidyr)

sel.clust = "RNA_snn_res.0.25" 

alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)

# plot this clustering
plot_grid(ncol = 3, DimPlot(alldata, label = T) + NoAxes(), DimPlot(alldata, group.by = "orig.ident", shuffle = TRUE) + NoAxes())


FeaturePlot(alldata, reduction = "umap", features = c("SOX2","CD3D","CD14","PTPRC","DMRTA2","RGS5","PDGFRB"), order = T, slot = "data", combine = T)


macExp <- alldata@assays$RNA@data

metadata <- cbind(alldata@meta.data$RNA_snn_res.1,alldata@reductions$umap@cell.embeddings,alldata@meta.data$orig.ident)
colnames(metadata) <- c("clusters","UMAP1","UMAP2","patID")
metadata <- data.frame(metadata)
levels(metadata$clusters) <- sort(unique(metadata$clusters))
metadata$clusters <- as.factor(metadata$clusters)
levels(metadata$clusters) <- as.character(sort(as.numeric(levels(metadata$clusters))))
metadata$UMAP1 <- as.numeric(metadata$UMAP1)
metadata$UMAP2 <- as.numeric(metadata$UMAP2)

saveRDS(metadata,"metadata_notIntegrated.rds")
saveRDS(macExp,"macExpn_notIntegrated.rds")
