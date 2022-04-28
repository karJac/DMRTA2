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



################## INTERPRETATION ##################

library(ggpubr)
library(viridis)

theme_feature<-theme_classic()+
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(hjust=0), legend.title = element_blank(),
          legend.position = "right", legend.text=element_text(size=14),
          title = element_text(size=40), axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16))
YlOrRd<-colorRampPalette(c("Yellow","Orange", "Red"))

macExp <- data.frame(macExp)
macExp <- t(macExp)

macExp <- cbind(metadata,macExp)

macExp <- macExp[,which(unique(colnames(macExp)) %in% colnames(macExp))]

#UMAP WITH CLUSTERS_ID
png("meh.png",width = 1400, height = 1400)
ggplot(metadata, aes(x=UMAP1, y=UMAP2, color=clusters), )+
    geom_jitter(size=0.7, alpha=0.5)+
    geom_text(aes(label=clusters), size=4, alpha=0.5, check_overlap=T, color="black")+
    #facet_grid(.~day)+
    guides(colour = guide_legend(override.aes = list(size=7))) +
    scale_color_discrete(breaks = levels(metadata$clusters))
dev.off()



#GENE EXPRESSION PLOTS

genes = c("SOX2","CD3D","CD14","PTPRC","DMRTA2","RGS5","PDGFRB")
graphs<-list()
graphs <- mclapply(1:length(genes), function(i) {
    if (i %in% colnames(macExp)){
        plot <- ggplot(macExp, aes(x = UMAP1, y = UMAP2)) +
            geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name(genes[[i]]))) + #czemu okregi a nie koła?
            labs(title = genes[[i]]) +
            scale_color_gradientn(colors = (c(
                viridis(50),
                YlOrRd(50)
            ))) + 
            # breaks=c(min,max))+
            # , labels=c("min", "max"))+
            xlab("UMAP_1") +
            ylab("UMAP_2") +
            theme_feature +
            theme(legend.text = element_text(size = 20), legend.position = "right") 
        return(plot)}
    else {
        return(NULL)
    }
}, mc.cores = length(genes))
#png(file = "exp.png", width = 4000, height = 2000)
ggarrange(plotlist = graphs, ncol = 3, nrow = 2, common.legend = F)
#dev.off()



colnames(macExp)[which(grepl("tdTomato",colnames(macExp)))]


FeaturePlot(alldata, reduction = "umap", features = c("SOX2","CD3D","CD14","PTPRC","DMRTA2","RGS5","SOX10"), order = T, slot = "data", combine = T)

FeaturePlot(alldata, reduction = "umap", features = c("FGF12","SOX10","MPZL1","SMOC1","TCF7L2","OLIG1","CSPG4"), order = T, slot = "data", combine = T)




############ DIFFERENTIONAL ANALYSIS #############

alldata <- readRDS("afterClustering.rds")
suppressPackageStartupMessages({
    library(Seurat)
    library(venn)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(enrichR)
    library(rafalib)
    library(clusterProfiler)
})



plot_grid(ncol = 2, DimPlot(alldata, label = T) + NoAxes(), DimPlot(alldata, group.by = "orig.ident", shuffle = TRUE) + NoAxes())

enrichR::listEnrichrDbs()

a = alldata@meta.data$RNA_snn_res.0.25
cols <- colnames(alldata)[which(a %in% c(3,12,8,4,6,11,5,10))]
sox2 <- subset(alldata, cells = cols)
WhichCells(sox2, expression = DMRTA2 > 0)

sox2$expression <- sox2@assays$RNA@data["DMRTA2",] > 0
Idents(sox2) <- sox2$expression 

marks <- FindAllMarkers(sox2,only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
library(clusterProfiler)
#1. podaj wszystkie geny z marks (jeśli mają sensowne log2change)
#2. zaktualizuj R lub
#3. uzyj enrichR

enrich_results <- enrichr(genes = rownames(marks), databases = "GO_Biological_Process_2017b")

neg <- subset(alldata,subset = DMRTA2 == 0)

enrich_results <- enrichr(genes = rownames(dmrt), databases = "GO_Biological_Process_2017b")

par(mfrow = c(1, 1), mar = c(3, 25, 2, 1))
barplot(height = -log10(enrich_results$GO_Biological_Process_2017b$P.value)[30:1], names.arg = enrich_results$GO_Biological_Process_2017b$Term[30:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6)
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)



enrich_results2 <- enrichr(genes = rownames(dmrt), databases = "KEGG_2019_Human")

par(mfrow = c(1, 1), mar = c(3, 25, 2, 1))
barplot(height = -log10(enrich_results2$KEGG_2019_Human$P.value)[30:1], names.arg = enrich_results2$KEGG_2019_Human$Term[30:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6)
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)



enrich_results3 <- enrichr(genes = rownames(dmrt), databases = "WikiPathways_2019_Human")

par(mfrow = c(1, 1), mar = c(3, 25, 2, 1))
barplot(height = -log10(enrich_results3$WikiPathways_2019_Human$P.value)[30:1], names.arg = enrich_results3$WikiPathways_2019_Human$Term[30:1], 
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6)
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)

biomart
#biological proceses, molecular functions
library("org.Hs.eg.db")
library()

alldata <- readRDS("alldata.rds")
marks <- FindAllMarkers(tumor,only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
genes <- rownames(marks)
#geneE <- bitr(genes, "GENENAME", "ENTREZID", org.Hs.eg.db, drop = TRUE)

mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'useast.ensembl.org')
mart <- useDataset('hsapiens_gene_ensembl', mart)
geneE <- getBM(values = genes,
                       filters = "hgnc_symbol",
                       mart = mart,
                       attributes = c("entrezgene_id"))


geneE[] <- lapply(geneE, as.character)
geneE <- as.list(geneE)
geneE <- unlist(geneE)
geneE <- as.character(geneE)

ggoMF <- groupGO(gene   = geneE,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

universe <- rownames(alldata)
universe <- getBM(values = universe,
                  filters = "hgnc_symbol",
                  mart = mart,
                  attributes = c("entrezgene_id"))
universe[] <- lapply(universe, as.character)
universe <- as.list(universe)
universe <- unlist(universe)
universe <- as.character(universe)

egoBP <- enrichGO(gene        = geneE,
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

png("egoBP.png", width=800,height=900)
dotplot(egoBP, showCategory=30) + ggtitle("GO Biological Processes") + theme(title = element_text(size=25))
dev.off()


kk  <- gseKEGG(geneList     = geneE,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

png("gseKK.png", width=800,height=900)
dotplot(kk, showCategory=30)  + ggtitle("KEGG") + theme(title = element_text(size=25))
dev.off()

map04370 <- pathview(gene.data  = geneE,
                     pathway.id = "map04370",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneE)), cpd=1))

library(ReactomePA)
react <- enrichPathway(gene=geneE, pvalueCutoff = 0.05, readable=TRUE)

png("react.png", width=1000,height=900)
dotplot(react, showCategory=30) + ggtitle("ReactomePA") + theme(title = element_text(size=25))
dev.off()

############## other stuff for Marta ###############
 
a = alldata@meta.data$RNA_snn_res.1
cols <- colnames(alldata)[-(which(a %in% 25))]
tumor <- subset(alldata, cells = cols)

tumor <- readRDS("tumor.rds")

marks <- readRDS("marks.rds")
geneE <- readRDS("geneE.rds")
universe <- readRDS("universe.rds")

macExp <- alldata@assays$RNA@data

png("lol.png")
DimPlot(alldata , group.by = "RNA_snn_res.1", label = T) + NoAxes()
dev.off()


dmrt <- subset(alldata, cells = WhichCells(tumor, expression = DMRTA2 > 0))
dmrt2 <- subset(dmrt, subset = orig.ident == "BT322")
macExp <- dmrt@assays$RNA@data

genes=c("DMRTA2","PDGFRB","ACTA2","CSPG4","DES","HIGD1B","S1PR3","MCAM","CD248")
graphs <- lapply(genes, umap)
png("peri_dmrt.png",width=1200,height = 500)
ggarrange(plotlist = graphs, ncol = 5, nrow = 2, common.legend = F)
dev.off()

png("test.png",width=1500,height=250)
ggplot(macExp, aes(DMRTA2,ACTA2))+
     geom_jitter(size = 0.5, alpha = 0.5)  +
     facet_grid(~patID)
dev.off()

png("test.png")
VlnPlot(dmrt2,features = "ACTA2")
dev.off()

png("test.png",width=1500,height=250)
ggplot(macExp, aes(UMAP1,UMAP2))+
     geom_jitter(size = 0.5, alpha = 0.5)  +
     facet_grid(~patID)
dev.off()

png("test.png")
DimPlot(alldata, group.by = "RNA_snn_res.1", label = T)
dev.off()


clu19 <- subset(alldata, subset = RNA_snn_res.1 == 19)
png("test.png")
DimPlot(clu19)
dev.off()


umap <- function(i){
     plot <- ggplot(macExp, aes(x = UMAP1, y = UMAP2)) +
          geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name(i))) + #czemu okregi a nie koła?
          labs(title = i) +
          scale_colour_gradientn(colours = rev(inferno(50))) +
          # breaks=c(min,max))+
          # , labels=c("min", "max"))+
          xlab("UMAP_1") +
          ylab("UMAP_2") +
          theme_feature +
          theme(legend.text = element_text(size = 20), legend.position = "right")
     return(plot)
}

genes=c("DMRTA2","OLIG2","CSPG4","NES","SOX2")
test(genes)
graphs <- lapply(genes, umap)
png("peri_dmrt.png",width=3400,height = 1500)
ggarrange(plotlist = graphs, ncol = 5, nrow = 3, common.legend = F)
dev.off()

set.seed(42)
rows <- sample(nrow(macExp))
macExp <- macExp[rows, ]

dmrt <- subset(alldata,subset = DMRTA2 > 0)

genes = c("DMRTA2","GFAP","PDGFRA","OLIG1","OLIG2","RBFOX3","TUBB3","NES","PECAM1","ENG","PDGFRB","ACTA2","CSPG4","PTPRC","SOX2") 

pos=list()
neg=list()
pieList <- lapply(genes,function(i){
     pos <- length(which(macExp[,i] > 0))/dim(macExp)[[1]]
     neg <- 1 - pos
     return(print(pie(x=c(pos,neg),labels = c("POSITIVE","NEGATIVE"), col = c("dodgerblue","firebrick2"), main = i))) #kolory jakos zgrane z umapem
     })


graphs <- lapply(genes,function(i){
dp <- which(macExp[,"DMRTA2"] > 0 & macExp[,i] > 0)
pp <- which(macExp[,"DMRTA2"] > 0 & macExp[,i] == 0)
ip <- which(macExp[,"DMRTA2"] == 0 & macExp[,i] > 0)
dn <- which(macExp[,"DMRTA2"] == 0 & macExp[,i] == 0)

macExp$colors <- c(rep("black",dim(macExp)[[1]]))
macExp$colors[dp] <- "dp"
macExp$colors[pp] <- "pp"
macExp$colors[ip] <- "ip"
macExp$colors[dn] <- "dn"

return(ggplot() +
     geom_jitter(data = macExp, aes(x = UMAP1, y = UMAP2, color = colors),size=1) +
     scale_color_manual(values = c("dp" = "limegreen", "pp" = "firebrick1", "ip" = "khaki1", "dn" = "ivory4"), name = "Combinations", labels = c(paste("DMRTA2+ ",i,"+",sep=""), paste("DMRTA2+ ",i,"-",sep=""), paste("DMRTA2- ",i,"+",sep=""),paste("DMRTA2- ",i,"-",sep=""))) +
     ggtitle(i) +
     theme_feature +
     theme(legend.text = element_text(size = 20), legend.position = "right")
)
})

umap <- function(i){
     plot <- ggplot(macExp, aes(x = UMAP1, y = UMAP2)) +
          geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name(i))) + #czemu okregi a nie koła?
          labs(title = i) +
          scale_colour_gradientn(colours = rev(inferno(50))) +
          # breaks=c(min,max))+
          # , labels=c("min", "max"))+
          xlab("UMAP_1") +
          ylab("UMAP_2") +
          theme_feature +
          theme(legend.text = element_text(size = 20), legend.position = "right")
     return(plot)
}
test <- function(j){
     for (i in j){
          if (identical(colnames(macExp)[which(grepl(i,colnames(macExp)))], character(0))){
               print(paste(i," ","is not in genes pool"))
          }
          else{
               print(colnames(macExp)[which(grepl(i,colnames(macExp)))])
          }
     }
}



glio <- subset(alldata, idents = c(3,4,5,6,8,10,11,12))



