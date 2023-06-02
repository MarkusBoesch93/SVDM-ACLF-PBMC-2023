#Reproducing Figure 1 and Supplementary Figure 1
#Created by Markus Boesch
#loaded in all necessary packages and colours from scRNA_prep
#load in full PBMC sample
#ProgRec contains ACLF-R (R) and ACLF-NR (NR)

PBMC$disease<-factor(x=PBMC$disease, levels = c("Healthy", "AD", "ACLF"))
PBMC$ProgRec<-factor(x=PBMC$ProgRec, levels = c("H", "AD", "R", "NR"))
PBMC$Clusters<-factor(x=PBMC$Clusters, levels = c("Mon", "CD4 T cells", "CD8 T cells", "MAIT", "NK", "B cells", "cDC", "pDC"))
Idents(PBMC)<-PBMC$Clusters
DefaultAssay(PBMC)<-"RNA"
Figure1A<-DimPlot(PBMC, reduction = "umap",  group.by = "Clusters", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_PBMC, order = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=Figure1A, filename="1A UMAP PBMC.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS1A<-DimPlot(PBMC, reduction = "umap",  group.by = "disease", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_disease, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS1A, filename="S1A UMAP PBMC disease.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS1B<-DimPlot(PBMC, reduction = "umap",  group.by = "orig.ident", repel = F, label = F, label.size = 6, pt.size = 0.5, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=FigureS1B, filename="S1B UMAP PBMC orig.svg" ,height=5.2, width=8, units="in", dpi=320)


#S1D
Vlnplotstats<-c("percent.mito", "nFeature_RNA", "nCount_RNA", "S.Score","G2M.Score")
for(i in Vlnplotstats){ 
  Figure<-VlnPlot(PBMC, group.by= "Clusters",features = i , assay = "RNA", cols = colours_PBMC, pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))+NoLegend()
  ggsave(plot=Figure, filename=paste(getwd(),"/S1C_",i,"_PBMC stats Violin.svg", sep = "") ,height=4, width=4, units="in", dpi=320)
}

#1C
Idents(PBMC)<-PBMC@meta.data$orig.ident
table<-table(PBMC@meta.data$Clusters, PBMC@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF","ACLF","ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("Mon", "CD4 T cells", "CD8 T cells", "MAIT", "NK", "B cells", "cDC", "pDC"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_PBMC<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), 
                      axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

##statistics
box_clusters_PBMC2<-box_clusters_PBMC+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = F) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=box_clusters_PBMC, filename="1C Percentage_PBMC.svg" ,height=4, width=5, units="in", dpi=320)
ggsave(plot=box_clusters_PBMC2, filename="1C Percentage_PBMC_stats.svg" ,height=4, width=5, units="in", dpi=320)


#miloR
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
sub<-PBMC
Idents(sub)<-sub$Clusters
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = 50, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=50, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(MILO_obj) ##need distribution peak betwwen 50 and 100, otherwise up k and lower prop
MILO_obj <- countCells(MILO_obj, meta.data = data.frame(colData(MILO_obj)), sample="orig.ident")
MILO_obj <- calcNhoodDistance(MILO_obj, d=50, reduced.dim = "PCA") #~5min
milo.design <- as.data.frame(xtabs(~disease + orig.ident, data=data.frame(colData(MILO_obj))))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$orig.ident
milo.design <- milo.design[colnames(nhoodCounts(MILO_obj)),]
milo.res <- testNhoods(MILO_obj, design=~disease, design.df=milo.design)
MILO_obj <- buildNhoodGraph(MILO_obj)
milo.res <- annotateNhoods(MILO_obj, milo.res, coldata_col = "Clusters")##change to metadata
milo.res$Clusters<-factor(x=milo.res$Clusters, levels = c("Mon", "CD4 T cells", "CD8 T cells", "NK", "B cells", "cDC", "pDC"))
milo_plot<-plotDAbeeswarm(milo.res, group.by = "Clusters", alpha = 0.05)+xlab("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize,
                                                                                           legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) 
ggsave(plot=milo_plot, filename="S1E milo_plot_PBMC_0_05.svg" ,height=5, width=6, units="in", dpi=320)

#1B
Idents(PBMC)<-PBMC$Clusters
PBMC.markers <- FindAllMarkers(object = PBMC, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top20<-PBMC.markers %>% group_by(cluster)%>% top_n(20, avg_log2FC)
PBMC <- ScaleData(PBMC, features = as.character(unique(top20$gene)), assay = "RNA")
Figure1E<-DoHeatmap(subset(PBMC, downsample=100), features = top20$gene, group.by = "Clusters", assay = "RNA", disp.min = -2, group.colors =  colours_PBMC, 
                    size=6, label = T)+theme(axis.text.y.left =element_blank())+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")
ggsave(plot=Figure1E, filename="1E PBMC TOP20 Heatmap.png" ,height=12, width=10, units="in", dpi=320)

#S1E
library(RColorBrewer)
DefaultAssay(PBMC)<-"ADT"
cite_genes<-c("0081-CD14","0083-CD16","0072-CD4","0046-CD8","0050-CD19","0160-CD1c","0047-CD56-NCAM","0034-CD3")
for(i in cite_genes){
  feature_cite<-FeaturePlot(PBMC, label = F, features = i, order = T, pt.size = 1.2, max.cutoff =4, min.cutoff = 0 )+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))+theme( 
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(),axis.line = element_blank(), axis.text = element_blank(),
    axis.ticks = element_blank(), axis.title = element_blank())+NoLegend()+ggtitle("")
  ggsave(plot=feature_cite, filename=paste(getwd(),"/",i,"_Cite_PBMC.png", sep = "") ,height=3, width=3, units="in", dpi=320)
  }


