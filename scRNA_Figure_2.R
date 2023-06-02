#Reproducing Figure 2 and Supplementary Figure 2
#Created by Markus Boesch
#loaded in all necessary packages and colours from scRNA_prep
#subset the Monocytes
#ProgRec contains ACLF-R (R) and ACLF-NR (NR)

Monocytes$disease<-factor(x=Monocytes$disease, levels = c("Healthy", "AD", "ACLF"))
Monocytes$ProgRec<-factor(x=Monocytes$ProgRec, levels = c("H", "AD", "R", "NR"))
Monocytes$subcluster<-factor(x=Monocytes$subcluster, levels = c("cMon", "intMon", "ncMon", "cDC"))

DefaultAssay(Monocytes)<-"RNA"
Figure2A<-DimPlot(Monocytes, reduction = "umap",  group.by = "subcluster", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_Mon, order = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=Figure2A, filename="2A UMAP Mons.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS2A<-DimPlot(Monocytes, reduction = "umap",  group.by = "disease", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_disease, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS2A, filename="S2A UMAP Mons disease.svg" ,height=5.2, width=7, units="in", dpi=320)

Idents(Monocytes)<-Monocytes@meta.data$orig.ident
table<-table(Monocytes@meta.data$subcluster, Monocytes@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF", "ACLF", "ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("cMon", "intMon", "ncMon", "cDC"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_Monocytes<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), 
                            axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

##statistics
box_clusters_Monocytes2<-box_clusters_Monocytes+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = F) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=box_clusters_Monocytes, filename="2B Percentage_Monocytes.svg" ,height=4, width=5, units="in", dpi=320)
ggsave(plot=box_clusters_Monocytes2, filename="2B Percentage_Monocytes_stats.svg" ,height=4, width=5, units="in", dpi=320)

#VlnPlots
Vlnplotgenes<-c("CD14", "S100A12", "VCAN", "FCGR3A", "HLA-DPB1", "HLA-DPA1", "HES4", "CDKN1C", "SIGLEC10", "CD1C", "FCER1A", "CLEC10A")
for(i in Vlnplotgenes){ 
  Figure<-VlnPlot(Monocytes, group.by= "subcluster",features = i , assay = "RNA", cols = colours_Mon, pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))+NoLegend()
  ggsave(plot=Figure, filename=paste(getwd(),"/2D_",i,"_Mon Violin.svg", sep = "") ,height=2, width=2, units="in", dpi=320)
}

##MultiBarHeatmap
Idents(Monocytes)<-Monocytes$subcluster
Mons.markers <- FindAllMarkers(object = Monocytes, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top15<-Mons.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)

RandNRcMon<-subset(Monocytes, subset= (ProgRec %in% c("R", "NR") & subcluster =="cMon" ))
Idents(RandNRcMon)<-RandNRcMon$ProgRec
cMonR_NR.markers <- FindAllMarkers(object = RandNRcMon, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
cMonR_NR_genes<-cMonR_NR.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
genes_mono_deg<-unique(c(cMonR_NR_genes$gene,top15$gene))

genes_mono_deg<-genes_mono_deg[! genes_mono_deg %in% c("XIST", "RPS4Y1", "RPS26", "RPS19", "RPS8", "RPS5", "RPS4X")]
genes_mono_deg<-c(genes_mono_deg,"IFIT2", "IFIT3", "SOCS3", "COMMD9", "ISG15", "ISG20", "OAS3", "MT-1A", "SOD2", "F13A1", "GSDMD", "HSPB1", "ENO1", "LDHA", "NLRP12", "MEFV", "CD74", "TBXAS1", "PECAM1")

Monocytes <- ScaleData(Monocytes, features = genes_mono_deg, assay = "RNA")
Figure2C<-DoMultiBarHeatmap(Monocytes, features = genes_mono_deg, group.by = "subcluster", additional.group.by = "ProgRec", additional.group.sort.by = "ProgRec",assay = "RNA", 
                            disp.min = -2, disp.max = 2, cols.use = list(subcluster=colours_Mon,ProgRec=colours_ProgRec), label = F)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize,
                                y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=Figure2C, filename="2C Heatmap Mons.png" ,height=18, width=10, units="in", dpi=320)

#DotPlot R vs NR cMon signature genes
Rgenes<-c("RETN", "S100A8", "LAP3", "LGALS1", "TMSB10", "PLAC8", "S100P", "CES1", "IFITM3", "SMIM25")
NRgenes<-c("NCF1","VIM", "LGALS2", "TKT", "CRIP1", "CD44", "FTH1", "EVI2B", "TAGLN2", "SGK1")
DotplotRvsNR<-c(Rgenes, NRgenes)
RandNRcMon$ProgRec<-factor(x=RandNRcMon$ProgRec, levels = c("R", "NR"))
Figure2E<-DotPlot(RandNRcMon, group.by= "ProgRec",features = DotplotRvsNR , assay = "RNA", cols = c("steelblue2", "red"), col.min = 0, dot.min = 0.1)+theme(axis.text.x = element_text(angle = 60, hjust=1), 
                  axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) +xlab("")+ylab("")+scale_y_discrete(position = "right")+theme(legend.position = "top")+guides(size="none", colour="none" )
ggsave(plot=Figure2E, filename="2E RNR cMon dotplot.svg" ,height=3, width=10.5, units="in", dpi=320)

#R vs NR cMon violinPlots
other_vlnplotgenes<-c("VCAN", "TKT", "FABP5", "IDH1", "TALDO1", "ITGB2", "B2M", "TREM1","ENO1", "LDHA", "PECAM1", "SOD2","LTA4H", "TBXAS1", "OAS3")
other_vlnplotgenes<-c("SOCS3", "COMMD9", "ISG15", "ISG20", "OAS3", "MT2A", "F13A1", "GSDMD","HSPB1", "MEFV", "NLRP12", "CD74","LYZ", "SELL", "PECAM1", "ITGB2")
other_vlnplotgenes<-c("SREBF1", "KLF6", "ITGAX", "CD83")

for(i in other_vlnplotgenes){ 
  Figure<-VlnPlot(RandNRcMon, group.by= "ProgRec",features = i , assay = "RNA", cols = c("sienna3"  ,  "red4"  ), pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))+NoLegend()
  ggsave(plot=Figure, filename=paste(getwd(),"/2F",i,"_cMon RvsNR new_Violin.svg", sep = "") ,height=2, width=2, units="in", dpi=320)
}


#stats for violinplots
modulescore_statistics_subcluster<-function(signature){
  Idents(RandNRcMon)<-RandNRcMon$ProgRec
  DefaultAssay(RandNRcMon)<-"RNA"
  sub<-RandNRcMon
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-as.matrix(sub@assays$RNA[signature])
  rownames(ahsc_only)<-signature
  cluster.group <- as.character(df.group$cluster) 
  cluster.group[which(cluster.group == "R")] <- 0
  cluster.group[which(cluster.group == "NR")] <- 1
  sub$ProgRec<-factor(x=sub$ProgRec, levels=c( "R" ,"NR"))
  wtDes<-as.data.frame(sub@meta.data) 
  wtDesMat <- model.matrix(~0 + ProgRec, wtDes )
  con <- makeContrasts(RvsNR=  ProgRecR -ProgRecNR,levels=wtDesMat)
  fit <- lmFit(ahsc_only, wtDesMat)
  fit2 <- contrasts.fit(fit, con)
  fit3 <- eBayes(fit2)
  sig_table<-topTable(fit3, adjust.method = "BH", n=Inf)
  return(sig_table)
}

datalist = list()
for(i in other_vlnplotgenes){
  datalist[[i]]<- modulescore_statistics_subcluster(signature = i)
}
big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/2F_cMon RvsNR Violin_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)


#Volcano Plot R vs NR cMon
Idents(RandNRcMon)<-RandNRcMon$ProgRec
cMonRvsNR_volc <- FindMarkers(object = RandNRcMon,ident.1 = "R", ident.2 = "NR",  min.pct = 0.10, logfc.threshold = 0,  slot="data", test.use = "LR", latent.vars = c("nCount_RNA","percent.mito"))
cMonRvsNR_volc = cMonRvsNR_volc[which(cMonRvsNR_volc$p_val_adj<0.05),]
cMonRvsNR_volc<-cMonRvsNR_volc %>% dplyr::rename(cMon_R=pct.1, cMon_NR=pct.2)
openxlsx::write.xlsx(cMonRvsNR_volc, file=paste(getwd(),"/S2D_LR cMonRvsNR_volc",".xlsx",sep = ""), col.names=TRUE, row.names=TRUE)
keyvals <- rep("grey", nrow(cMonRvsNR_volc))
# set the base name/label as "NS"
names(keyvals) <- rep("NS", nrow(cMonRvsNR_volc))
# modify keyvals for variables with fold change >= 1
keyvals[which(cMonRvsNR_volc$avg_log2FC >= 0.5)] <- "sienna3"
names(keyvals)[which(cMonRvsNR_volc$avg_log2FC >= 0.5)] <- "Up"
# modify keyvals for variables with fold change <= -1
keyvals[which(cMonRvsNR_volc$avg_log2FC <= -0.5)] <- "red4"
names(keyvals)[which(cMonRvsNR_volc$avg_log2FC <= -0.5)] <- "Down"
FigureS2D<-EnhancedVolcano(cMonRvsNR_volc,
                           lab = rownames(cMonRvsNR_volc),
                           x = 'avg_log2FC',
                           y = 'p_val_adj', 
                           title = "", subtitle=NULL, gridlines.major = F, gridlines.minor = F, labCol =  NA,
                           colCustom = keyvals, FCcutoff = 0.5, pCutoff = NA, xlim =c(-2,2), xlab = "", ylab = "" , caption = "")+NoLegend()+  
  FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))

ggsave(plot=FigureS2D, filename="S2D cMon R vs NR volcano.svg" ,height=5, width=6, units="in", dpi=320)

##MS1 score
MS1_genes<-c("S100A8", "S100A12", "RETN", "CLU", "MCEMP1", "IL1R2", "CYP1B1", "SELL", "ALOX5AP", "SLC39A8", "PLAC8","ACSL1", "CD163", "VCAN", "HP", "CTSD", "LGALS1", "THBS1", "CES1", "S100P")
RandNRcMon<-AddModuleScore(RandNRcMon, features = list(MS1_genes), name = "MS1")
FigureS2E<-VlnPlot(RandNRcMon, "MS11", group.by = "ProgRec", assay = "RNA", split.by = "ProgRec", pt.size = 0, cols=c("sienna3", "red4"))+ggtitle("")+xlab("")+NoLegend()+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))
ggsave(plot=FigureS2E, filename="S2E Mon MS1 score.svg" ,height=4, width=2, units="in", dpi=320)

modulescore_statistics_subcluster<-function(m){
  Idents(RandNRcMon)<-RandNRcMon$ProgRec
  DefaultAssay(RandNRcMon)<-"RNA"
  sub<-RandNRcMon
  sub=AddModuleScore(object = sub, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-t(as.matrix(sub[[paste(names[m],"1", sep = "")]]))
  rownames(ahsc_only)<-names[m]
  cluster.group <- as.character(df.group$cluster) 
  cluster.group[which(cluster.group == "R")] <- 0
  cluster.group[which(cluster.group == "NR")] <- 1
  sub$ProgRec<-factor(x=sub$ProgRec, levels=c( "R" ,"NR"))
  wtDes<-as.data.frame(sub@meta.data) 
  wtDesMat <- model.matrix(~0 + ProgRec, wtDes )
  con <- makeContrasts(RvsNR=  ProgRecR -ProgRecNR,levels=wtDesMat)
  fit <- lmFit(ahsc_only, wtDesMat)
  fit2 <- contrasts.fit(fit, con)
  fit3 <- eBayes(fit2)
  sig_table<-topTable(fit3, adjust.method = "BH", n=Inf)
  return(sig_table)
}
datalist = list()
scores<-"MS11"
for(m in 1:length(scores)){
  datalist[[scores[m]]]<- modulescore_statistics_subcluster(m = m)
}
big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/S2E MS1 score RvsNR stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)

#MHC score
MHC = c('HLA-DMA','HLA-DMB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB5')
Monocytes<-AddModuleScore(Monocytes, features = list(MHC), name = "MHC", assay = "RNA")
FigureS2C<-VlnPlot(subset(Monocytes, subset=subcluster %in% "cMon"), "MHC1", group.by = "ProgRec", assay = "RNA", split.by = "ProgRec", pt.size = 0,cols=colours_ProgRec)+
  ggtitle("")+xlab("")+NoLegend()+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))
ggsave(plot=FigureS2C, filename="S2C Mon MHC score.svg" ,height=4, width=2, units="in", dpi=320)

#R and NR gene score
Monocytes<-AddModuleScore(Monocytes, features = list(Rgenes), name = "Rgenes", assay = "RNA")
Monocytes<-AddModuleScore(Monocytes, features = list(NRgenes), name = "NRgenes", assay = "RNA")
FigureS2Y<-VlnPlot(subset(Monocytes, subset=subcluster %in% "cMon"), "Rgenes1", group.by = "ProgRec", assay = "RNA", split.by = "ProgRec", pt.size = 0,cols=colours_ProgRec)+
  ggtitle("")+xlab("")+NoLegend()+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))
ggsave(plot=FigureS2Y, filename="S2C Mon Rgenes score.svg" ,height=4, width=2, units="in", dpi=320)
FigureS2Z<-VlnPlot(subset(Monocytes, subset=subcluster %in% "cMon"), "NRgenes1", group.by = "ProgRec", assay = "RNA", split.by = "ProgRec", pt.size = 0,cols=colours_ProgRec)+
  ggtitle("")+xlab("")+NoLegend()+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size=fontsize))
ggsave(plot=FigureS2Z, filename="S2C Mon NRgenes score.svg" ,height=4, width=2, units="in", dpi=320)

#Stats for scoring
modulescore_stats<-function(signature){ 
  sub<-subset(Monocytes, subset=subcluster %in% "cMon")
  Idents(sub)<-sub$ProgRec
  DefaultAssay(sub)<-"RNA"
  sub<-sub
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-t(as.matrix(sub[[paste(signature, sep = "")]]))
  rownames(ahsc_only)<-paste(signature, sep = "")
  cluster.group <- as.character(df.group$cluster) 
  cluster.group[which(cluster.group == "H")] <- 0
  cluster.group[which(cluster.group == "AD")] <- 1
  cluster.group[which(cluster.group == "R")] <- 2
  cluster.group[which(cluster.group == "NR")] <- 3
  sub$ProgRec<-factor(x=sub$ProgRec, levels=c( "H" ,"AD","R", "NR"))
  wtDes<-as.data.frame(sub@meta.data) 
  wtDesMat <- model.matrix(~0 + ProgRec, wtDes )
  colnames(wtDesMat)<- gsub(" ", "_", colnames(wtDesMat))
  cont<-combn(colnames(wtDesMat), 2, FUN = function(x){paste0(x[1], "-", x[2])})
  con <- makeContrasts(levels = wtDesMat, contrasts = cont )
  fit <- lmFit(ahsc_only, wtDesMat)
  fit2 <- contrasts.fit(fit, con)
  fit3 <- eBayes(fit2)
  df<-data.frame()
  for(i in 1:length(colnames(con))){
    sig_table<-topTable(fit3, coef = i, adjust.method = "BH", n=Inf)
    df<-qpcR:::rbind.na(df, as.data.frame(c(colnames(con))[i]),sig_table)}
  return(df)
}
library(limma)
datalist = list()
score_signatures<-c("MHC1", "Rgenes1", "NRgenes1")
for(i in score_signatures){
  datalist[[i]]<- modulescore_stats(signature = i)
}
big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/S2C_cMon MHC, Rgene, NRgene score Violin_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)

#miloR
##for Monocytes
Monocytes@meta.data$subcluster <- droplevels(Monocytes@meta.data$subcluster)
sub<-Monocytes
Idents(sub)<-sub$subcluster
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = 25, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=25, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(MILO_obj) ##need distribution peak betwwen 50 and 100, otherwise up k and lower prop
MILO_obj <- countCells(MILO_obj, meta.data = data.frame(colData(MILO_obj)), sample="orig.ident")
MILO_obj <- calcNhoodDistance(MILO_obj, d=25, reduced.dim = "PCA") #~5min
milo.design <- as.data.frame(xtabs(~disease + orig.ident, data=data.frame(colData(MILO_obj))))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$orig.ident
milo.design <- milo.design[colnames(nhoodCounts(MILO_obj)),]
milo.res <- testNhoods(MILO_obj, design=~disease, design.df=milo.design)
MILO_obj <- buildNhoodGraph(MILO_obj)
milo.res <- annotateNhoods(MILO_obj, milo.res, coldata_col = "subcluster")##change to metadata
milo.res$subcluster<-factor(x=milo.res$subcluster, levels = c("cMon", "intMon", "ncMon", "cDC"))
milo_plot<-plotDAbeeswarm(milo.res, group.by = "subcluster", alpha = 0.05)+xlab("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize,legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) 

ggsave(plot=milo_plot, filename="S2B milo_plot_Mon_0_05.svg" ,height=5, width=6, units="in", dpi=320)






