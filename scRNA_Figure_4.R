#Reproducing Figure 4 and Supplementary Figure 4-6
#Created by Markus Boesch
#loaded in all necessary packages and colours from scRNA_prep
#subset the Tcells, NKs and B cells
#ProgRec contains ACLF-R (R) and ACLF-NR (NR)

##Tcells
#4A
Tcells$disease<-factor(x=Tcells$disease, levels = c("Healthy", "AD", "ACLF"))
Tcells$ProgRec<-factor(x=Tcells$ProgRec, levels = c("H", "AD", "R", "NR"))
Tcells$subcluster<-factor(x=Tcells$subcluster, levels = c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL", "CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"))

DefaultAssay(Tcells)<-"RNA"
Figure4A<-DimPlot(Tcells, reduction = "umap",  group.by = "subcluster", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_Tcells, order = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=Figure4A, filename="4A UMAP Tcells.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS4A<-DimPlot(Tcells, reduction = "umap",  group.by = "disease", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_disease, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS4A, filename="S4A UMAP Tcells disease.svg" ,height=5.2, width=7, units="in", dpi=320)

Idents(Tcells)<-Tcells@meta.data$orig.ident
table<-table(Tcells@meta.data$subcluster, Tcells@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF","ACLF","ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL", "CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_Tcells<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,60))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

#Progrec
Idents(Tcells)<-Tcells@meta.data$orig.ident
table<-table(Tcells@meta.data$subcluster, Tcells@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF", "ACLF", "ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL", "CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_Tcells_prog<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,70))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")
ggsave(plot=box_clusters_Tcells_prog, filename="4D Percentage_Tcells_progRec.png" ,height=4, width=5, units="in", dpi=320)


##statistics
box_clusters_Tcells2<-box_clusters_Tcells+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = F) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=box_clusters_Tcells, filename="4D Percentage_Tcells.svg" ,height=4, width=5, units="in", dpi=320)
ggsave(plot=box_clusters_Tcells2, filename="4D Percentage_Tcells_stats.svg" ,height=4, width=5, units="in", dpi=320)

#4B CD4 Heatmap
CD4_cells<-subset(Tcells, subset=subcluster %in% c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL"))
Idents(CD4_cells)<-CD4_cells$subcluster
CD4_cells.markers <- FindAllMarkers(object = CD4_cells, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top15<-CD4_cells.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)

CD4_cells$subcluster<-droplevels(CD4_cells$subcluster)
CD4_cells <- ScaleData(CD4_cells, features = c("CD4",top15$gene), assay = "RNA")
Figure4B<-DoMultiBarHeatmap(CD4_cells, features = c("CD4",top15$gene), group.by = "subcluster", additional.group.by = "ProgRec", additional.group.sort.by = "ProgRec",assay = "RNA", disp.min = -2, disp.max = 2, cols.use = list(subcluster=colours_Tcells[1:4],ProgRec=colours_ProgRec), label = F)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=Figure4B, filename="4B Heatmap CD4.png" ,height=18, width=10, units="in", dpi=320)

#4C CD8 Heatmap
CD8_cells<-subset(Tcells, subset=subcluster %in% c("CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"))
Idents(CD8_cells)<-CD8_cells$subcluster
CD8_cells.markers <- FindAllMarkers(object = CD8_cells, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top15<-CD8_cells.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
CD8_cells$subcluster<-droplevels(CD8_cells$subcluster)

CD8_cells <- ScaleData(CD8_cells, features = c("CD8B", "CD8A",top15$gene), assay = "RNA")
Figure4C<-DoMultiBarHeatmap(CD8_cells, features = c("CD8B", "CD8A",top15$gene), group.by = "subcluster", additional.group.by = "ProgRec", additional.group.sort.by = "ProgRec",assay = "RNA", disp.min = -2, disp.max = 2, cols.use = list(subcluster=colours_Tcells[5:9],ProgRec=colours_ProgRec), label = F)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=Figure4C, filename="4C Heatmap CD8.png" ,height=18, width=10, units="in", dpi=320)

#Tcells GO scoring
#4E
Activation_T=list(c("ABL1", "ADA", "ADAM17", "ADAM8", "ADORA2A", "AGER", "AIF1", "AIRE", "AKT1", "AMBRA1", "ANXA1", "AP3B1", "AP3D1", "APBB1IP", "ARG1", "ARG2", "ARMC5", "ATF2", "ATG5", "ATP7A", "AZI2", "B2M", "BAD", "BATF", "BAX", "BCL10", "BCL11B", "BCL2", "BCL3", "BCL6", "BMP4", "BTN2A2", "BTN3A1", "CADM1", "CAMK4", "CARD11", "CASP3", "CASP8", "CAV1", "CBFB", "CBLB", "CCDC88B", "CCL19", "CCL2", "CCL21", "CCL5", "CCND3", "CCR2", "CCR6", "CCR7", "CCR9", "CD151", "CD160", "CD1C", "CD1D", "CD2", "CD209", "CD24", "CD27", "CD274", "CD276", "CD28", "CD300A", "CD3D", "CD3E", "CD3G", "CD4", "CD40LG", "CD44", "CD46", "CD47", "CD48", "CD5", "CD55", "CD6", "CD7", "CD70", "CD74", "CD80", 
                    "CD81", "CD83", "CD84", "CD86", "CD8A", "CD8B", "CDH26", "CDK6", "CEACAM1", "CEBPB", "CGAS", "CHD7", "CLC", "CLEC4A", "CLEC4D", "CLEC4E", "CLEC4G", "CLEC7A", "CLECL1", "CLPTM1", "CORO1A", "CR1", "CRACR2A", "CRTAM", "CSK", "CTLA4", "CTNNB1", "CTPS1", "CTSG", "CTSL", "CXADR", "CYP26B1", "CYRIB", "DDOST", "DHPS", "DLG1", "DLG5", "DLL4", "DNAJA3", "DOCK2", "DOCK8", "DPP4", "DROSHA", "DUSP10", "DUSP22", "DUSP3", "EBI3", "EFNB1", 
                    "EFNB2", "EFNB3", "EGR1", "EGR3", "EIF2AK4", "ELF4", "ENTPD7", "EOMES", "EPO", "ERBB2", "F2RL1", "FADD", "FANCA", "FANCD2", "FCER1G", "FCGR2B", "FCHO1", "FGL1", "FGL2", "FKBP1A", "FKBP1B", "FLOT2", "FOXJ1", "FOXN1", "FOXO3", "FOXP1", "FOXP3", "FUT7", "FYN", "FZD5", "FZD7", "FZD8", "GATA3", "GBA", "GJA1", "GLI2", "GLI3", "GLMN", "GPAM", "GPNMB", "GPR18", "GPR183", "GPR89A", "GPR89B", "GSN", "HAVCR2", "HES1", "HFE", "HHLA2", "HLA-A", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-E", "HLA-G", "HLX", "HMGB1", "HSH2D", "HSPD1", "HSPH1", "ICAM1", "ICOS", "ICOSLG", "IDO1", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IFNE", "IFNG", "IFNK", "IFNL1", "IFNW1", "IGF1", "IGF2", "IGFBP2", "IHH", "IL10", "IL12A", "IL12B",
                    "IL12RB1", "IL15", "IL18", "IL18R1", "IL1A", "IL1B", "IL1RL2", "IL2", "IL20RB", "IL21", "IL23A", "IL23R", "IL27", "IL27RA", "IL2RA", "IL36B", "IL4", "IL4I1", "IL4R", "IL6", "IL6R", "IL6ST", "IL7", "IL7R", "ILDR2", "INS", "IRF1", "IRF4", "ITCH", "ITGAL", "ITK", "ITPKB", "JAG2", "JAK2", "JAK3", "JAML", "JMJD6", "KAT2A", "KDELR1", "KIF13B", "KIT", "KLRC1", "KLRC4", "KLRK1", "LAG3", "LAPTM5", "LAT", "LAX1", "LCK", "LCP1", "LEF1", "LEP", "LEPR", "LFNG", "LGALS1", "LGALS3", "LGALS7B", "LGALS9", "LGALS9B", "LGALS9C", "LIG4", "LILRB1", "LILRB2", "LILRB4", "LMBR1L", "LMO1", "LOXL3", "LRRC32", "LY9", "LYN", "MAD1L1", "MAFB", "MALT1", "MAP3K8", "MAPK8IP1", "MARCHF7", "MDK", "METTL3", "MICA", "MICB", "MIR181C", "MIR21", "MIR27A", "MIR30B", "MR1", "MSN", "MTOR", "MYB", "MYH9", "NCAPH2", "NCK1", "NCK2", "NCKAP1L", "NCSTN", "NDFIP1", "NFKBID", "NFKBIZ", "NHEJ1", "NKAP", "NKG7", "NKX2-3", "NLRC3", "NLRP3", "NOD2", "NRARP", "P2RX7", "PAG1", "PATZ1", "PAWR", "PCK1", "PDCD1LG2", "PDE5A", "PELI1", "PIK3CA", "PIK3CD", "PIK3CG", "PIK3R6", "PKNOX1", "PLA2G2D", "PLA2G2F", "PNP", "PPP3CA", "PPP3CB", "PRDM1", "PRDX2", "PRELID1", "PREX1", "PRKAR1A", "PRKCQ", "PRKCZ", "PRKDC", "PRNP", "PRR7", "PSEN1", "PSG9", "PSMB10", "PSMB11",
                    "PTGER4", "PTPN11", "PTPN2", "PTPN22", "PTPN6", "PTPRC", "PYCARD", "RAB27A", "RAB29", "RABL3", "RAC2", "RAG1", "RAG2", "RARA", "RASAL3", "RASGRP1", "RC3H1", "RC3H2", "RELB", "RHOA", "RHOH", "RIPK2", "RIPK3", "RIPOR2", "RORA", "RORC", "RPL22", "RPS3", "RPS6", "RSAD2", "RUNX1", "RUNX2", "RUNX3", "SART1", "SASH3", "SCGB1A1", "SCRIB", "SDC4", "SELENOK", "SEMA4A", "SFTPD", "SH3RF1", "SHH", "SIRPA", "SIRPB1", "SIRPG", "SIT1", "SLA2", "SLAMF1", "SLAMF6", "SLAMF7", "SLAMF9", "SLC11A1", "SLC46A2", "SLC7A1", "SMAD3", "SMAD7", "SOCS1", "SOCS5", "SOCS6", "SOD1", "SOS1", "SOS2", "SOX12", "SOX13", "SOX4", "SP3", "SPI1", "SPINK5", "SPN", "SPTA1", "SRC", "SRF", "STAT3", "STAT5B", "STAT6", "STK11", "STOML2", "SYK", "TARM1", "TBX21", "TCF7", "TCIRG1", "TESPA1", "TFRC", "TGFBR2", "THEMIS", "THY1", "TIGIT", "TMEM131L", "TMEM98", "TMIGD2", "TNFAIP8L2", "TNFRSF13C", "TNFRSF14", "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFSF11", "TNFSF13B", "TNFSF14", "TNFSF18", "TNFSF4", "TNFSF8", "TNFSF9", "TOX", "TP53", "TRAF3IP2", "TRAF6", "TREML2", "TREX1", "TSC1", "TWSG1", "TYK2", "VAV1", "VCAM1", "VNN1", "VSIG4", "VSIR", "VTCN1", "WAS", "WDFY4", "WNT1", "WNT4", "XBP1", "XCL1", "YES1", "ZAP70", "ZBTB1", "ZBTB16", "ZBTB7B", "ZC3H12A", "ZC3H8", "ZEB1", "ZFP36L1", "ZFP36L2", "ZFPM1", "ZMIZ1", "ZNF683", "ZP3","ZP4"))

Apoptosis= list(c('TNFSF10', 'TNFRSF10A','TNFRSF10B','FASLG','FAS','FADD','TNF','TNFRSF1A','TRADD','CFLAR','CASP8','CASP10','CASP6','CASP3','CASP7','BID','BAX','BAK1','DIABLO','SEPTIN4','HTRA2','CYCS','APAF1','CASP9','PRF1','GZMB','TUBA1B','TUBA4A','TUBA3C','TUBA1A','TUBA1C','TUBA8','TUBA3E','TUBA3D','TUBAL3','MCL1','ACTG1','ACTB','SPTA1','SPTAN1','LMNA','LMNB1','LMNB2','PARP1','PARP2','PARP3','PARP4','DFFA','DFFB','ENDOG','AIFM1','ERN1','TRAF2','ITPR1','ITPR2','ITPR3','CAPN1','CAPN2','CASP12','EIF2AK3','EIF2S1','ATF4','DDIT3','CTSB','CTSC','CTSD','CTSF','CTSH','CTSK','CTSL','CTSO','CTSS','CTSV','CTSW','CTSZ','BIRC2','BIRC3','XIAP','BIRC5','BCL2L11','BCL2L1','BCL2','DAXX','RIPK1','DAB2IP','MAP3K5','MAPK8','MAPK10','MAPK9','BAD','JUN','FOS','TP53','HRK','MAP3K14','CHUK','IKBKB','IKBKG','NFKBIA','NFKB1','RELA','PTPN13','GADD45A','GADD45B','GADD45G','TRAF1','BCL2A1','ATM','PIDD1','TP53AIP1','BBC3','PMAIP1','CASP2','NGF','NTRK1','IL3','IL3RA','CSF2RB','PIK3CA','PIK3CD','PIK3CB','PIK3R1','PIK3R2','PIK3R3','PDPK1','AKT1','AKT2','AKT3','HRAS','KRAS','NRAS','RAF1','MAP2K1','MAP2K2','MAPK1','MAPK3'))

Cytokines=list(c("CCL2", "CCL3", "CCL4", "CCL5", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL8", "CXCL9", "FGF2", "IFNG", "IL10", "IL12A", "IL15RA", "IL18", "IL1B", "IL1RN", "IL2", "IL4", "IL6", "IL7", "LTA", "PDGFB", "TNF", "VEGFA"))

S100= list(c('S100A1','S100A2','S100A3','S100A4','S100A5','S100A6','S100A7','S100A7A','S100A7L2','S100A7P1','S100A7P2','S100A8','S100A9','S100A10','S100A11','S100A12','S100A13','S100A14','S100A15A','S100A16','S100B','S100G','S100P','S100Z'))

cytotoxicity=list(c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7"))

types=c(Activation_T, Apoptosis, Cytokines, S100, cytotoxicity)
names=c("Activation_T", "Apoptosis", "Cytokines", "S100", "cytotoxicity")

cells<-c("CD4_TEM", "CD4_CTL", "CD8_GZMK", "CD8_CTL")
for(i in 1:length(cells)){
  subc<-subset(Tcells, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){ tryCatch({
    Idents(subc)<-subc$ProgRec
    subc=AddModuleScore(object = subc, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
    Figure<-VlnPlot(subc, group.by= "ProgRec",features = paste(names[m],"1", sep = "") , assay = "RNA", cols = colours_ProgRec[3:4], pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_blank())+NoLegend()
    
    ggsave(plot=Figure, filename=paste(getwd(),"/4E_",cellname,names[m], "_score_Violin.svg", sep = "") ,height=3, width=3, units="in", dpi=320)
  },error=function(e){})}}

modulescore_statistics_subcluster<-function(m){
  Idents(sub)<-sub$ProgRec
  sub=AddModuleScore(object = sub, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-t(as.matrix(sub[[paste(names[m],"1", sep = "")]]))
  rownames(ahsc_only)<-paste0(names[m],"_", cellname, sep="")
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
library(limma)
datalist = list()
for(i in 1:length(cells)){
  sub<-subset(Tcells, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){
    datalist[[paste0(names[m],"_", cellname, sep="")]]<- modulescore_statistics_subcluster(m = m)
  }}

big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/4E_Tcells RvsNR Violin_scores_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)

#4F metabolic comparisons
library(stats)
a=read.table('metabolism.term.txt',sep='\t')
b=read.table('KEGG.txt',sep='\t')
colnames(a)=c('ID','Term','Category')
colnames(b)=c('ID','Term','GeneID','Gene')
term=unique(a$ID)
uniqueb<-unique(b$ID)
uniqueterm<-a[a$ID %in% uniqueb,]$ID
cells_suba<-c("CD4_TN", "CD4_TCM","CD4_TEM", "CD4_CTL", "CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT")

suba<-subset(Tcells, subset=subcluster %in% cells_suba)
suba$subcluster<-droplevels(suba$subcluster)
Idents(suba)<-suba$subcluster
DefaultAssay(suba)='RNA'
for(i in 1:length(uniqueterm)){
  su=subset(b,ID==as.character(uniqueterm[i]))
  suba=AddModuleScore(object = suba,features = list(su$Gene),name=uniqueterm[i],assay = 'RNA',ctrl = 10)
}

score=suba@meta.data
types=paste(uniqueterm,'1',sep='')
cells=unique(score$subcluster)
ps=matrix(nrow=0,ncol=5)
for(i in c(1:length(cells))){tryCatch({
  tmp=subset(score,subcluster==cells[i])
  h=subset(tmp,ProgRec=="NR")
  m=subset(tmp,ProgRec=='R')
  for(j in 1:length(types)){
    pv1=t.test(m[,types[j]],h[,types[j]])
    pvd1=data.frame(cells[i],types[j],pv1$p.value,(pv1$estimate[1]-pv1$estimate[2]))
    pvd1$Compare='R'
    colnames(pvd1)=c('cells','types','pvalue','diff','Compare')
    pvd=rbind(pvd1)
    ps=rbind(ps,pvd)
  }
},error=function(e){})}

ps$group[ps$diff>0]='up'
ps$group[ps$diff<0]='down'
ps$group[ps$diff==0]='none'
ps=subset(ps,pvalue<0.05)
ps$groups[ps$diff>0]=1
ps$groups[ps$diff<0]=-1
ps$groups[ps$diff==0]=0
ps$logp=log(-log10(ps$pvalue))
ps$logps=(ps$groups)*(ps$logp)
ps$ID=substr(ps$types,1,8)
psall=merge(ps,a,by='ID')

svg('4F score.point_R_NR_Tcells_flipped.svg',width=9,height=10)
p2=ggplot(psall,aes(x=Term,y=cells,fill=logps,shape=Compare))+geom_point(stroke=0.3,size=3.5)+scale_fill_gradient2(low='navy',mid='white',high='#8B0000',midpoint=0)+coord_flip()
p2=p2+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(color='black'))+scale_shape_manual(values=c(21,22,24))+facet_grid(Category~group, scales = "free", space = "free")
p2=p2+theme(axis.text.x = element_text(angle=90,vjust=1,hjust=1, size = fontsize), axis.text = element_text(size=fontsize), legend.text = element_text(size=fontsize), strip.text = element_text(size=fontsize))+labs(x='',y='')+scale_color_manual(values=c('black','black'))
p2$labels$fill<-"-log10(P)"
print(p2)
dev.off()

#S4F
Cytokines=list(c("CCL2", "CCL3", "CCL4", "CCL5", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL8", "CXCL9", "FGF2", "IFNG", "IL10", "IL12A", "IL15RA", "IL18", "IL1B", "IL1RN", "IL2", "IL4", "IL6", "IL7", "LTA", "PDGFB", "TNF", "VEGFA"))

DefaultAssay(PBMC)<-"RNA"
subc<-PBMC
subc$subcluster<-factor(x=subc$subcluster, levels = c("cMon", "intMon", "ncMon", "cDC","CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL", "CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT", "CD56low NK", "CD56high NK", "B cells", "Plasma cells", "pDC"))
subc=AddModuleScore(object = subc, features = Cytokines,name = "Cytokines", assay = "RNA",ctrl = 10)
Figure<-VlnPlot(subc, group.by= "subcluster",features = "Cytokines1" , assay = "RNA", cols = colours_allsub, pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_blank())+NoLegend()
ggsave(plot=Figure, filename=paste(getwd(),"/S4F_cytokine_score_allPBMC_Violin.svg", sep = "") ,height=4, width=6, units="in", dpi=320)
#S4D miloR
#for CD4 Tcells
CD4_cells<-subset(Tcells, subset=subcluster %in% c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL"))
CD4_cells@meta.data$subcluster <- droplevels(CD4_cells@meta.data$subcluster)
sub<-CD4_cells
Idents(sub)<-sub$subcluster
DefaultAssay(sub)<-"RNA"
sub<-FindVariableFeatures(object =sub, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
DefaultAssay(sub)<-"integrated"
sub <- ScaleData(object = sub,
                 vars.to.regress = c("nCount_RNA", "percent.mito"))
sub <- RunPCA(object = sub, npcs = 30)
sub <- RunUMAP(sub, dims = 1:30, reduction = "pca", n.components = 5)
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = 30, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
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
milo.res$subcluster<-factor(x=milo.res$subcluster, levels = c("CD4_TN", "CD4_TCM","CD4_TEM","CD4_CTL"))
milo_plot<-plotDAbeeswarm(milo.res, group.by = "subcluster", alpha = 0.05)+xlab("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize,legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) 

ggsave(plot=milo_plot, filename="S4B milo_plot_CD4_0_05.svg" ,height=5, width=6, units="in", dpi=320)

#for CD8 Tcells
CD8_cells<-subset(Tcells, subset=subcluster %in% c("CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"))
CD8_cells@meta.data$subcluster <- droplevels(CD8_cells@meta.data$subcluster)
sub<-CD8_cells
Idents(sub)<-sub$subcluster
DefaultAssay(sub)<-"RNA"
sub<-FindVariableFeatures(object =sub, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
DefaultAssay(sub)<-"integrated"
sub <- ScaleData(object = sub,
                 vars.to.regress = c("nCount_RNA", "percent.mito"))
sub <- RunPCA(object = sub, npcs = 30)
sub <- RunUMAP(sub, dims = 1:30, reduction = "pca", n.components = 5)
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = 30, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = "PCA")
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
milo.res$subcluster<-factor(x=milo.res$subcluster, levels = c("CD8_TN","CD8_GZMK", "CD8_CTL", "gd_T", "MAIT"))
milo_plot<-plotDAbeeswarm(milo.res, group.by = "subcluster", alpha = 0.05)+xlab("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize,legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) 

ggsave(plot=milo_plot, filename="S4B milo_plot_CD8_0_05.svg" ,height=5, width=6, units="in", dpi=320)

##B cells
#S5A
Bcells$disease<-factor(x=Bcells$disease, levels = c("Healthy", "AD", "ACLF"))
Bcells$ProgRec<-factor(x=Bcells$ProgRec, levels = c("H", "AD", "R", "NR"))
Bcells$subcluster<-factor(x=Bcells$subcluster, levels = c("B cells", "Plasma cells"))
Bcells$subcluster<-droplevels(Bcells$subcluster)

DefaultAssay(Bcells)<-"RNA"
FigureS5A<-DimPlot(Bcells, reduction = "umap",  group.by = "subcluster", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_Bcell, order = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS5A, filename="S5A UMAP Bcells.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS5B<-DimPlot(Bcells, reduction = "umap",  group.by = "disease", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_disease, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS5B, filename="S5B UMAP Bcells disease.svg" ,height=5.2, width=7, units="in", dpi=320)

Idents(Bcells)<-Bcells@meta.data$orig.ident
table<-table(Bcells@meta.data$subcluster, Bcells@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF","ACLF","ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("B cells", "Plasma cells"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_Bcells<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

##statistics
box_clusters_Bcells2<-box_clusters_Bcells+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = F) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=box_clusters_Bcells, filename="S5D Percentage_Bcells.svg" ,height=4, width=5, units="in", dpi=320)
ggsave(plot=box_clusters_Bcells2, filename="S5D Percentage_Bcells_stats.svg" ,height=4, width=5, units="in", dpi=320)

#S5C Heatmap
Idents(Bcells)<-Bcells$subcluster
Bcells.markers <- FindAllMarkers(object = Bcells, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top15<-Bcells.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
Bcells <- ScaleData(Bcells, features = top15$gene, assay = "RNA")
FigureS5C<-DoMultiBarHeatmap(Bcells, features = top15$gene, group.by = "subcluster", additional.group.by = "ProgRec", additional.group.sort.by = "ProgRec",assay = "RNA", disp.min = -2, disp.max = 2, cols.use = list(subcluster=colours_Bcell,ProgRec=colours_ProgRec), label = F)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=FigureS5C, filename="S5C Heatmap Bcells.png" ,height=10, width=10, units="in", dpi=320)

#S5E GO scoring B cells
activation_B=list(c("ABL1", "ADA", "ADAM17", "ADGRG3", "AHR", "AICDA", "AKAP17A", "APLF", "ATAD5", "ATM", "BAD", "BAK1", "BANK1", "BATF", "BAX", "BCL2", "BCL3", "BCL6", "BLK", "BLNK", "BST1", "BST2", "BTK", "C17orf99", "CARD11", "CASP3", "CASP8", "CCR6", "CD180", "CD19", "CD22", "CD27", "CD28", "CD300A", "CD320", "CD38", "CD40", "CD40LG", "CD70", "CD74", "CD79A", "CD79B", "CD81", "CD86", "CDH17", "CDKN1A", "CEBPG", "CHRNA4", "CHRNB2", "CLCF1", "CMTM7", "CR1", "CR2", "CTLA4", "CTPS1", "CXCR5", "DCAF1", "DCLRE1C", "DLL1", "DNAJB9", "DOCK10", "DOCK11", "EP300", "EPHB2", "ERCC1", "EXO1", "EXOSC3", "EXOSC6", "EZH2", "FCGR2B", "FCRL1", "FCRL3", "FLT3", "FNIP1", "FOXJ1", "FOXP3", "FZD9", "GAPT", "GON4L", "GPR183", "GPS2", "HDAC4", "HDAC5", "HDAC9", "HHEX", "HMCES", "HMGB3", "HSPD1", "ICOSLG", "ID2", "IFNA1", "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNB1", "IFNE", "IFNK", "IFNW1", "IGBP1", "IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGHV1-18", "IGHV1-24", "IGHV1-3", "IGHV1-45", "IGHV1-58", "IGHV1-69", "IGHV1-69-2", "IGHV1-69D", "IGHV1OR15-1", "IGHV2-26", "IGHV2-5", "IGHV2-70", "IGHV2-70D", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-16", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-35", "IGHV3-38", "IGHV3-43", "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64", "IGHV3-64D", "IGHV3-66", "IGHV3-7", "IGHV3-72", "IGHV3-73", "IGHV3-74", "IGHV4-28", "IGHV4-31", "IGHV4-34", "IGHV4-39", "IGHV4-4", "IGHV4-59", "IGHV4-61", "IGHV5-10-1", "IGHV5-51", "IGHV6-1", "IGHV7-4-1", "IGHV7-81", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGLL1", "IGLL5", "IKZF3", "IL10", "IL11", "IL13", "IL2", "IL21", "IL27RA", "IL4", "IL4I1", "IL5", "IL6", "IL7", "IL7R", "IL9", "INHA", "INHBA", "INPP5D", "IRF2BP2", "IRF8", "IRS2", "ITFG2", "ITGA4", "ITGB1", "ITM2A", "JAK3", "KIT", "KLF6", "KMT5B", "KMT5C", "LAPTM5", "LAT2", "LAX1", "LEF1", "LFNG", "LGALS1", "LIG4", "LRRC8A", "LYL1", "LYN", "MAD2L2", "MALT1", "MEF2C", "MFNG", "MIF", "MIR17HG", "MIR185", "MIR19A", "MLH1", "MMP14", "MNDA", "MS4A1", "MSH2", "MSH6", "MZB1", "NBN", "NCKAP1L", "NDFIP1", "NFAM1", "NFATC2", "NHEJ1", "NKX2-3", "NOD2", "NOTCH2", "NSD2", "NTRK1", "ONECUT1", "PARP3", "PAWR", "PAXIP1", "PCID2", "PELI1", "PHB", "PHB2", "PIK3CD", "PKN1", "PLCG2", "PLCL2", "POU2AF1", "PPP2R3C", "PRKCB", "PRKCD", "PRKDC", "PTK2B", "PTPN2", "PTPN6", "PTPRC", "PTPRJ", "RABL3", "RAG1", "RAG2", "RASGRP1", "RBPJ", "RC3H1", "RIF1", "RNF168", "RNF8", "SAMSN1", "SASH3", "SFRP1", "SH3KBP1", "SHLD1", "SHLD2", "SHLD3", "SKAP2", "SLAMF8", "SLC15A4", "SLC25A5", "SLC39A10", "SP3", "SPI1", "ST3GAL1", "STAT5B", "STAT6", "SUPT6H", "SWAP70", "SYK", "SYVN1", "TBC1D10C", "TBX21", "TCF3", "TCIRG1", "TFRC", "TGFB1", "THEMIS2", "THOC1", "TICAM1", "TIRAP", "TLR4", "TLR9", "TNFAIP3", "TNFRSF13B", "TNFRSF13C", "TNFRSF21", "TNFRSF4", "TNFSF13", "TNFSF13B", "TNFSF4", "TNIP2", "TP53", "TP53BP1", "TPD52", "TRAF3IP2", "TRBC1", "TRBC2", "TRDC", "TXLNA", "TYROBP", "UNG", "VAV3", "VCAM1", "WNT3A", "XBP1", "YY1", "ZAP70", "ZBTB7A", "ZFP36L1", "ZFP36L2"))#GO:0042113

types=c(activation_B)
names=c("activation_B")

cells<-c("B cells", "Plasma cells")
for(i in 1:length(cells)){
  subc<-subset(Bcells, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){ tryCatch({
    Idents(subc)<-subc$ProgRec
    subc=AddModuleScore(object = subc, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
    Figure<-VlnPlot(subc, group.by= "ProgRec",features = paste(names[m],"1", sep = "") , assay = "RNA", cols = colours_ProgRec[3:4], pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_blank())+NoLegend()
    
    ggsave(plot=Figure, filename=paste(getwd(),"/S5E_",cellname,names[m], "_score_Violin.svg", sep = "") ,height=3, width=3, units="in", dpi=320)
  },error=function(e){})}}

modulescore_statistics_subcluster<-function(m){
  Idents(sub)<-sub$ProgRec
  sub=AddModuleScore(object = sub, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-t(as.matrix(sub[[paste(names[m],"1", sep = "")]]))
  rownames(ahsc_only)<-paste0(names[m],"_", cellname, sep="")
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
library(limma)
datalist = list()
for(i in 1:length(cells)){
  sub<-subset(Bcells, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){
    datalist[[paste0(names[m],"_", cellname, sep="")]]<- modulescore_statistics_subcluster(m = m)
  }}

big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/S5E_Bcells RvsNR Violin_scores_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)


##NKs
#S6A
NKs$disease<-factor(x=NKs$disease, levels = c("Healthy", "AD", "ACLF"))
NKs$ProgRec<-factor(x=NKs$ProgRec, levels = c("H", "AD", "R", "NR"))
NKs$subcluster<-factor(x=NKs$subcluster, levels = c("CD56low NK", "CD56high NK"))
NKs$subcluster<-droplevels(NKs$subcluster)

DefaultAssay(NKs)<-"RNA"
FigureS6A<-DimPlot(NKs, reduction = "umap",  group.by = "subcluster", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_NK, order = T, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS6A, filename="S6A UMAP NKs.svg" ,height=5.2, width=7, units="in", dpi=320)

FigureS6B<-DimPlot(NKs, reduction = "umap",  group.by = "disease", repel = F, label = F, label.size = 6, pt.size = 0.5, cols=colours_disease, order = F, shuffle = T)+ggtitle("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+NoLegend()
ggsave(plot=FigureS6B, filename="S6B UMAP NKs disease.svg" ,height=5.2, width=7, units="in", dpi=320)

#S6D
Idents(NKs)<-NKs@meta.data$orig.ident
table<-table(NKs@meta.data$subcluster, NKs@active.ident)
table<-as.data.frame.matrix(table)
table<-table[order(row.names(table)),]
table <- table[ order(as.numeric(row.names(table))), ]
table.perc<-apply(table[],2,function (x){(x/sum(x))*100})
table.perc.round<-as.data.frame(t(round(table.perc, digits = 2))) ##rounded number
table.perc.round$patient<-rownames(table.perc.round)
table.perc.round$disease<-c("ACLF","ACLF","ACLF","ACLF","ACLF","AD","AD","AD", "Healthy","Healthy","Healthy","ACLF","ACLF","ACLF","ACLF","ACLF")
table.perc.round$disease<-factor(table.perc.round$disease, levels = c("Healthy", "AD", "ACLF"))

longer_table<-table.perc.round %>% tidyr::pivot_longer(cols =c("CD56low NK", "CD56high NK"),  names_to = "Cluster", values_to = "Percentage")

box_clusters_NKs<-ggboxplot(longer_table, x="Cluster", y="Percentage", add = "jitter", color="disease", palette=colours_disease, ylim=c(0,100))+theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank())+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))+xlab("")+ylab("")

##statistics
box_clusters_NKs2<-box_clusters_NKs+geom_pwc(
  aes(group = disease), tip.length = 0,
  method = "t_test", label = "p.adj.format",
  bracket.nudge.y = -0.08,  hide.ns = F) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(plot=box_clusters_NKs, filename="S6D Percentage_NKs.svg" ,height=4, width=5, units="in", dpi=320)
ggsave(plot=box_clusters_NKs2, filename="S6D Percentage_NKs_stats.svg" ,height=4, width=5, units="in", dpi=320)

#S6C Heatmap

Idents(NKs)<-NKs$subcluster
NKs.markers <- FindAllMarkers(object = NKs, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
top15<-NKs.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
NKs <- ScaleData(NKs, features =c("NCAM1", "FCGR3A", top15$gene), assay = "RNA")
FigureS6C<-DoMultiBarHeatmap(NKs, features = c("NCAM1", "FCGR3A", top15$gene), group.by = "subcluster", additional.group.by = "ProgRec", additional.group.sort.by = "ProgRec",assay = "RNA", disp.min = -2, disp.max = 2, cols.use = list(subcluster=colours_NK,ProgRec=colours_ProgRec), label = F)+ scale_fill_gradientn(colors = c("steelblue2", "white", "red2"))+guides(colour="none")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize, legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize ))
ggsave(plot=FigureS6C, filename="S6C Heatmap NKs.png" ,height=10, width=10, units="in", dpi=320)

#S6F GO scoring
Activation_NK=list(c("CD300A", "ELF4", "HSPH1", "IL12A", "IL12B", "IL15", "IL18", "IL23A", "IL23R", "JAK2", "RASAL3", "TYK2", "ZBTB7B"))#GO:0051132
Apoptosis= list(c('TNFSF10', 'TNFRSF10A','TNFRSF10B','FASLG','FAS','FADD','TNF','TNFRSF1A','TRADD','CFLAR','CASP8','CASP10','CASP6','CASP3','CASP7','BID','BAX','BAK1','DIABLO','SEPTIN4','HTRA2','CYCS','APAF1','CASP9','PRF1','GZMB','TUBA1B','TUBA4A','TUBA3C','TUBA1A','TUBA1C','TUBA8','TUBA3E','TUBA3D','TUBAL3','MCL1','ACTG1','ACTB','SPTA1','SPTAN1','LMNA','LMNB1','LMNB2','PARP1','PARP2','PARP3','PARP4','DFFA','DFFB','ENDOG','AIFM1','ERN1','TRAF2','ITPR1','ITPR2','ITPR3','CAPN1','CAPN2','CASP12','EIF2AK3','EIF2S1','ATF4','DDIT3','CTSB','CTSC','CTSD','CTSF','CTSH','CTSK','CTSL','CTSO','CTSS','CTSV','CTSW','CTSZ','BIRC2','BIRC3','XIAP','BIRC5','BCL2L11','BCL2L1','BCL2','DAXX','RIPK1','DAB2IP','MAP3K5','MAPK8','MAPK10','MAPK9','BAD','JUN','FOS','TP53','HRK','MAP3K14','CHUK','IKBKB','IKBKG','NFKBIA','NFKB1','RELA','PTPN13','GADD45A','GADD45B','GADD45G','TRAF1','BCL2A1','ATM','PIDD1','TP53AIP1','BBC3','PMAIP1','CASP2','NGF','NTRK1','IL3','IL3RA','CSF2RB','PIK3CA','PIK3CD','PIK3CB','PIK3R1','PIK3R2','PIK3R3','PDPK1','AKT1','AKT2','AKT3','HRAS','KRAS','NRAS','RAF1','MAP2K1','MAP2K2','MAPK1','MAPK3'))

Cytokines=list(c("CCL2", "CCL3", "CCL4", "CCL5", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL8", "CXCL9", "FGF2", "IFNG", "IL10", "IL12A", "IL15RA", "IL18", "IL1B", "IL1RN", "IL2", "IL4", "IL6", "IL7", "LTA", "PDGFB", "TNF", "VEGFA"))

S100= list(c('S100A1','S100A2','S100A3','S100A4','S100A5','S100A6','S100A7','S100A7A','S100A7L2','S100A7P1','S100A7P2','S100A8','S100A9','S100A10','S100A11','S100A12','S100A13','S100A14','S100A15A','S100A16','S100B','S100G','S100P','S100Z'))

cytotoxicity=list(c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7"))

types=c(Activation_NK,Apoptosis,Cytokines,S100,cytotoxicity)
names=c("Activation_NK","Apoptosis","Cytokines","S100","cytotoxicity")

cells<-c("CD56low NK", "CD56high NK")
for(i in 1:length(cells)){
  subc<-subset(NKs, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){ tryCatch({
    Idents(subc)<-subc$ProgRec
    subc=AddModuleScore(object = subc, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
    Figure<-VlnPlot(subc, group.by= "ProgRec",features = paste(names[m],"1", sep = "") , assay = "RNA", cols = colours_ProgRec[3:4], pt.size=0)+theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), plot.title = element_blank())+NoLegend()
    
    ggsave(plot=Figure, filename=paste(getwd(),"/S6F_",cellname,names[m], "_score_Violin.svg", sep = "") ,height=3, width=3, units="in", dpi=320)
  },error=function(e){})}}

modulescore_statistics_subcluster<-function(m){
  Idents(sub)<-sub$ProgRec
  sub=AddModuleScore(object = sub, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
  df.group <- data.frame(umi = names(Idents(sub)), 
                         cluster = as.character(sub@meta.data$ProgRec), 
                         stringsAsFactors = F)
  ahsc_only<-t(as.matrix(sub[[paste(names[m],"1", sep = "")]]))
  rownames(ahsc_only)<-paste0(names[m],"_", cellname, sep="")
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
library(limma)
datalist = list()
for(i in 1:length(cells)){
  sub<-subset(NKs, subset= (ProgRec %in% c("R", "NR") & subcluster ==cells[i] ))
  cellname<-cells[i]
  for(m in 1:length(names)){
    datalist[[paste0(names[m],"_", cellname, sep="")]]<- modulescore_statistics_subcluster(m = m)
  }}

big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/S6F_NKs RvsNR Violin_scores_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)

#S6G metabolic comparison
library(stats)
a=read.table('metabolism.term.txt',sep='\t')
b=read.table('KEGG.txt',sep='\t')
colnames(a)=c('ID','Term','Category')
colnames(b)=c('ID','Term','GeneID','Gene')
term=unique(a$ID)
uniqueb<-unique(b$ID)
uniqueterm<-a[a$ID %in% uniqueb,]$ID
cells_suba<-c("CD56low NK", "CD56high NK")

suba<-subset(NKs, subset=subcluster %in% cells_suba)
suba$subcluster<-droplevels(suba$subcluster)
Idents(suba)<-suba$subcluster
DefaultAssay(suba)='RNA'
for(i in 1:length(uniqueterm)){
  su=subset(b,ID==as.character(uniqueterm[i]))
  suba=AddModuleScore(object = suba,features = list(su$Gene),name=uniqueterm[i],assay = 'RNA',ctrl = 10)
}

score=suba@meta.data
types=paste(uniqueterm,'1',sep='')
cells=unique(score$subcluster)
ps=matrix(nrow=0,ncol=5)
for(i in c(1:length(cells))){tryCatch({
  tmp=subset(score,subcluster==cells[i])
  h=subset(tmp,ProgRec=="NR")
  m=subset(tmp,ProgRec=='R')
  for(j in 1:length(types)){
    pv1=t.test(m[,types[j]],h[,types[j]])
    pvd1=data.frame(cells[i],types[j],pv1$p.value,(pv1$estimate[1]-pv1$estimate[2]))
    pvd1$Compare='R'
    colnames(pvd1)=c('cells','types','pvalue','diff','Compare')
    pvd=rbind(pvd1)
    ps=rbind(ps,pvd)
  }
},error=function(e){})}

ps$group[ps$diff>0]='up'
ps$group[ps$diff<0]='down'
ps$group[ps$diff==0]='none'
ps=subset(ps,pvalue<0.05)
ps$groups[ps$diff>0]=1
ps$groups[ps$diff<0]=-1
ps$groups[ps$diff==0]=0
ps$logp=log(-log10(ps$pvalue))
ps$logps=(ps$groups)*(ps$logp)
ps$ID=substr(ps$types,1,8)
psall=merge(ps,a,by='ID')

svg('S6G score.point_R_NR_NKs_flipped.svg',width=8,height=10)
p2=ggplot(psall,aes(x=Term,y=cells,fill=logps,shape=Compare))+geom_point(stroke=0.3,size=3.5)+scale_fill_gradient2(low='navy',mid='white',high='#8B0000',midpoint=0)+coord_flip()
p2=p2+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(color='black'))+scale_shape_manual(values=c(21,22,24))+facet_grid(Category~group, scales = "free", space = "free")
p2=p2+theme(axis.text.x = element_text(angle=90,vjust=1,hjust=1, size = fontsize), axis.text = element_text(size=fontsize), legend.text = element_text(size=fontsize), strip.text = element_text(size=fontsize))+labs(x='',y='')+scale_color_manual(values=c('black','black'))
p2$labels$fill<-"-log10(P)"
print(p2)
dev.off()

#S6E miloR
sub<-NKs
Idents(sub)<-sub$subcluster
DefaultAssay(sub)<-"RNA"
MILO  <- as.SingleCellExperiment(sub)
reducedDim(MILO, "PCA", withDimnames=TRUE) <- sub[['pca']]@cell.embeddings
reducedDim(MILO, "UMAP", withDimnames=TRUE) <- sub[['umap']]@cell.embeddings
MILO_obj <- Milo(MILO)
MILO_obj <- buildGraph(MILO_obj, k = 30, d = 10, reduced.dim = "PCA")
MILO_obj <- makeNhoods(MILO_obj, prop = 0.2, k = 30, d=10, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(MILO_obj) ##need distribution peak betwwen 50 and 100, otherwise up k and lower prop
MILO_obj <- countCells(MILO_obj, meta.data = data.frame(colData(MILO_obj)), sample="orig.ident")
MILO_obj <- calcNhoodDistance(MILO_obj, d=10, reduced.dim = "PCA") #~5min
milo.design <- as.data.frame(xtabs(~disease + orig.ident, data=data.frame(colData(MILO_obj))))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$orig.ident
milo.design <- milo.design[colnames(nhoodCounts(MILO_obj)),]
milo.res <- testNhoods(MILO_obj, design=~disease, design.df=milo.design)
MILO_obj <- buildNhoodGraph(MILO_obj)
milo.res <- annotateNhoods(MILO_obj, milo.res, coldata_col = "subcluster")##change to metadata
milo.res$subcluster<-factor(x=milo.res$subcluster, levels = c("CD56low NK", "CD56high NK"))
milo_plot<-plotDAbeeswarm(milo.res, group.by = "subcluster", alpha = 0.7)+xlab("")+FontSize(x.text = fontsize, y.text = fontsize, x.title = fontsize, y.title = fontsize,                                                        legend.text=element_text(size=fontsize), legend.title=element_text(size=fontsize )) 
ggsave(plot=milo_plot, filename="S6E milo_plot_NKs_0_7.svg" ,height=5, width=6, units="in", dpi=320)
