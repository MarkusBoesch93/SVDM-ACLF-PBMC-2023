#Reproducing Figure 3 and Supplementary Figure 3
#Created by Markus Boesch
#loaded in all necessary packages and colours from scRNA_prep
#subset the Monocytes
#ProgRec contains ACLF-R (R) and ACLF-NR (NR)

subc=subset(Monocytes, subset= (ProgRec %in% c("R", "NR") & subcluster =="cMon" ))
Idents(subc)<-subc$ProgRec

#Scoring of GO pathways
Phagocytosis=list(c("ABCA1", "ABCA7", "ABL1", "ADGRB1", "ADIPOQ", "ADORA1", "ADORA2A", "AHSG", "AIF1", "ALOX15", "ANO6", "ANXA1", "ANXA11", "ANXA3", "APOA1", "APOA2", "APPL1", "APPL2", "ARHGAP12", "ARHGAP25", "ARL8B", "ATG3", "ATG5", "AXL", "AZU1", "BCR", "BECN1", "BIN2", "C1orf43", "C2", "C3", "C4A", "C4B", "C4BPA", "C4BPB", "CALR", "CAMK1D", "CCL2", "CCR2", "CD14", "CD300A", "CD300LF", "CD302", "CD36", "CD47", "CD93", "CDC42", "CDC42SE1", "CDC42SE2", "CEACAM4", "CEBPE", "CLCN3", "CLEC7A", "CLN3", "CNN2", "COLEC10", "COLEC11", "COLEC12", "CORO1A", "CORO1C", "CRP", "CRYBA1", "CSK", "CYBA", "DNM2", "DOCK1", "DOCK2", "DYSF", "EIF2AK1", "ELANE", "ELMO1", "ELMO2", "ELMO3", "F2RL1", "FCER1G", "FCGR1A", "FCGR2B", "FCN1", "FCN2", "FCN3", "FER1L5", "FGR", "FPR2", "FYN", "GAS6", "GATA2", "GSN", "GULP1", "HAVCR1", "HCK", "HMGB1", "ICAM3", "ICAM5", "IFNG", "IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGHV1-18", "IGHV1-24", "IGHV1-3", "IGHV1-45", "IGHV1-58", "IGHV1-69", "IGHV1-69-2", "IGHV1-69D", "IGHV1OR15-1", "IGHV2-26", "IGHV2-5", "IGHV2-70", "IGHV2-70D", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-16", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-35", "IGHV3-38", "IGHV3-43", "IGHV3-48", "IGHV3-49", "IGHV3-53", "IGHV3-64", "IGHV3-64D", "IGHV3-66", "IGHV3-7", "IGHV3-72", "IGHV3-73", "IGHV3-74", "IGHV4-28", "IGHV4-31", "IGHV4-34", "IGHV4-39", "IGHV4-4", "IGHV4-59", "IGHV4-61", "IGHV5-10-1", "IGHV5-51", "IGHV6-1", "IGHV7-4-1", "IGHV7-81", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7", "IGLL1", "IGLL5", "IL15", "IL15RA", "IL1B", "IL2RB", "IL2RG", "IRF8", "ITGA2", "ITGAL", "ITGAM", "ITGAV", "ITGB1", "ITGB2", "ITGB3", "JMJD6", "KIAA1109", "LBP", "LDLR", "LEP", "LEPR", "LIMK1", "LMAN2", "LRP1", "LYAR", "LYN", "LYST", "MARCO", "MBL2", "MEGF10", "MERTK", "MESD", "MET", "MFGE8", "MIR17", "MIR181B1", "MIR183", "MIR20A", "MSR1", "MST1R", "MYD88", "MYH9", "MYO18A", "MYO1G", "MYO7A", "NCF2", "NCF4", "NCKAP1L", "NOD2", "NR1H3", "OLFM4", "P2RX7", "P2RY6", "PAK1", "PEAR1", "PECAM1", "PIK3CA", "PIKFYVE", "PIP4P2", "PIP5K1A", "PIP5K1C", "PLA2G5", "PLA2G6", "PLCG2", "PLD2", "PLD4", "PLPP4", "PLSCR1", "PRKCD", "PRKCE", "PRKCG", "PRTN3", "PTK2", "PTPRC", "PTPRJ", "PTX3", "PYCARD", "RAB11FIP2", "RAB14", "RAB20", "RAB27A", "RAB31", "RAB34", "RAB39A", "RAB5A", "RAB7A", "RAB7B", "RAC1", "RAC2", "RACK1", "RARA", "RHOBTB1", "RHOBTB2", "RHOG", "RHOH", "RUBCN", "SCARB1", "SFTPA1", "SFTPD", "SH3BP1", "SIRPA", "SIRPB1", "SIRPG", "SLAMF1", "SLC11A1", "SNX3", "SOD1", "SPACA3", "SPG11", "SPHK1", "SPON2", "SRC", "SRPX", "STAP1", "SYK", "SYT11", "SYT7", "TAFA4", "TGM2", "THBS1", "TICAM2", "TLR2", "TLR4", "TM9SF4", "TMEM175", "TNF", "TRBC1", "TRBC2", "TRDC", "TREM2", "TREX1", "TUB", "TULP1", "TUSC2", "TXNDC5", "TYRO3", "TYROBP", "UNC13D", "VAMP7", "VAV1", "VAV2", "VAV3", "XKR4", "XKR5", "XKR6", "XKR7", "XKR8", "XKR9", "YES1"))#GO:0006909

antigen_pres=list(c("ABCB9", "ACE", "AP3B1", "AP3D1", "ARL8B", "ATG5", "AZGP1", "B2M", "CALR", "CCL19", "CCL21", "CCR7", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD209", "CD68", "CD74", "CD8A", "CLEC4A", "CLEC4M", "CTSD", "CTSE", "CTSF", "CTSH", "CTSL", "CTSS", "CTSV", "DNM2", "ERAP1", "ERAP2", "EXT1", "FCER1G", "FCGR2B", "FGL2", "GBA", "HFE", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "ICAM1", "IDE", "IFI30", "IKBKB", "KDM5D", "LGMN", "LILRB2", "LNPEP", "MARCHF1", "MARCHF8", "MFSD6", "MR1", "NOD1", "NOD2", "PDIA3", "PIKFYVE", "PSMB8", "PSME1", "PYCARD", "RAB10", "RAB27A", "RAB32", "RAB33A", "RAB34", "RAB35", "RAB3B", "RAB3C", "RAB4A", "RAB5B", "RAB6A", "RAB8B", "RELB", "RFTN1", "SAR1B", "SLC11A1", "TAP1", "TAP2", "TAPBP", "TAPBPL", "THBS1", "TRAF6", "TREM2", "TREX1", "WAS", "WDFY4", "YTHDF1"))#GO:0019882

migration =list(c('PECAM1','GAS6','IL4','PLCB1','APP','APOD','PTK2B','HMGB1','TMEM102','GATA3','CCL17','C1QBP','GCSAM','CD99L2','MSMP','ARHGEF5','JAM2','MSTN','TNFSF11','CXCR4','XG','ANXA1','CCL7','CCL8','DUSP1','ZAP70','S100A14','SLC8B1','MAPK1','C5AR1','CKLF','CCL4','PDGFB','PTK2','CCL3L1','CCL22','CMKLR1','CCN3','PIK3CG','TNFRSF14','AKT1','CCL19','KLRK1','MDK','DEFB104A','CXCL17','TNF','EPS8','XCL1','S1PR1','RARRES2','CCL21','SLC12A2','AIF1','FPR2','CXCR2','CXCR1','ITGA4','CCL18','CCL23','ADAM10','CALCA','CCL13','ANO6','CCR2','IL12A','IL34','WNK1','CD99','ARTN','MYO1G','PDGFD','FADD','AKIRIN1','CXCL10','PLA2G7','GREM1','PADI2','CDC42','THBS1','RET','LYN','ICAM1','FLT1','LGALS3','CSF1','CRKL','CRK','CD200R1','TRPM4','KARS1','NBL1','LGALS9','MADCAM1','CCL16','CCL25','RHOA','WASL','GPR15L','JAML','LGMN','IL27RA','CCL24','WNT5A','SPNS2','CD200','CX3CL1','MSN','ITGB7','CXCL11','ALOX5','PIK3CD','SLIT2','STK10','CCR6','CCR5','CCL11','S100A7','AGER','CXCL12','TRPM2','ADAM8','SIRPA','CCL20','ADAM17','SAA1','SERPINE1','TBX21','CCL14','ECM1','RIPOR2','HSD3B7','CCL1','CXCL16','F11R','CD47','PYCARD','IL6R','CCL5','CCL2','S100A12','GPR15','CXCR3','FUT7','CXCL14','CXCL13','AIRE','OXSR1','SPN','ITGAL','SELENOK','CRTAM','CREB3','CCL26','CSF1R','MAPK14','PLG','TNFRSF11A','RPS19','CH25H','DEFB124','CCL4L1','BMP5','SLAMF8','CCL15','XCL2','ADTRP','IL6','CYP7B1','MAPK3','CCL27','RIPK3','MOSPD2','MIA3','GPR183','CCR7','CCR1','STK39','GCSAML','DEFA1','CALR','NLRP12','C3AR1','LRCH1','CCL3','DOCK8','PTPRO'))

Cytokines=list(c("CCL2", "CCL3", "CCL4", "CCL5", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL8", "CXCL9", "FGF2", "IFNG", "IL10", "IL12A", "IL15RA", "IL18", "IL1B", "IL1RN", "IL2", "IL4", "IL6", "IL7", "LTA", "PDGFB", "TNF", "VEGFA"))

S100= list(c('S100A1','S100A2','S100A3','S100A4','S100A5','S100A6','S100A7','S100A7A','S100A7L2','S100A7P1','S100A7P2','S100A8','S100A9','S100A10','S100A11','S100A12','S100A13','S100A14','S100A15A','S100A16','S100B','S100G','S100P','S100Z'))

IFNresponse =list(c('ADAR','APOBEC3','BST2','CD74','MB21D1','DDIT4','DDX58','DDX60','EIF2AK2','GBP1','GBP2','HPSE','IFI44L','IFI6','IFIH1','IFIT1','IRF1','IRF7','ISG15','ISG20','MAP3K14','MOV10','MS4A4A','MX1','MX2','NAMPT','NT5C3','OAS1','OAS2','OAS3','OASL','P2RY6','PHF15','PML','RSAD2','RTP4','SLC15A3','SLC25A28','SSBP3','TREX1','TRIM5','TRIM25','SUN2','ZC3HAV1','IFITM1','IFITM2','IFITM3'))

Inflammation = list(c('IFNG', 'IL10', 'IL12A', 'IL13', 'IL17A', 'IL18','IL1A', 'IL1B', 'IL2', 'IL21', 'IL22', 'IL23A', 'IL4', 'IL5', 'IL6','TNF','CXCL8'))
MHC = list(c('HLA-DMA','HLA-DMB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB5'))
Activation_M=list(c("ADAM10", "ADAM9", "AZU1", "CD33", "CSF1", "DYSF", "FER1L5", "FN1", "FOXP1", "HYAL2", "LILRB4", "MT1G", "NPY", "SPACA3"))#GO:0042117
Differentiation_M=list(c("ACIN1", "APCS", "BMP4", "CD4", "CD74", "CDK6", "CSF1", "CSF1R", "CSF2", "CTNNBIP1", "DCSTAMP", "FASN", "FES", "FOXP1", "GPR68", "HLA-DRB1", "HOXA7", "IFI16", "IL31RA", "IL34", "INPP5D", "IRF7", "JUN", "MED1", "MIR125B1", "MT1G", "MYC", "MYH9", "PDE1B", "PIR", "PPARG", "SP3", "THOC5", "VEGFA", "ZBTB46", "ZFP36L1"))#GO:0030224

types=c(Phagocytosis,antigen_pres ,Cytokines,migration, IFNresponse, S100,Differentiation_M)
names=c("Phagocytosis","Ag_Presentation",  "Cytokine", "Migration", "IFN_response", "S100_Family","Differentiation_M")

for(m in 1:length(names)){ tryCatch({
  subc=AddModuleScore(object = subc, features = types[m],name = names[m], assay = "RNA",ctrl = 10)
  Figure<-VlnPlot(subc, group.by= "ProgRec",features = paste(names[m],"1", sep = "") , assay = "RNA", cols = c("sienna3"  ,  "red4" ), pt.size=0)+theme(axis.title = element_blank(), axis.text.x = element_text(size = fontsize), plot.title = element_blank())+NoLegend()
  ggsave(plot=Figure, filename=paste(getwd(),"/3A_",names[m], "_ProgRec_Violin_score.svg", sep = "") ,height=2, width=2, units="in", dpi=320)
},error=function(e){})}

#stats for GO scoring
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
library(limma)
datalist = list()
for(m in 1:length(names)){
  datalist[[names[m]]]<- modulescore_statistics_subcluster(m = m)
}

big_data = do.call(rbind, datalist)
openxlsx::write.xlsx(big_data, file=paste(getwd(),"/3A_cMon score RvsNR Violin_stats.xlsx",sep = ""), col.names=TRUE, row.names=TRUE)


##3C gProfiler
library(gprofiler2)
library(openxlsx)
Idents(RandNRcMon)<-RandNRcMon$ProgRec
genes <-FindAllMarkers(RandNRcMon, min.pct = 0.25 , test.use="MAST", latent.vars = c("nCount_RNA","percent.mito"), assay = "RNA")
genes2 <- genes[genes$p_val_adj<0.05,]
genes3 <- genes2[genes2$avg_log2FC >0.25,]
genes3 <- genes3[order(genes3$avg_log2FC, decreasing = T),]
R_go<-genes3[which(genes3$cluster=="R"),]$gene 
NR_go<-genes3[which(genes3$cluster=="NR"),]$gene 

sub_gprofiler_genes<-qpcR:::cbind.na(R_go, NR_go)
wb<-openxlsx::createWorkbook()
addWorksheet(wb, sheetName = "genes")
writeData(wb=wb,sheet="genes", x=sub_gprofiler_genes)
gostR_go<-gost(query=R_go, organism = "hsapiens", sources="GO:BP", evcodes = T,ordered_query=T)
gostNR_go<-gost(query=NR_go, organism = "hsapiens", sources="GO:BP", evcodes = T, ordered_query = T)

addWorksheet(wb, sheetName = "R_go")
writeData(wb=wb,sheet="R_go", x=as.data.frame(gostR_go$result))
addWorksheet(wb, sheetName = "NR_go")
writeData(wb=wb,sheet="NR_go", x=as.data.frame(gostNR_go$result))

openxlsx::saveWorkbook(wb, paste(getwd(),"/3C_cMon_R_vs_NR_gprofiler_genes.xlsx",sep = ""), overwrite = T)

#3H
library(stats)
a=read.table('metabolism.term.txt',sep='\t')
b=read.table('KEGG.txt',sep='\t')
colnames(a)=c('ID','Term','Category')
colnames(b)=c('ID','Term','GeneID','Gene')
term=unique(a$ID)
uniqueb<-unique(b$ID)
uniqueterm<-a[a$ID %in% uniqueb,]$ID

cells_suba<-c("cMon")
suba<-subset(Monocytes, subset=subcluster %in% cells_suba)
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

svg('3I score.point_R_NR_cMon_flipped.svg',width=8,height=10)
p2=ggplot(psall,aes(x=Term,y=cells,fill=logps,shape=Compare))+geom_point(stroke=0.3,size=3.5)+scale_fill_gradient2(low='navy',mid='white',high='#8B0000',midpoint=0)+coord_flip()
p2=p2+theme_bw()+theme(panel.grid = element_blank(),axis.text = element_text(color='black'))+scale_shape_manual(values=c(21,22,24))+facet_grid(Category~group, scales = "free", space = "free")
p2=p2+theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1, size = fontsize), axis.text = element_text(size=fontsize), legend.text = element_text(size=fontsize), strip.text = element_text(size=fontsize))+labs(x='',y='')+scale_color_manual(values=c('black','black'))
p2$labels$fill<-"-log10(P)"
print(p2)
dev.off()
