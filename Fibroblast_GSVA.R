############################Limma analysis for GSVA value##################################################

library(limma)
library(pheatmap)
library(ggplot2)
library("reshape2")

res <- readRDS("/home/data/human/ESCC/Linna-ESCC/Fibro/ESCC_Fibro_GSVA_MarkerGenes.rds")
metainfo <- readRDS("/home/data/human/ESCC/Linna-ESCC/Fig/Cellinfo_20190924/Cell.info/Fibroblast_info.rds")

metainfo=metainfo[row.names(metainfo) %in% colnames(res),]
metainfo=metainfo[colnames(res),]
metainfo$cluster=metainfo$celltype

res=res[,row.names(metainfo)]

##################################################
#####Calculate mean GSVA score for each group
##################################################
index=colnames(res)

result=NULL
TPM=t(apply(res,1,function(x){
  re<-NULL
  for (i in unique(metainfo$cluster)){
    tmp<-rownames(subset(metainfo,cluster==i))
    re<-c(re,mean(x[index%in%tmp]))
  }
  return(re)
}))
colnames(TPM)=unique(metainfo$cluster)

row.names(TPM)=gsub("HALLMARK_","",row.names(TPM))

setwd("/home/data/human/ESCC/Linna-ESCC/Fig/Figure_20191004/SI/S5/Fibro_GSVA/")
P1=pheatmap(TPM[,c("NMF","NAF1",  "NAF2", "CAF1","CAF2",  "CAF3",  "CAF4" )],
            cluster_cols = F,
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            border_color=c("white"),
            treeheight_col=5,
            treeheight_row=0,
            cellwidth = 15, 
            cellheight = 15,
            angle_col = "0",
            main="Hallmark pathways in fibroblast cells"
            
)

ggsave("GSVA_FIBRO_markergenes.pdf",P1,height = 11,width = 6)


######canculate P and T #######

change=metainfo$celltype

lev=unique(metainfo$celltype)
est=res

pw<-lapply(lev,function(i){
  design<-ifelse(change==i,'A','other')
  design <- model.matrix(~ 0+factor(design))
  head(design)
  colnames(design) <- c('A','ohter')
  head(design)
  contrast.matrix <- makeContrasts(A-ohter, levels=design)
  
  fit <- lmFit(est, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  sigPathways <- topTable(fit2,  number=Inf,  adjust="BH")
  sigPathways$pathway<-rownames(sigPathways)
  sigPathways$cluster<-i
  return(sigPathways)
})


temp=do.call(rbind.data.frame,pw)


##select pathway
temp1=temp[temp$adj.P.Val<0.05 & abs(temp$t)>45,]
temp2=temp1[,c( "pathway", "cluster","t")]
temp3=dcast(temp2,pathway~cluster)
temp3[is.na(temp3)]=0
row.names(temp3)=temp3$pathway
temp3=temp3[,-1]

pathway=row.names(temp3)
df1=temp[temp$pathway %in%pathway,c( "pathway", "cluster","t")]
df2=dcast(df1,pathway~cluster)
row.names(df2)=df2$pathway
df2=df2[,-1]
row.names(df2)=gsub("HALLMARK_","",row.names(df2))



P3=pheatmap(df2[,c("NMF","NAF1",  "NAF2", "CAF1","CAF2",  "CAF3",  "CAF4" ,  "VSMC","Pericyte")],
            cluster_cols = F,
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            cellwidth = 15, 
            cellheight = 15,
            border_color=c("white"),
            treeheight_col=5,
            treeheight_row=0,
            main="Hallmark pathways in fibroblast cells"
            
)
ggsave("GSVA_Fibroblast_markergenes_tvalue_selected.pdf",P3,height = 8,width = 8.5)


P4=pheatmap(df2[,c("NMF","NAF1",  "NAF2", "CAF1","CAF2",  "CAF3",  "CAF4" ,  "VSMC","Pericyte")], scale="row",
            cluster_cols = F,
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            cellwidth = 15, 
            cellheight = 15,
            border_color=c("white"),
            treeheight_col=5,
            treeheight_row=0,
            main="Hallmark pathways in fibroblast cells"
            
)
ggsave("GSVA_Fibroblast_markergenes_zscore_selected1.pdf",P4,height = 8,width = 8.5)

order=P4$tree_row$order
df3=df2
df3$pathway=rownames(df3)
rownames(df3)=c(1:nrow(df3))
pathway=df3[order,'pathway']
edit(pathway)

pp=c("COMPLEMENT", "COAGULATION", "PEROXISOME","APOPTOSIS", "XENOBIOTIC_METABOLISM", "IL2_STAT5_SIGNALING", 
     "IL6_JAK_STAT3_SIGNALING",  "TGF_BETA_SIGNALING", 
     "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "INFLAMMATORY_RESPONSE", 
     "UNFOLDED_PROTEIN_RESPONSE", "MTORC1_SIGNALING","PI3K_AKT_MTOR_SIGNALING","GLYCOLYSIS", "APICAL_JUNCTION", "EPITHELIAL_MESENCHYMAL_TRANSITION", 
     "ANGIOGENESIS", "OXIDATIVE_PHOSPHORYLATION", "HEDGEHOG_SIGNALING")



PP=pheatmap(df2[pp,c("NMF","NAF1",  "NAF2", "CAF1","CAF2",  "CAF3",  "CAF4" ,  "VSMC","Pericyte")], scale="row",
            cluster_cols = F,cluster_rows  = F,
            color = colorRampPalette(c("#004A72", "#027FB6","#01AECA", "#7ED2E3", 
                                       "#D1E7F2", "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50),
            fontsize = 12,
            fontsize_col = 12,
            fontsize_row =12,
            cellwidth = 15, 
            cellheight = 15,
            border_color=c("white"),
            treeheight_col=5,
            treeheight_row=0,
            main="Hallmark pathways in fibroblast cells"
            
)

ggsave("GSVA_Fibroblast_markergenes_zscore_selected2.pdf",PP,height = 8,width = 8.5)


