


# ---------------------- TF ----------------------
tf_all <- read.csv("/Users/niexiner/Documents/Phd/STTT/data/results/20260724/DatabaseExtract_v_1.01.csv", header = T, row.names = 1)
tf_all <- tf_all %>% filter(Is.TF. == "Yes")
head(tf_all)
dim(tf_all)

tflist <- genelist %>% lapply(function(i){
  intersect(i, unique(tf_all$HGNC.symbol))
})
lengths(tflist)


max_len <- max(sapply(tflist, length))
tf.mat<- as.data.frame(lapply(tflist, function(x) {
  length(x) <- max_len
  return(x)
}))
colnames(tf.mat) <- gsub("\\.", " ", colnames(tf.mat))
write.csv(tf.mat, "Supplementary TFs of four groups.csv", row.names = TRUE, na = "")

# old.deg <- read.table("../../BMK_5_DEG_Analysis/BMK_1_All_DEG/All.DEG_final.xls", header = F)
# head(old.deg)
# intersect(old.tfs, old.deg$V2)
# old.tfs <- c("ZNF641", "KIAA1549", "CEBPA", "KLF15", "ATOH8", "FOSL2", "STAT5B",
#              "IKZF2", "HAND2", "TWIST1", "WT1", "SMAD9", "OSR2", "GCFC2", "ID4",
#              "STAT5A", "HLX", "ZNF704", "LBX2", "MITF", "CUX1", "PROX1", "BACH2",
#              "ZFP36L1", "DACH2", "BHLHE41", "BHLHE40", "FOSL1", "PKNOX2", "POU2F2",
#              "GLI1", "E2F7", "JDP2", "GATA6", "SP110", "MKX", "CENPA", "CDX1",
#              "ZBTB46", "MSC", "PITX1", "RELB", "CREB5", "ZFP69B", "MYBL1", "TOX",
#              "HMGA1", "HMGA2", "MYC", "ETV4", "FOXF1", "HES1", "HIF1A", "TFAP2C", "NR4A2")

# old.crs <- c("TRPC3", "CRLF1", "FZD3", "ADGRF4", "GPR162", "PTPRQ", #"ADGRL3", "RXFP1",
#              "NPY1R", "KDR", "PTPN5", "LPAR3", "GABARAPL1", "P2RY1", "NPR3", "TRPC6",
#              "ERBB3", "TGFBR1", "ITPR1", "TRPC4", "GPR155", "CNR1", "IL6R", "ROR2",
#              "S1PR3", "P2RY14", "PRLR", "RARRES1", "SCARA5", "INSR", "EPHA4", "GPRC5B",
#              "FZD4", "RGMA", "IGF1R", "PTCH1")

# heatmap: group1 + group3
annotation_colors.tmp <- list(
  Condition = c("NCs" = "blue", "DCs" = "red", "F2tDCs" = "green"),
  Group = c("Group 1" = "#E64B35", "Group 3" = "#00A087")
)
tf_ph <- pheatmap::pheatmap(
  mat = pdata[unique(unlist(tflist[c(1,3)])), ],#as.matrix(pdata[unique(unlist(tflist_13)), ]),
  scale = "row",
  clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
  clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
  clustering_method = "complete",   # 聚类方法：最长距离法
  cluster_rows = TRUE,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  legend_labels = FALSE,
  annotation_col = col_group_df,
  annotation_colors = annotation_colors.tmp,
  annotation_row = deg_groups_df,
  annotation_names_row = F,annotation_names_col = F,cellwidth = 20, cellheight = 7,treeheight_row = 0,
  show_rownames = T,            # 如果基因太多，建议隐藏行名
  show_colnames = T,             # 显示样本名
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30), # 经典红白蓝配色
  clustering_callback = callback, 
  fontsize = 6,#, fontsize_row = 8, fontsize_col = 8, 
  name = NULL,border_color = NA
)
tf_ph
ggsave(tf_ph, file = "The heatmap of tf in group1 and group3.pdf", width = 6, height = 15)


# tf_ph <- pheatmap::pheatmap(
#   mat = pdata[unique(unlist(tflist)), ],#as.matrix(pdata[unique(unlist(tflist_13)), ]), 
#   scale = "row",                    
#   clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
#   clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
#   clustering_method = "complete",   # 聚类方法：最长距离法
#   cluster_rows = FALSE,              # 开启基因聚类
#   cluster_cols = FALSE,              # 开启样本聚类
#   legend_labels = FALSE,
#   annotation_col = col_group_df, 
#   annotation_colors = annotation_colors,
#   annotation_row = deg_groups_df,
#   annotation_names_row = F,annotation_names_col = F,cellwidth = 20, cellheight = 6,
#   show_rownames = TRUE,            # 如果基因太多，建议隐藏行名
#   show_colnames = F,             # 显示样本名
#   angle_col = 0,
#   color = colorRampPalette(c("navy", "white", "firebrick3"))(30), # 经典红白蓝配色
#   clustering_callback = callback, fontsize_row = 6, fontsize = 12,name = NULL,border_color = NA
# )
# tf_ph
# ggsave(tf_ph, file = "The heatmap of tf in four groups.pdf", width = 6, height = 24)

# # ---------------------- membrane protein ----------------------
CCchildren <- as.list(GOCCCHILDREN)
MFchildren <- as.list(GOMFCHILDREN)
CCchildren[["GO:0016020"]]
filter(egolist[[1]]@result, ID %in% CCchildren[["GO:0016020"]])
# 
# old.crs %>% length
# genelist %>% lapply(function(i){
#   intersect(i, old.crs)
# })
# gene_ids <- bitr(old.crs,
#                  fromType = "SYMBOL",
#                  toType = "ENTREZID",
#                  OrgDb = org.Hs.eg.db)
# ego<-enrichGO(
#   gene=gene_ids$ENTREZID,
#   OrgDb=org.Hs.eg.db,
#   keyType="ENTREZID",
#   ont="all",
#   pAdjustMethod="none",
#   pvalueCutoff=0.05,
#   qvalueCutoff=0.2,
#   minGSSize = 10,
#   maxGSSize = 50000,
#   readable = TRUE
# )

cr_ph <- pheatmap::pheatmap(
  mat = pdata[old.crs, ],#as.matrix(pdata[unique(unlist(tflist_13)), ]), 
  scale = "row",                    
  clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
  clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
  clustering_method = "complete",   # 聚类方法：最长距离法
  cluster_rows = FALSE,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  legend_labels = FALSE,
  annotation_col = col_group_df, 
  annotation_colors = annotation_colors,
  annotation_row = deg_groups_df,
  annotation_names_row = F,annotation_names_col = F,cellwidth = 20, cellheight = 7,
  show_rownames = T,            # 如果基因太多，建议隐藏行名
  show_colnames = T,             # 显示样本名
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30), # 经典红白蓝配色
  clustering_callback = callback, 
  fontsize = 6,#, fontsize_row = 8, fontsize_col = 8, 
  name = NULL,border_color = NA
)
cr_ph
ggsave(cr_ph, file = "The heatmap of membrane protein genes.pdf", width = 6, height = 6)


old.crs %>% length
# tt <- sapply(egolist[[1]]@result$ID, function(i){
#   aa <- filter(egolist[[1]]@result, ID == i) %>% .$geneID
#   sum(!is.na(match(old.crs, unlist(strsplit(aa, "/")))))
# })
# 
# tt[tt > 30]
# filter(egolist[[1]]@result, ID %in% names(tt[tt > 30]))


# # ---------------------- Inflammatory ----------------------

library(msigdbr)

cytokine_genes <- sapply(c("TNFSF", "IL", "CCL", "CXCL"), function(i){
  ind <- grepl(paste0("^", i), unique(unlist(genelist)))
  unique(unlist(genelist))[ind]
}) %>% unlist %>% unique
cytokine_genes <- cytokine_genes[-3]
cytokine_genes <- deg_groups_df[cytokine_genes, , drop = F] %>% mutate(Cytokine = rownames(.)) %>% group_by(Group) %>% arrange(Group, .by_group = F) %>% .$Cytokine


# # heatmap: group1 + group3
# annotation_colors.tmp <- list(
#   Condition = c("NCs" = "blue", "DCs" = "red", "F2tDCs" = "green"),
#   Group = c("Group 1" = "#E64B35", "Group 3" = "#00A087")
# )
cytokine_ph <- pheatmap::pheatmap(
  mat = pdata[cytokine_genes, ],#as.matrix(pdata[unique(unlist(tflist_13)), ]), 
  scale = "row",                    
  clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
  clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
  clustering_method = "complete",   # 聚类方法：最长距离法
  cluster_rows = F,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  legend_labels = FALSE,
  annotation_col = col_group_df,
  annotation_colors = annotation_colors,
  annotation_row = deg_groups_df,
  annotation_names_row = F,annotation_names_col = F,cellwidth = 20, cellheight = 7,
  show_rownames = T,            # 如果基因太多，建议隐藏行名
  show_colnames = T,             # 显示样本名
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30), # 经典红白蓝配色
  clustering_callback = callback, 
  fontsize = 6,
  name = NULL,border_color = NA
)
cytokine_ph
ggsave(cytokine_ph, file = "The heatmap of cytokines.pdf", width = 6, height = 6)