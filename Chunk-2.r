

# ---------------------- Gene cluster ----------------------
# library(pheatmap)
# library(tidyr)
# library(ggplot2)
# library(ggh4x)
# library(dplyr)


# tpm
tpm <- sweep(fpkm[, -c(1:2)], 2, colSums(fpkm[, -c(1:2)]), "/") * 1e6
head(tpm)
tpm <- cbind(fpkm[, c(1:2)], tpm)
head(tpm)
saveRDS(tpm, file = "All_gene_tpm.rds")

# log
log_tpm <- cbind(tpm[, c(1:2)], log2(tpm[, -c(1:2)] + 1))
# mean
pdata <- log_tpm %>% dplyr::filter(ID %in% unique(unlist(gene_list))) %>% dplyr::select(-ID) 
head(pdata)
rownames(pdata) <- pdata$Symbol
pdata <- pdata[, -1]
head(pdata)



callback = function(hc, ...){
  library(dendsort)
  dendsort(hc)}

# callback = function(hc, mat) {
#   dend <- as.dendrogram(hc)
#   dend <- rev(dend)
#   as.hclust(dend)
# }
# 
# callback = function(hc, mat) {
#   w <- rowMeans(mat, na.rm = TRUE)
#   
#   dend <- reorder(
#     as.dendrogram(hc),
#     wts = w,
#     agglo.FUN = mean
#   )
#   
#   as.hclust(dend)
# }

p <- pheatmap(
  mat = pdata, 
  scale = "row",                    
  clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
  clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
  clustering_method = "complete",   # 聚类方法：最长距离法
  cluster_rows = TRUE,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  show_rownames = FALSE,            # 如果基因太多，建议隐藏行名
  show_colnames = FALSE,             # 显示样本名
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30), # 经典红白蓝配色
  # main = "DEGs Clustering Heatmap",  # 标题
  clustering_callback = callback,
  angle_col = 0,cutree_rows  = 4
)
p
row_tree <- p$tree_row
deg_groups <- cutree(row_tree, k = 4)
deg_groups_df <- data.frame(
  row.names = names(deg_groups),
  Group = factor(
    paste0("Group ", deg_groups),
    levels = paste0("Group ", 1:4))
)
levels(deg_groups_df$Group) <- c("Group 4", "Group 2", "Group 3", "Group 1")
deg_groups_df$Group <- as.factor(as.character(droplevels(deg_groups_df$Group)))

col_group_df <- data.frame(row.names = c("C4003C-1", "C4003C-2", "C4003C-3", "C4003M-1", "C4003M-2", "C4003M-3", "C4003MB-1", "C4003MB-2",
                                         "C4003MB-3"), Condition = factor(c("NCs", "NCs", "NCs", "DCs", "DCs", "DCs", "F2tDCs", "F2tDCs", "F2tDCs"), levels = c("NCs", "DCs", "F2tDCs")))
annotation_colors <- list(
  Condition = c("NCs" = "blue", "DCs" = "red", "F2tDCs" = "green"),
  Group = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488")
  )

p <- pheatmap::pheatmap(
  mat = pdata, 
  scale = "row",  
  clustering_distance_rows = "euclidean", # 纵向（基因）聚类距离度量：欧氏距离
  clustering_distance_cols = "euclidean", # 横向（样本）聚类距离度量
  clustering_method = "complete",   # 聚类方法：最长距离法
  cluster_rows = TRUE,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  annotation_col = col_group_df,
  annotation_row = deg_groups_df,
  annotation_colors = annotation_colors,
  annotation_names_row = F,annotation_names_col = F,cellwidth = 20, #cellheight = 6,
  show_rownames = FALSE,            # 如果基因太多，建议隐藏行名
  show_colnames = FALSE,             # 显示样本名
  breaks = seq(-1, 1, length.out = 101),
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 经典红白蓝配色
  # main = "DEGs Clustering Heatmap",  # 标题
  clustering_callback = callback,
  angle_col = 0,
  cutree_rows  = 4,
  fontsize = 12,
  # filename = "The heatmap of DEG groups.pdf",
  width = 8, height = 8,border_color = NA
)
p
ggsave(p, file = "The heatmap of DEG groups.pdf", width = 6, height = 6)


pdata_merge <- data.frame("NCs" = rowMeans(pdata[, c(1:3)]), "DCs" = rowMeans(pdata[, c(4:6)]), "F2tDCs" = rowMeans(pdata[, c(7:9)]))
pdata_merge <- pdata_merge[deg_groups_df %>% arrange(Group) %>% rownames(), ]
head(pdata_merge)
p_merge <- pheatmap::pheatmap(
  mat = pdata_merge, 
  scale = "row",  
  cluster_rows = FALSE,              # 开启基因聚类
  cluster_cols = FALSE,              # 开启样本聚类
  annotation_col = data.frame(Condition = c("NCs", "DCs", "F2tDCs"), row.names = c("NCs", "DCs", "F2tDCs")),
  annotation_row = deg_groups_df,
  annotation_colors = annotation_colors,
  annotation_names_row = F,annotation_names_col = F,
  show_rownames = F,            # 如果基因太多，建议隐藏行名
  show_colnames = FALSE,             # 显示样本名
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 经典红白蓝配色
  clustering_callback = callback,
  angle_col = 0,
  fontsize = 12,
  width = 8, height = 8,border_color = NA
)
ggsave(p_merge, file = "The heatmap of DEG groups merged.pdf", width = 5, height = 5)

# line
pdata_z <- t(scale(t(pdata))) %>% as.data.frame()
deg_dat <- merge(x = deg_groups_df %>% mutate(Symbol = rownames(.)), y = pdata_z %>% mutate(Symbol = rownames(.)), by = "Symbol") %>%
  pivot_longer(cols = -c("Symbol", "Group"), names_to = "Sample", values_to = "Value") %>% 
  mutate(Condition = case_when(grepl("C-", Sample)~"NCs",
                               grepl("M-", Sample)~"DCs",
                               grepl("MB-", Sample)~"F2tDCs")) %>%
  group_by(Group, Symbol, Condition) %>% summarise(Value = mean(Value), .groups = "drop") %>%
  as.data.frame() 
deg_dat$Condition <- factor(deg_dat$Condition, levels = c("NCs", "DCs", "F2tDCs"))
head(deg_dat)

# p1 <- ggplot(deg_dat,
#              aes(
#                x = Condition,
#                y = Value,
#                col = Group,
#                group = Group
#              )) +
#   stat_summary(fun = "mean", geom = "line", linewidth = 0.5) +
#   stat_summary(fun = "mean", geom = "point", size = 2) +
#   labs(x = NULL, y = "Mean scaled expression") +
#   
#   scale_color_manual(values = c( "Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488"))+
#   theme_classic() +
#   theme(axis.text = element_text(size = 14),
#         axis.title  = element_text(size = 14),
#         legend.position = "top", 
#         legend.title = element_blank(), 
#         legend.text = element_text(size = 10)
#         ) 
# p1
# ggsave(p1, file = "The expression pattern of four DEG groups.pdf", width = 5, height = 4)


p2 <- ggplot(deg_dat,
       aes(
         x = Condition,
         y = Value,
         col = Group,
         fill = Group,
         group = Group
       )) +
  stat_summary(
    fun.data = ggpubr::mean_sd,
    geom = "ribbon",
    alpha = 0.18,
    color = NA
  ) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5) +
  labs(x = NULL, y = "Mean scaled expression") +
  theme_test() + 
  scale_color_manual(values = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488"))+
  scale_fill_manual(values = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488"))+
  theme(
        legend.position = "none", legend.title = element_blank(), 
        strip.background = element_blank(),strip.text = element_text(face = "bold"),
        text = element_text(size = 12)) +
  facet_wrap( ~ Group, ncol = 2)
p2
ggsave(p2, file = "The expression pattern of four DEG groups split.pdf", width = 6, height = 6)
# ggsave(p2, file = "The expression pattern of four DEG groups split.tif", width = 5, height = 4, dpi = 300)

head(deg_dat)
deg_dat %>% 
  group_by(Group, Condition) %>% summarise(mean(Value)) %>% 
  tidyr::pivot_wider(id_cols = Group, names_from = Condition, values_from = `mean(Value)`) %>% 
  write.csv(file = "The table of mean scaled expression value.csv")


deg_group_info<- as.data.frame(lapply(genelist, function(x) {
  length(x) <- max(lengths(genelist))
  return(x)
}))
colnames(deg_group_info) <- gsub("\\.", " ", colnames(deg_group_info))
write.csv(deg_group_info, "Degs of four groups.csv", row.names = TRUE, na = "")