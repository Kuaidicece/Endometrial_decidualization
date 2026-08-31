library(dplyr)
library(DESeq2)
library(openxlsx)
library(VennDiagram)
library(grid)
library(GseaVis)
source("modified_function.R")
library(ggplot2)
library(pheatmap)
library(ggh4x)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
seed <- 2026

# ---------------------- read matrix ----------------------
mat <- read.table("All_gene_counts.list", header = TRUE, sep = "\t", comment.char = "", check.names = F)
colnames(mat)[1] <- "ID"
head(mat)

counts <- mat
rownames(counts) <- counts$ID
counts <- counts[, -c(1, 2)]
head(counts)
dim(counts)

# filter genes
counts <- counts[rowSums(counts) != 0, ]
dim(counts)

# ---------------------- DESeq2 ----------------------
# sample information
sample_info <- data.frame(
  row.names = c("C4003C-1", "C4003C-2", "C4003C-3", "C4003M-1", "C4003M-2", "C4003M-3", "C4003MB-1", "C4003MB-2", "C4003MB-3"),
  condition = factor(c("NCs", "NCs", "NCs", "DCs", "DCs", "DCs", "F2tDCs", "F2tDCs", "F2tDCs"))
)
# sample_info

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~condition)

# run deg analysis
dds <- DESeq(dds)

# extract NCs vs DCs
res_NCs_vs_DCs <- results(dds, contrast = c("condition", "NCs", "DCs"))
res_NCs_vs_DCs <- as.data.frame(res_NCs_vs_DCs)
head(res_NCs_vs_DCs)

# extract DCs vs DCs+FGF2
res_DCs_vs_DCs_FGF2 <- results(dds, contrast = c("condition", "DCs", "F2tDCs"))
res_DCs_vs_DCs_FGF2 <- as.data.frame(res_DCs_vs_DCs_FGF2)

# extract DCs vs NCs+FGF2
res_NCs_vs_DCs_FGF2 <- results(dds, contrast = c("condition", "NCs", "F2tDCs"))
res_NCs_vs_DCs_FGF2 <- as.data.frame(res_NCs_vs_DCs_FGF2)

deg <- list(res_NCs_vs_DCs, res_DCs_vs_DCs_FGF2, res_NCs_vs_DCs_FGF2) %>%
  lapply(function(i) {
    i %>%
      mutate(ID = row.names(.), tendency = case_when(
        padj < 0.01 & log2FoldChange > 1 ~ "Up",
        padj < 0.01 & log2FoldChange < -1 ~ "Down",
        TRUE ~ "ns"
      ))
  })
names(deg) <- c("NCs_vs_DCs", "DCs_vs_F2tDCs", "NCs_vs_F2tDCs")
deg <- lapply(deg, function(i) {
  merge(x = mat[, 1:2], y = i, by = "ID") %>%
    arrange(tendency)
})

# save
write.xlsx(
  deg,
  file = "The comparison of all genes among NCs, DCs and F2tDCs.xlsx",
  asTable = FALSE,
  rowNames = FALSE,
  keepNA = TRUE,
  firstRow = FALSE
)


# visualization Venn
gene_list <- lapply(deg, function(i) {
  filter(i, tendency != "ns") %>% .$ID
})
venn_obj <- venn.diagram(
  x = gene_list,
  filename = NULL,
  fill = c("#88a0dc", "#ed968c", "#f9d14a"),
  alpha = c(0.5, 0.5, 0.5),
  col = c("#88a0dc", "#ed968c", "#f9d14a"),
  category.names = c("NCs vs DCs", "DCs vs F2tDCs", "NCs vs F2tDCs"),
  cex = 1.2, cat.cex = 1.2, cat.fontface = "bold", cat.dist = c(0.1, 0.1, 0.05), margin = 0.1
)

pdf("Venn_of_DEGs.pdf", width = 6, height = 6)
grid.newpage()
grid.draw(venn_obj)
dev.off()


deg_genes <- unique(unlist(gene_list))
# read fpkm
fpkm <- read.table("All_gene_fpkm.list",
  header = TRUE, sep = "\t", comment.char = "", check.names = F
)
colnames(fpkm)[1] <- "ID"
head(fpkm)

dat <- filter(fpkm, ID %in% deg_genes) %>% select(-ID)
colnames(dat) <- c("Symbol name", "NCs-1", "NCs-2", "NCs-3", "DCs-1", "DCs-2", "DCs-3", "F2tDCs-1", "F2tDCs-2", "F2tDCs-3")
head(dat)

# save
write.xlsx(list("all deg" = dat),
  file = "The fpkm of DEGs among NCs, DCs and F2tDCs.xlsx", asTable = FALSE, rowNames = FALSE, keepNA = TRUE,
  firstRow = FALSE
)


# ---------------------- GSEA ----------------------

make_geneList <- function(res) {
  res$SYMBOL <- mat[match(rownames(res), mat$ID), "Symbol"]
  gene_ids <- bitr(res$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  res_m <- merge(gene_ids, res, by = "SYMBOL")
  res_m <- res_m[order(abs(res_m$stat), decreasing = TRUE), ]
  x <- res_m$stat
  names(x) <- res_m$ENTREZID

  x <- x[!is.na(x)]
  x <- x[!duplicated(names(x))]
  sort(x, decreasing = TRUE)
}
# res_NCs_vs_DCs, res_DCs_vs_DCs_FGF2, res_NCs_vs_DCs_FGF2
geneList_NCs_vs_DCs <- make_geneList(res_NCs_vs_DCs)
geneList_DCs_vs_DCs_FGF2 <- make_geneList(res_DCs_vs_DCs_FGF2)
geneList_NCs_vs_DCs_FGF2 <- make_geneList(res_NCs_vs_DCs_FGF2)

all_glist <- list(geneList_NCs_vs_DCs, geneList_DCs_vs_DCs_FGF2, geneList_NCs_vs_DCs_FGF2)

# loop to enrich
m_gsea_list <- lapply(all_glist, function(x) {
  gsea_kegg <- clusterProfiler::gseKEGG(
    geneList = x,
    organism = "hsa",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = FALSE,
    eps = 0,
    seed = TRUE
  )
  gsea_kegg
})
names(m_gsea_list) <- c("NCs_vs_DCs", "DCs_vs_F2tDCs", "NCs_vs_F2tDCs")


geneset.info <- m_gsea_list[[1]] %>%
  data.frame() %>%
  filter(Description %in% c(
    "FoxO signaling pathway",
    "MAPK signaling pathway",
    "JAK-STAT signaling pathway",
    "PI3K-Akt signaling pathway",
    "cAMP signaling pathway"
  )) %>%
  dplyr::select(Description)
geneset.info
geneset.id <- rownames(geneset.info)
names(geneset.id) <- geneset.info$Description
geneset.id

for (i in 1:length(geneset.id)) {
  gsea_plot <- GSEAmultiGP_modified(
    gsea_list = m_gsea_list,
    geneSetID = geneset.id[[i]],
    exp_name = gsub("_", " ", names(m_gsea_list)),
    curve.col = c("#88a0dc", "#ed968c", "#f9d14a"),
    addPval = T, pvalX = 0.01, pvalY = 0.01, pDigit = 3, base_size = 12
  )
  # gsea_plot
  ggsave(gsea_plot, file = paste0("./GSEA/The GSEA plot of ", names(geneset.id)[[i]], ".pdf"), width = 8, height = 6, create.dir = T)
  print(i)
}


# ---------------------- Gene cluster ----------------------
# tpm
tpm <- sweep(fpkm[, -c(1:2)], 2, colSums(fpkm[, -c(1:2)]), "/") * 1e6
head(tpm)
tpm <- cbind(fpkm[, c(1:2)], tpm)
head(tpm)
saveRDS(tpm, file = "All_gene_tpm.rds")

# log
log_tpm <- cbind(tpm[, c(1:2)], log2(tpm[, -c(1:2)] + 1))
# mean
pdata <- log_tpm %>%
  dplyr::filter(ID %in% unique(unlist(gene_list))) %>%
  dplyr::select(-ID)
head(pdata)
rownames(pdata) <- pdata$Symbol
pdata <- pdata[, -1]
head(pdata)


callback <- function(hc, ...) {
  library(dendsort)
  dendsort(hc)
}

p <- pheatmap(
  mat = pdata,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
  # main = "DEGs Clustering Heatmap",
  clustering_callback = callback,
  angle_col = 0, cutree_rows = 4
)
p
row_tree <- p$tree_row
deg_groups <- cutree(row_tree, k = 4)
deg_groups_df <- data.frame(
  row.names = names(deg_groups),
  Group = factor(
    paste0("Group ", deg_groups),
    levels = paste0("Group ", 1:4)
  )
)
levels(deg_groups_df$Group) <- c("Group 4", "Group 2", "Group 3", "Group 1")
deg_groups_df$Group <- as.factor(as.character(droplevels(deg_groups_df$Group)))

col_group_df <- data.frame(row.names = c(
  "C4003C-1", "C4003C-2", "C4003C-3", "C4003M-1", "C4003M-2", "C4003M-3", "C4003MB-1", "C4003MB-2",
  "C4003MB-3"
), Condition = factor(c("NCs", "NCs", "NCs", "DCs", "DCs", "DCs", "F2tDCs", "F2tDCs", "F2tDCs"), levels = c("NCs", "DCs", "F2tDCs")))
annotation_colors <- list(
  Condition = c("NCs" = "blue", "DCs" = "red", "F2tDCs" = "green"),
  Group = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488")
)

p <- pheatmap::pheatmap(
  mat = pdata,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_group_df,
  annotation_row = deg_groups_df,
  annotation_colors = annotation_colors,
  annotation_names_row = F, annotation_names_col = F, cellwidth = 20, # cellheight = 6,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(-1, 1, length.out = 101),
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  # main = "DEGs Clustering Heatmap",
  clustering_callback = callback,
  angle_col = 0,
  cutree_rows = 4,
  fontsize = 12,
  # filename = "The heatmap of DEG groups.pdf",
  width = 8, height = 8, border_color = NA
)
p
ggsave(p, file = "The heatmap of DEG groups.pdf", width = 6, height = 6)


pdata_merge <- data.frame("NCs" = rowMeans(pdata[, c(1:3)]), "DCs" = rowMeans(pdata[, c(4:6)]), "F2tDCs" = rowMeans(pdata[, c(7:9)]))
pdata_merge <- pdata_merge[deg_groups_df %>%
  arrange(Group) %>%
  rownames(), ]
head(pdata_merge)
p_merge <- pheatmap::pheatmap(
  mat = pdata_merge,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = data.frame(Condition = c("NCs", "DCs", "F2tDCs"), row.names = c("NCs", "DCs", "F2tDCs")),
  annotation_row = deg_groups_df,
  annotation_colors = annotation_colors,
  annotation_names_row = F, annotation_names_col = F,
  show_rownames = F,
  show_colnames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  clustering_callback = callback,
  angle_col = 0,
  fontsize = 12,
  width = 8, height = 8, border_color = NA
)
ggsave(p_merge, file = "The heatmap of DEG groups merged.pdf", width = 5, height = 5)

# line
pdata_z <- t(scale(t(pdata))) %>% as.data.frame()
deg_dat <- merge(x = deg_groups_df %>% mutate(Symbol = rownames(.)), y = pdata_z %>% mutate(Symbol = rownames(.)), by = "Symbol") %>%
  pivot_longer(cols = -c("Symbol", "Group"), names_to = "Sample", values_to = "Value") %>%
  mutate(Condition = case_when(
    grepl("C-", Sample) ~ "NCs",
    grepl("M-", Sample) ~ "DCs",
    grepl("MB-", Sample) ~ "F2tDCs"
  )) %>%
  group_by(Group, Symbol, Condition) %>%
  summarise(Value = mean(Value), .groups = "drop") %>%
  as.data.frame()
deg_dat$Condition <- factor(deg_dat$Condition, levels = c("NCs", "DCs", "F2tDCs"))
head(deg_dat)


p2 <- ggplot(
  deg_dat,
  aes(
    x = Condition,
    y = Value,
    col = Group,
    fill = Group,
    group = Group
  )
) +
  stat_summary(
    fun.data = ggpubr::mean_sd,
    geom = "ribbon",
    alpha = 0.18,
    color = NA
  ) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5) +
  labs(x = NULL, y = "Mean scaled expression") +
  theme_test() +
  scale_color_manual(values = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488")) +
  scale_fill_manual(values = c("Group 1" = "#E64B35", "Group 2" = "#4DBBD5", "Group 3" = "#00A087", "Group 4" = "#3C5488")) +
  theme(
    legend.position = "none", legend.title = element_blank(),
    strip.background = element_blank(), strip.text = element_text(face = "bold"),
    text = element_text(size = 12)
  ) +
  facet_wrap(~Group, ncol = 2)
p2
ggsave(p2, file = "The expression pattern of four DEG groups split.pdf", width = 6, height = 6)

head(deg_dat)
deg_dat %>%
  group_by(Group, Condition) %>%
  summarise(mean(Value)) %>%
  tidyr::pivot_wider(id_cols = Group, names_from = Condition, values_from = `mean(Value)`) %>%
  write.csv(file = "The table of mean scaled expression value.csv")


deg_group_info <- as.data.frame(lapply(genelist, function(x) {
  length(x) <- max(lengths(genelist))
  return(x)
}))
colnames(deg_group_info) <- gsub("\\.", " ", colnames(deg_group_info))
write.csv(deg_group_info, "Degs of four groups.csv", row.names = TRUE, na = "")


# ---------------------- GO ----------------------
genelist <- split(deg_dat$Symbol, deg_dat$Group) %>% lapply(unique)

dir.create("./Enrich/Go/", showWarnings = F, recursive = T)


egolist <- list()
dir.create("./Enrich/Go", showWarnings = F, recursive = T)
for (i in names(genelist)) {
  genes <- genelist[[i]]
  # genes2 <- filter(mat, Symbol %in% genes) %>% .$ID
  print(paste0(i, " : ", length(genes)))
  gene_ids <- bitr(genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  ego <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "all",
    pAdjustMethod = "none",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 50000,
    readable = TRUE
  )

  write.csv(ego@result, file = paste0("./Enrich/Go/Go_", gsub(" ", "_", i), "_results.csv"))
  egolist[[i]] <- ego
}


BPancestor <- as.list(GOBPANCESTOR)
BPancestor[["GO:0007169"]]
BPancestor[["GO:0035556"]]

select(GO.db,
  keys = head(BPancestor[["GO:0007169"]], -1),
  columns = c("TERM", "ONTOLOGY")
)
goids <- keys(GO.db, keytype = "GOID")
goids_filtered <- do.call(rbind, lapply(goids, function(x) {
  if (length(BPancestor[[x]]) >= 11 && length(BPancestor[[x]]) <= 13) {
    return(select(GO.db, keys = x, columns = c("TERM")))
  }
}))

top_term_1 <- filter(egolist[[1]]@result, ONTOLOGY == "BP") %>% top_n(n = 200, wt = Count)
top_term_1 <- top_term_1[!is.na(match(rownames(top_term_1), goids_filtered$GOID)), ]
dim(top_term_1)

top_term_dot_1 <- ggplot(top_term_1, aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T))) +
  geom_point(aes(size = Count, color = pvalue)) +
  # scale_color_gradientn(colors = c("#BF1E27", "#FEB466","#F9FCCB","#6296C5","#38489D"))+
  scale_color_gradientn(colors = rev(c("navy", "white", "firebrick3"))) +
  guides(color = guide_colorbar(reverse = T)) +
  labs(x = "-log10(Pvalue)", y = NULL, color = "Pvalue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(lineheight = 0.8), legend.position = "right",
    text = element_text(size = 12), # panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10), legend.title = element_text(size = 10)
  )
top_term_dot_1
ggsave(top_term_dot_1, file = "The dotplot of significant GO BP in Group1.pdf", width = 6, height = 6)


top_term_3 <- filter(egolist[[3]]@result, ONTOLOGY == "BP") %>% top_n(n = 150, wt = Count)
top_term_3 <- top_term_3[!is.na(match(rownames(top_term_3), goids_filtered$GOID)), ]
dim(top_term_3)

top_term_dot_3 <- ggplot(top_term_3, aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T))) +
  geom_point(aes(size = Count, color = pvalue)) +
  # scale_color_gradientn(colors = c("#BF1E27", "#FEB466","#F9FCCB","#6296C5","#38489D"))+
  scale_color_gradientn(colors = rev(c("navy", "white", "firebrick3"))) +
  guides(color = guide_colorbar(reverse = T)) +
  labs(x = "-log10(Pvalue)", y = NULL, color = "Pvalue") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(lineheight = 1),
    text = element_text(size = 12), # panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10), legend.title = element_text(size = 10)
  )
top_term_dot_3
ggsave(top_term_dot_3, file = "The dotplot of significant GO BP in Group3.pdf", width = 6, height = 6)

# ---------------------- KEGG ----------------------
# kegg
dir.create("./Enrich/KEGG/", showWarnings = F, recursive = T)
eklist <- list()
for (i in names(genelist)) {
  genes <- genelist[[i]]
  # genes2 <- filter(mat, Symbol %in% genes) %>% .$ID
  print(paste0(i, " : ", length(genes)))
  gene_ids <- bitr(genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  ekegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  ekegg2 <- setReadable(ekegg,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID"
  )

  write.csv(ekegg2@result, file = paste0("./Enrich/KEGG/KEGG_", gsub(" ", "_", i), "_results.csv"))
  eklist[[i]] <- ekegg2
}

# signal transduction
kegg_1 <- filter(eklist[[1]]@result, pvalue < 0.05, subcategory == "Signal transduction", Count >= 15) %>% mutate(group = "Group 1")
kegg_3 <- filter(eklist[[3]]@result, pvalue < 0.05, subcategory == "Signal transduction") %>% mutate(group = "Group 3")
kegg_proc <- rbind(kegg_1, kegg_3[!is.na(match(kegg_3$Description, kegg_1$Description)), ])


kegg_proc$Description <- factor(kegg_proc$Description, levels = filter(kegg_proc, group == "Group 1") %>% arrange(desc(pvalue)) %>% .$Description)
# kegg_proc <- kegg_proc %>% filter(Description != "Phospholipase D signaling pathway")

bar_kegg <- ggplot(kegg_proc, aes(x = -log10(pvalue), y = Description, fill = group, group = group)) +
  geom_bar(alpha = 0.8, stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(position = position_dodge(0.7), aes(label = Count), size = 3.5, hjust = -0.3, color = "black") +
  scale_fill_manual(values = c("Group 1" = "#E64B35", "Group 3" = "#00A087"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = NULL) +
  theme_test() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 14, color = "black"),
    # axis.title.x = element_text(size = 12, face = "bold"),
    # axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"
  )
bar_kegg
ggsave(bar_kegg, file = "The barplot of significant signal transduction pathways between Group1 and Group3.pdf", width = 6, height = 6)


# ---------------------- TF ----------------------
tf_all <- read.csv("DatabaseExtract_v_1.01.csv", header = T, row.names = 1)
tf_all <- tf_all %>% filter(Is.TF. == "Yes")
head(tf_all)
dim(tf_all)

tflist <- genelist %>% lapply(function(i) {
  intersect(i, unique(tf_all$HGNC.symbol))
})
lengths(tflist)


max_len <- max(sapply(tflist, length))
tf.mat <- as.data.frame(lapply(tflist, function(x) {
  length(x) <- max_len
  return(x)
}))
colnames(tf.mat) <- gsub("\\.", " ", colnames(tf.mat))
write.csv(tf.mat, "Supplementary TFs of four groups.csv", row.names = TRUE, na = "")


crs <- c(
  "TRPC3", "CRLF1", "FZD3", "ADGRF4", "GPR162", "PTPRQ",
  "NPY1R", "KDR", "PTPN5", "LPAR3", "GABARAPL1", "P2RY1", "NPR3", "TRPC6",
  "ERBB3", "TGFBR1", "ITPR1", "TRPC4", "GPR155", "CNR1", "IL6R", "ROR2",
  "S1PR3", "P2RY14", "PRLR", "RARRES1", "SCARA5", "INSR", "EPHA4", "GPRC5B",
  "FZD4", "RGMA", "IGF1R", "PTCH1"
)

# heatmap: group1 + group3
annotation_colors.tmp <- list(
  Condition = c("NCs" = "blue", "DCs" = "red", "F2tDCs" = "green"),
  Group = c("Group 1" = "#E64B35", "Group 3" = "#00A087")
)
tf_ph <- pheatmap::pheatmap(
  mat = pdata[unique(unlist(tflist[c(1, 3)])), ], # as.matrix(pdata[unique(unlist(tflist_13)), ]),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  legend_labels = FALSE,
  annotation_col = col_group_df,
  annotation_colors = annotation_colors.tmp,
  annotation_row = deg_groups_df,
  annotation_names_row = F, annotation_names_col = F, cellwidth = 20, cellheight = 7, treeheight_row = 0,
  show_rownames = T,
  show_colnames = T,
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
  clustering_callback = callback,
  fontsize = 6, # , fontsize_row = 8, fontsize_col = 8,
  name = NULL, border_color = NA
)
tf_ph
ggsave(tf_ph, file = "The heatmap of tf in group1 and group3.pdf", width = 6, height = 15)


# # ---------------------- membrane protein ----------------------

cr_ph <- pheatmap::pheatmap(
  mat = pdata[crs, ], # as.matrix(pdata[unique(unlist(tflist_13)), ]),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend_labels = FALSE,
  annotation_col = col_group_df,
  annotation_colors = annotation_colors,
  annotation_row = deg_groups_df,
  annotation_names_row = F, annotation_names_col = F, cellwidth = 20, cellheight = 7,
  show_rownames = T,
  show_colnames = T,
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
  clustering_callback = callback,
  fontsize = 6, # , fontsize_row = 8, fontsize_col = 8,
  name = NULL, border_color = NA
)
cr_ph
ggsave(cr_ph, file = "The heatmap of membrane protein genes.pdf", width = 6, height = 6)


# # ---------------------- Inflammatory ----------------------

library(msigdbr)

cytokine_genes <- sapply(c("TNFSF", "IL", "CCL", "CXCL"), function(i) {
  ind <- grepl(paste0("^", i), unique(unlist(genelist)))
  unique(unlist(genelist))[ind]
}) %>%
  unlist() %>%
  unique()
cytokine_genes <- cytokine_genes[-3]
cytokine_genes <- deg_groups_df[cytokine_genes, , drop = F] %>%
  mutate(Cytokine = rownames(.)) %>%
  group_by(Group) %>%
  arrange(Group, .by_group = F) %>%
  .$Cytokine


cytokine_ph <- pheatmap::pheatmap(
  mat = pdata[cytokine_genes, ], # as.matrix(pdata[unique(unlist(tflist_13)), ]),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  cluster_rows = F,
  cluster_cols = FALSE,
  legend_labels = FALSE,
  annotation_col = col_group_df,
  annotation_colors = annotation_colors,
  annotation_row = deg_groups_df,
  annotation_names_row = F, annotation_names_col = F, cellwidth = 20, cellheight = 7,
  show_rownames = T,
  show_colnames = T,
  angle_col = 45,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
  clustering_callback = callback,
  fontsize = 6,
  name = NULL, border_color = NA
)
cytokine_ph
ggsave(cytokine_ph, file = "The heatmap of cytokines.pdf", width = 6, height = 6)
