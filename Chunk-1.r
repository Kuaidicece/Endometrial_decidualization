


# ---------------------- read matrix ----------------------
mat <- read.table("/Users/niexiner/Documents/Phd/STTT/data/BMK_4_geneExpression/BMK_3_Expression_Statistics/All_gene_counts.list", header = TRUE, sep = "\t", comment.char = "",
    check.names = F)
colnames(mat)[1] <- "ID"
head(mat)

counts <- mat
rownames(counts) <- counts$ID
counts <- counts[, -c(1, 2)]
head(counts)
dim(counts)

# filter genes rowSums(counts) %>% head
counts <- counts[rowSums(counts) != 0, ]
dim(counts)

# ---------------------- DESeq2 ----------------------
# sample information
sample_info <- data.frame(row.names = c("C4003C-1", "C4003C-2", "C4003C-3", "C4003M-1", "C4003M-2", "C4003M-3", "C4003MB-1", "C4003MB-2", "C4003MB-3"), 
                          condition = factor(c("NCs", "NCs", "NCs", "DCs", "DCs", "DCs", "F2tDCs", "F2tDCs", "F2tDCs")))
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
        mutate(ID = row.names(.), tendency = case_when(padj < 0.01 & log2FoldChange > 1 ~ "Up", 
                                                       padj < 0.01 & log2FoldChange < -1 ~ "Down", 
                                                       TRUE ~ "ns"))
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



# tt <- read.table('/Users/niexiner/Documents/Phd/STTT/data/BMK_5_DEG_Analysis/BMK_1_All_DEG/All.DEG_final.xls', header = F) 
# tt2 <-deg %>% lapply(function(i){ filter(i, tendency != 'ns') %>% .$ID }) %>% unlist %>% unique 
# intersect(tt2, tt$V1) %>%length #3778个重叠, 0.9737113.


# visualization Venn
gene_list <- lapply(deg, function(i) {
    filter(i, tendency != "ns") %>%.$ID
})
venn_obj <- venn.diagram(
  x = gene_list,
  filename =NULL,
  fill =c("#88a0dc","#ed968c","#f9d14a"),
  alpha =c(0.5, 0.5, 0.5),
  col =c("#88a0dc","#ed968c","#f9d14a"),
  category.names = c("NCs vs DCs", "DCs vs F2tDCs", "NCs vs F2tDCs"),
  cex =1.2, cat.cex =1.2, cat.fontface ="bold", cat.dist = c(0.1, 0.1, 0.05), margin =0.1)

pdf("Venn_of_DEGs.pdf", width = 6, height = 6)
grid.newpage()
grid.draw(venn_obj)
dev.off()

# table s5
library(tidyr)
deg_genes <- unique(unlist(gene_list))

# read fpkm
fpkm <- read.table("/Users/niexiner/Documents/Phd/STTT/data/BMK_4_geneExpression/BMK_3_Expression_Statistics/All_gene_fpkm.list", 
                   header = TRUE, sep = "\t", comment.char = "",check.names = F)
colnames(fpkm)[1] <- "ID"
head(fpkm)

dat <- filter(fpkm, ID %in% deg_genes) %>% select(-ID)
colnames(dat) <- c("Symbol name", "NCs-1", "NCs-2", "NCs-3", "DCs-1", "DCs-2", "DCs-3", "F2tDCs-1", "F2tDCs-2", "F2tDCs-3")
head(dat)

# save
write.xlsx(list("all deg" = dat), file = "The fpkm of DEGs among NCs, DCs and F2tDCs.xlsx", asTable = FALSE, rowNames = FALSE, keepNA = TRUE,
           firstRow = FALSE)