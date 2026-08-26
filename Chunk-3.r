


# ---------------------- GO ----------------------
genelist <- split(deg_dat$Symbol, deg_dat$Group) %>% lapply(unique)

dir.create("./Enrich/Go/", showWarnings = F, recursive = T)

# # group1 + group3
# genes <- unique(unlist(genelist[c(1, 3)]))
# gene_ids <- bitr(genes,
#                  fromType = "SYMBOL",
#                  toType = "ENTREZID",
#                  OrgDb = org.Hs.eg.db)
# ego <- enrichGO(gene=gene_ids$ENTREZID,
#                 OrgDb=org.Hs.eg.db,
#                 keyType="ENTREZID",
#                 ont="all",
#                 pAdjustMethod="BH",
#                 pvalueCutoff=1,
#                 qvalueCutoff=1,
#                 readable = TRUE)
# 
# write.csv(ego@result, file = paste0("./Enrich/Go/Go_mix_group1_and_group3_results.csv"))
# 
# 
# 
# mf_barplot <- barplot(ego,
#         showCategory = c("growth factor binding",                          
#                          "semaphorin receptor binding",                     
#                          "alcohol binding",                                 
#                          "calmodulin binding",                              
#                          "semaphorin receptor activity",                    
#                          "collagen binding",                                
#                          "heparin binding",                                 
#                          "transforming growth factor beta receptor binding",
#                          "integrin binding",                                
#                          "glycosaminoglycan binding",                       
#                          "SMAD binding",                                    
#                          "calcium ion transmembrane transporter activity",                         
#                          "PDZ domain binding",                                                     
#                          "insulin-like growth factor binding",                                     
#                          "heparan sulfate proteoglycan binding",                                   
#                          "insulin-like growth factor I binding",                                   
#                          "sterol binding",
#                          "metalloendopeptidase inhibitor activity"),
#         x = "Count",
#         color = "pvalue",
#         label_format = 50,
#         font.size = 12,
#         title = "")
# mf_barplot
# ggsave(mf_barplot, file = "./Enrich/Go/The barplot of Go MF pathways in group1 and group3 .pdf", width = 8, height = 6)

egolist <- list()
dir.create("./Enrich/Go", showWarnings = F, recursive = T)
for(i in names(genelist)){
  genes <- genelist[[i]]
  # genes2 <- filter(mat, Symbol %in% genes) %>% .$ID
  print(paste0(i, " : ", length(genes)))
  gene_ids <- bitr(genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  ego<-enrichGO(
    gene=gene_ids$ENTREZID,
    OrgDb=org.Hs.eg.db,
    keyType="ENTREZID",
    ont="all",
    pAdjustMethod="none",
    pvalueCutoff=0.05,
    qvalueCutoff=0.2,
    minGSSize = 10,
    maxGSSize = 50000,
    readable = TRUE
  )

  write.csv(ego@result, file = paste0("./Enrich/Go/Go_", gsub(" ", "_", i), "_results.csv"))
  egolist[[i]] <- ego
}
# for(i in names(egolist)){
#   ego <- egolist[[i]]
#   write.csv(ego@result, file = paste0("./Enrich/Go/Go_", gsub(" ", "_", i), "_results.csv"))
# }

# egolist[[1]] %>% 
library(GO.db)

BPancestor <- as.list(GOBPANCESTOR)
BPancestor[["GO:0007169"]]
BPancestor[["GO:0035556"]]

select(GO.db, keys = head(BPancestor[["GO:0007169"]],-1),
       columns = c("TERM", "ONTOLOGY"))
goids <- keys(GO.db, keytype = "GOID")
goids_filtered <- do.call(rbind, lapply(goids, function(x) {
  if (length(BPancestor[[x]]) >= 11 && length(BPancestor[[x]]) <= 13) {
    return(select(GO.db, keys = x, columns = c("TERM")))
  }
}))

top_term_1 <- filter(egolist[[1]]@result, ONTOLOGY == "BP") %>% top_n(n = 200, wt = Count)
top_term_1 <- top_term_1[!is.na(match(rownames(top_term_1), goids_filtered$GOID)), ]
dim(top_term_1)
# top_term_bar <- ggplot(top_term, aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T)))+
#   geom_bar(stat = "identity", aes(fill = pvalue))+
#   scale_fill_gradientn(colours = rev(c("navy", "white", "firebrick3")))+
#   guides(fill = guide_colorbar(reverse = T))+
#   labs(x = "-log10(Pvalue)", y = NULL, fill = "Pvalue")+
#   theme_bw()+
#   theme(text = element_text(size = 12))
# top_term_bar
# ggsave(top_term_bar, file = "The barplot of significant GO BP in Group1.pdf", width = 7, height = 5)

top_term_dot_1 <- ggplot(top_term_1, aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T)))+
  geom_point(aes(size = Count, color = pvalue))+
  # scale_color_gradientn(colors = c("#BF1E27", "#FEB466","#F9FCCB","#6296C5","#38489D"))+
  scale_color_gradientn(colors = rev(c("navy", "white", "firebrick3")))+
  guides(color = guide_colorbar(reverse = T))+
  labs(x = "-log10(Pvalue)", y = NULL, color = "Pvalue")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))+
  theme_bw()+
  theme(axis.text.y = element_text(lineheight = 0.8),legend.position = "right",
        text = element_text(size = 12),#panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10)
  )
top_term_dot_1 
ggsave(top_term_dot_1, file = "The dotplot of significant GO BP in Group1.pdf", width = 6, height = 6)


top_term_3 <- filter(egolist[[3]]@result, ONTOLOGY == "BP") %>% top_n(n = 150, wt = Count)
top_term_3 <- top_term_3[!is.na(match(rownames(top_term_3), goids_filtered$GOID)), ]
dim(top_term_3)

top_term_dot_3 <- ggplot(top_term_3, aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T)))+
  geom_point(aes(size = Count, color = pvalue))+
  # scale_color_gradientn(colors = c("#BF1E27", "#FEB466","#F9FCCB","#6296C5","#38489D"))+
  scale_color_gradientn(colors = rev(c("navy", "white", "firebrick3")))+
  guides(color = guide_colorbar(reverse = T))+
  labs(x = "-log10(Pvalue)", y = NULL, color = "Pvalue")+
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50))+
  theme_bw()+
  theme(axis.text.y = element_text(lineheight = 1),
        text = element_text(size = 12), #panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),legend.title = element_text(size = 10)
  )
top_term_dot_3 
ggsave(top_term_dot_3, file = "The dotplot of significant GO BP in Group3.pdf", width = 6, height = 6)

# ---------------------- KEGG ----------------------
# kegg
dir.create("./Enrich/KEGG/", showWarnings = F, recursive = T)
eklist <- list()
for(i in names(genelist)){
  genes <- genelist[[i]]
  # genes2 <- filter(mat, Symbol %in% genes) %>% .$ID
  print(paste0(i, " : ", length(genes)))
  gene_ids <- bitr(genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  ekegg <- enrichKEGG(
    gene = gene_ids$ENTREZID,
    organism = 'hsa',
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  ekegg2 <- setReadable(ekegg, 
                       OrgDb = "org.Hs.eg.db", 
                       keyType = "ENTREZID") 
  
  write.csv(ekegg2@result, file = paste0("./Enrich/KEGG/KEGG_", gsub(" ", "_", i), "_results.csv"))
  eklist[[i]] <- ekegg2
}
# for(i in names(eklist)){
#   ekegg <- eklist[[i]]
#   write.csv(ekegg@result, file = paste0("./Enrich/KEGG/KEGG_", gsub(" ", "_", i), "_results.csv"))
# }

# signal transduction
kegg_1 <- filter(eklist[[1]]@result, pvalue < 0.05, subcategory == "Signal transduction", Count >= 15) %>% mutate(group = "Group 1")
kegg_3 <- filter(eklist[[3]]@result, pvalue < 0.05, subcategory == "Signal transduction") %>% mutate(group = "Group 3")
kegg_proc <- rbind(kegg_1, kegg_3[!is.na(match(kegg_3$Description, kegg_1$Description)), ])


kegg_proc$Description <- factor(kegg_proc$Description, levels = filter(kegg_proc, group == "Group 1") %>% arrange(desc(pvalue)) %>% .$Description)
# kegg_proc <- kegg_proc %>% filter(Description != "Phospholipase D signaling pathway")

bar_kegg <- ggplot(kegg_proc, aes(x = -log10(pvalue), y = Description, fill = group, group = group)) +
  geom_bar(alpha = 0.8, stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "white", linewidth = 0.2) +
  geom_text(position = position_dodge(0.7), aes(label = Count), size = 3.5, hjust = -0.3, color = "black")+
  scale_fill_manual(values = c("Group 1" = "#E64B35", "Group 3" = "#00A087"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = NULL) +
  theme_test() +
  theme(text = element_text(size = 12),
    axis.text.y = element_text(size = 14, color = "black"),
    # axis.title.x = element_text(size = 12, face = "bold"),
    # axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"
  )
bar_kegg
ggsave(bar_kegg, file = "The barplot of significant signal transduction pathways between Group1 and Group3.pdf", width = 6, height = 6)

# group1_kegg <- ggplot(filter(kegg_proc, group == "Group 1"), aes(x = -log10(pvalue), y = reorder(Description, pvalue, decreasing = T)))+
#   geom_point(aes(size = Count, color = pvalue))+
#   # scale_color_gradientn(colors = c("#BF1E27", "#FEB466","#F9FCCB","#6296C5","#38489D"))+
#   scale_color_gradientn(colors = rev(c("navy", "white", "firebrick3")))+
#   guides(color = guide_colorbar(reverse = T))+
#   labs(x = "-log10(Pvalue)", y = NULL, color = "Pvalue")+
#   theme_bw()+
#   theme(text = element_text(size = 12)#,panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#         )
# group1_kegg
# ggsave(group1_kegg, file = "The barplot of significant signal transduction pathways in Group1.pdf", width = 6, height = 6)


# ---------------------- GSEA ----------------------

make_geneList <- function(res) {
  res$SYMBOL <- mat[match(rownames(res), mat$ID), "Symbol"]
  gene_ids <- bitr(res$SYMBOL,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
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
geneList_DCs_vs_DCs_FGF2  <- make_geneList(res_DCs_vs_DCs_FGF2)
geneList_NCs_vs_DCs_FGF2 <- make_geneList(res_NCs_vs_DCs_FGF2)

all_glist <-list(geneList_NCs_vs_DCs,geneList_DCs_vs_DCs_FGF2,geneList_NCs_vs_DCs_FGF2)

# loop to enrich
m_gsea_list <- lapply(all_glist, function(x){
  
  gsea_kegg <- clusterProfiler::gseKEGG(
    geneList = x,
    organism = "hsa",         # hsa 代表人类 (Homo sapiens)
    minGSSize = 10,           # 最小基因集大小
    maxGSSize = 500,          # 最大基因集大小
    pvalueCutoff = 1,      # p 值阈值
    verbose = FALSE,        # 不输出详细信息
    eps = 0,
    seed = TRUE
  )
  gsea_kegg
})  
names(m_gsea_list) <- c("NCs_vs_DCs", "DCs_vs_F2tDCs", "NCs_vs_F2tDCs")


geneset.info <- m_gsea_list[[1]] %>% data.frame() %>% filter(Description %in% c("FoxO signaling pathway",
                                                                              "MAPK signaling pathway",
                                                                              "JAK-STAT signaling pathway",
                                                                              "PI3K-Akt signaling pathway",
                                                                              "cAMP signaling pathway")) %>% dplyr::select(Description)
geneset.info
geneset.id <- rownames(geneset.info)
names(geneset.id) <- geneset.info$Description
geneset.id

for(i in 1:length(geneset.id)){
  
  gsea_plot <- GSEAmultiGP_modified(gsea_list = m_gsea_list,
                           geneSetID = geneset.id[[i]],
                           exp_name = gsub("_", " ", names(m_gsea_list)),
                           curve.col = c("#88a0dc","#ed968c","#f9d14a"),
                           addPval = T,pvalX = 0.01,pvalY = 0.01,pDigit = 3,base_size = 12) 
  # gsea_plot
  ggsave(gsea_plot, file = paste0("./GSEA/The GSEA plot of ", names(geneset.id)[[i]], ".pdf"), width = 8, height = 6, create.dir = T)
  print(i)
}