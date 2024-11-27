library(data.table)
library(tidyverse)

# 设置文件夹路径
result_dir <- "../07-separate_cell_type_DEG/results/"

# 获取所有CSV文件的文件名
file_list <- list.files(result_dir, pattern = "\\.csv$", full.names = TRUE) %>%
  .[!str_detect(., "combined")]


# 初始化一个空的数据框用于存储结果
OLAH <- data.frame()

# 遍历每个文件
for (file in file_list) {
  # 读取数据
  data <- fread(file)
  
  # 提取文件名（去掉路径和后缀）
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # 添加新列记录文件名信息
  data <- data %>%
    mutate(cell_cytokine = file_name)
  
  # 筛选 gene == "OLAH" 的行
  data_filtered <- data %>%
    filter(gene == "OLAH")
  
  # 合并结果
  OLAH <- bind_rows(OLAH, data_filtered)
}


cDC <- fread("../06-test_in_Jupyter/results/cDC_combined_OLAH_rows.csv")
CD8_Naive <- fread("../06-test_in_Jupyter/results/CD8_Naive_combined_OLAH_rows.csv")
OLAH <- OLAH[,-1]
OLAH <- bind_rows(OLAH, cDC)
OLAH <- bind_rows(OLAH, CD8_Naive)

NROW(unique(OLAH$gene))   # 1
NROW(unique(OLAH$cell_cytokine))   # 1620

# 使用正则表达式基于最后一个下划线分割
OLAH <- OLAH %>%
  mutate(
    cell_type = str_extract(cell_cytokine, ".*(?=_[^_]+$)"),  # 提取最后一个下划线之前的部分
    cytokine = str_extract(cell_cytokine, "[^_]+$")          # 提取最后一个下划线之后的部分
  )

NROW(unique(OLAH$cell_type))   # 18
NROW(unique(OLAH$cytokine))   # 90
unique(OLAH$cell_type)

# 如果需要保存结果为一个新的CSV文件
fwrite(OLAH, "OLAH_combined.csv")

colnames(OLAH)

### Plot
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

# 假设 OLAH 数据框中已有以下列：
# - logfoldchanges：横坐标数据
# - pvals_adj：调整后的 p 值（未取对数）
# - cell_type：颜色分组
# - cytokine：形状分组

# 计算 -log10 转换的 p 值
OLAH <- OLAH %>%
  mutate(
    neg_log_pvals = -log10(pvals_adj),
    point_color = case_when(
      cell_cytokine == "cDC_TSLP" ~ "red",  # 特殊点1
      cell_cytokine == "cDC_IL-4" ~ "blue", # 特殊点2
      TRUE ~ "black"                        # 其他点
    ),
    label = case_when(
      cell_cytokine == "cDC_TSLP" ~ "cDC_TSLP",
      cell_cytokine == "cDC_IL-4" ~ "cDC_IL-4",
      TRUE ~ NA_character_
    )
  )


# 绘制火山图
volcano_plot <- ggplot(OLAH, aes(x = logfoldchanges, y = neg_log_pvals)) +
  geom_point(aes(color = point_color), shape = 16, size = 3, alpha = 0.8) + # 散点图，颜色动态
  scale_color_identity() + # 直接使用指定的颜色
  geom_text_repel(aes(label = label), size = 4, color = "black", max.overlaps = 10) + # 标注特定点
  theme_pubr() +
  labs(
    x = "Log Fold Change",
    y = "-log10(Adjusted P-Value)",
    title = "OLAH in Cytokine vs. PBS"
  ) +
  theme(
    legend.position = "none", # 去掉图例
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

volcano_plot
ggsave("OLAH-volcano_plot.png", plot = volcano_plot, dpi = 300, units = "in", height = 4, width = 4)
ggsave("OLAH-volcano_plot.pdf", plot = volcano_plot, dpi = 300, units = "in", height = 4, width = 4)



### Extract data
library(data.table)
library(tidyverse)

OLAH <- fread("OLAH_combined.csv")

cDC <- OLAH %>% 
  filter(cell_type=="cDC") %>% 
  arrange(pvals, desc(logfoldchanges))

CD14M <- OLAH %>% 
  filter(cell_type=="CD14_Mono") %>% 
  arrange(pvals, desc(logfoldchanges))

CD16M <- OLAH %>% 
  filter(cell_type=="CD16_Mono") %>% 
  arrange(pvals, desc(logfoldchanges))

TNFa <- OLAH %>% 
  filter(cytokine=="TNF-alpha") %>% 
  arrange(pvals, desc(logfoldchanges))



### DEGs in cDC with TSLP vs. PBS
DEG <- fread("/home/jqu/project/parse_biosiences/06-test_in_Jupyter/results/cDC_TSLP.csv")

# Volcano
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

DEG <- DEG %>%
  mutate(
    neg_log_pvals = -log10(pvals_adj),
    point_color = case_when(
      logfoldchanges >= 1 & pvals_adj < 0.01 ~ "red",  # Up
      logfoldchanges <= -1 & pvals_adj < 0.01 ~ "blue", # Down
      TRUE ~ "black"                        # NS
    ),
    regulation = case_when(
      logfoldchanges >= 1 & pvals_adj < 0.01 ~ "Up",
      logfoldchanges <= -1 & pvals_adj < 0.01 ~ "Down",
      TRUE ~ "NS"
    ),
    label = case_when(
      gene == "OLAH" ~ "OLAH",
      TRUE ~ NA_character_
    )
  )

fwrite(DEG, "DEG-cDC_TSLP.csv")

NROW(DEG[DEG$regulation=="Up",])   # 392
NROW(DEG[DEG$regulation=="Down",])   # 199


volcano_plot2 <- ggplot(DEG, aes(x = logfoldchanges, y = neg_log_pvals)) +
  geom_point(aes(color = point_color), shape = 16, size = 1.5, alpha = 0.5) +
  scale_color_identity() +
  geom_text_repel(aes(label = label), size = 3, color = "black", max.overlaps = 10) +
  theme_pubr() +
  labs(
    x = "Log Fold Change",
    y = "-log10(Adjusted P-Value)",
    title = "DEGs in cDC with TSLP vs. PBS"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
volcano_plot2
ggsave("volcano-cDC_TSLP.png", plot = volcano_plot2, dpi = 300, units = "in", height = 4, width = 4)
ggsave("volcano-cDC_TSLP.pdf", plot = volcano_plot2, units = "in", height = 4, width = 4)


### clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)



NROW(unique(DEG$gene))   # 40352

gene_ID <- bitr(DEG$gene,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
NROW(unique(gene_ID$SYMBOL))   # 27345
NROW(unique(gene_ID$ENTREZID))   # 27348
gene_ID$SYMBOL[duplicated(gene_ID$SYMBOL)]   # "TEC"  "HBD"  "MMD2"

gene_ID[gene_ID$SYMBOL=="TEC",]
# SYMBOL  ENTREZID
# 1450    TEC      7006
# 1451    TEC 100124696

gene_ID[gene_ID$SYMBOL=="HBD",]
# SYMBOL  ENTREZID
# 22590    HBD      3045
# 22591    HBD 100187828

gene_ID[gene_ID$SYMBOL=="MMD2",]
# SYMBOL  ENTREZID
# 33736   MMD2    221938
# 33737   MMD2 100505381

# Remove 
# 1451    TEC 100124696
# 22591    HBD 100187828
# 33737   MMD2 100505381
gene_ID <- gene_ID[,-c(1451,22591,33737)]

DEG <- DEG %>% 
  left_join(gene_ID, by=c("gene"="SYMBOL"))
NROW(unique(DEG$gene))   # 40352

saveRDS(DEG, "DEG-cDC_TSLP-ENTREZID.rds")
fwrite(DEG, "DEG-cDC_TSLP-ENTREZID.csv")

DEG <- DEG %>% 
  filter(regulation != "NS")

mydir <- "cDC_TSLP_enrich/"
if (!file.exists(mydir)){
  dir.create(mydir)
}

# KEGG
KEGG <- enrichKEGG(gene = DEG$ENTREZID,
                   organism ='hsa',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   use_internal_data =FALSE)
KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
fwrite(as.data.frame(KEGG@result),paste0(mydir,"KEGG.csv"))
saveRDS(KEGG,paste0(mydir,"KEGG.rds"))
dot_KEGG <- dotplot(KEGG)
saveRDS(dot_KEGG, paste0(mydir,"KEGG.rds"))
ggsave(filename = paste0(mydir,"KEGG.pdf"), plot = dot_KEGG, width = 6, height = 6, units = "in")  
ggsave(filename = paste0(mydir,"KEGG.png"), plot = dot_KEGG, width = 6, height = 6, units = "in", dpi = 300)  

# GO
GO <- enrichGO(gene = DEG$ENTREZID,
               OrgDb = "org.Hs.eg.db",
               keyType = "ENTREZID",
               ont = "BP",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05,
               readable=TRUE)
GO <- simplify(GO,cutoff=0.7,
               by="p.adjust",
               select_fun=min,
               measure = "Wang")
fwrite(as.data.frame(GO@result),paste0(mydir,"GO_BP.csv"))
saveRDS(GO,paste0(mydir,"GO_BP.rds"))
dot_GO <- dotplot(GO)
saveRDS(dot_GO, paste0(mydir,"GO_BP.rds"))
ggsave(filename = paste0(mydir,"GO_BP.pdf"), plot = dot_GO, width = 6, height = 6, units = "in")  
ggsave(filename = paste0(mydir,"GO_BP.png"), plot = dot_GO, width = 6, height = 6, units = "in", dpi = 300)  
