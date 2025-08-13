library(KEGGREST)
library(dplyr)
library(readr)
library(openxlsx)

## 输入所有人类 KEGG 通路
kegg_all <- read.csv("./ClusterProfiler/kegg_category.csv", row.names = 1, stringsAsFactors = FALSE)
head(kegg_all)

## -------- 1) 根据课题定制：黑名单规则 --------
# A. 按 subcategory 批量加入黑名单（强烈建议）
blacklist_subcategory <- c(
  # Human Diseases - Infectious diseases 等
  "Cancer: specific types",
  "Infectious diseases: viral",
  "Infectious diseases: bacterial",
  "Infectious diseases: parasitic",
  "Immune disease",
  "Neurodegenerative disease",
  "Substance dependence",
  "Cardiovascular disease",
  "Endocrine and metabolic disease",
  "Drug resistance: antimicrobial",
  # Organismal Systems - 明显不相干的生理系统
  "Endocrine system",
  "Circulatory system",
  "Digestive system",
  "Sensory system",
  "Nervous system",
  "Digestive system",
  "Excretory system"
  # 注：避免把 "Immune system"一刀切删除
)

# B. 按通路名称关键词（不区分大小写）加入黑名单
#   这里列出常见“误富集”的传染病与病原体相关通路关键字
blacklist_name <- c(
  # 病毒相关
  "covid", "coronavirus", "influenza", "measles", "mumps", "hiv", "htlv", "papillomavirus",
  "hepatitis", "herpes", "epstein-barr", "adenovirus",
  # 细菌/寄生虫相关
  "tuberculosis", "cholera", "shigellosis", "salmonella", "pertussis", "legionellosis",
  "leishman", "trypanosom", "malaria", "toxoplasmosis", "amebiasis",
  # 其他
  "staphylococcus", "vibrio", "helicobacter", "clostridium", "borrelia", "treponema"
)

# C. 按 ID（纯数字，等价于 hsa04110→写 4110 或 04110）加入黑名单（可留空）
blacklist_ids <- c(
  # 把几个已知常误富集、与你课题不相干的通路直接列入
)

## -------- 2) 白名单兜底（避免误删关键机制通路）--------
# 与 Venetoclax 耐药/AML 强相关的通路名关键词（不区分大小写）
whitelist_name <- c(
  "apoptosis", "bcl", "p53", "mitophagy", "autophagy", "lysosome",
  "oxidative phosphorylation", "tca cycle", "fatty acid degradation", "fatty acid metabolism",
  "glutathione metabolism", "peroxisome", "foxo", "mTOR", "PI3K-Akt", "MAPK", "NF-kappa B",
  "HIF-1 signaling", "amino acid metabolism", "glycine", "serine", "threonine", "one carbon",
  "cysteine", "methionine", "arginine", "proline", "glutamine",
  "hematopoietic cell lineage", "leukemia"
)

# 如你手头掌握的关键通路 ID（可选，示例）
whitelist_ids <- c(
  4210,   # Apoptosis (hsa04210)
  4151,   # PI3K-Akt signaling
  4115,   # p53 signaling
  190,    # Oxidative phosphorylation (hsa00190)
  5410,   # HIF-1 signaling
  4330,   # Jak-STAT (常见于血液肿瘤信号)
  5200,   # Pathways in cancer
  5213,   # Endocrine resistance（与耐药可能相关，慎删）
  5202    # Transcriptional misregulation in cancer
  # 按需增减
)

## -------- 3) 生成黑名单（多条件并集），再用白名单兜底 --------
# 标记黑名单候选
cand_black <- with(kegg_all,
                   (subcategory %in% blacklist_subcategory) |
                     grepl(paste(blacklist_name, collapse = "|"), name, ignore.case = TRUE) |
                     (id %in% blacklist_ids)
)

blacklist <- kegg_all[cand_black, , drop = FALSE]

# 白名单保护：凡是命中白名单关键词或 ID 的，从黑名单里剔除
keep_mask <- with(blacklist,
                  grepl(paste(whitelist_name, collapse = "|"), name, ignore.case = TRUE) |
                    (id %in% whitelist_ids)
)

blacklist <- blacklist[!keep_mask, , drop = FALSE]

## -------- 4) 保存黑名单 --------
write.csv(blacklist, "kegg_blacklist.csv", row.names = FALSE)
cat("✅ 黑名单已生成，共", nrow(blacklist), "条通路。\n")

## （可选）预览每类被屏蔽数量，便于你校对是否过严
table_black_by_subcat <- sort(table(blacklist$subcategory), decreasing = TRUE)
print(table_black_by_subcat)
