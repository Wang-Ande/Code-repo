# HDMB database ----
library(readr)
hmdb <- read_tsv("./Metabolome/ID Conversion/01_data/hmdb_metabolites_parsed.tsv")
head(hmdb)

library(tidyr)
hmdb_long <- hmdb |> 
  separate_rows(synonyms, sep = ";") |> 
  dplyr::mutate(synonyms = trimws(synonyms))  # 去掉首尾空格

head(hmdb_long)
write.csv(hmdb_long, file = "./Metabolome/ID Conversion/03_result/hmdb_long.csv")

# 检查是否有希腊字母
greek_letters <- c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", 
                   "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", 
                   "ρ", "σ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω")
contains_greek <- sapply(hmdb_long$name, function(x) any(grepl(paste(greek_letters, collapse = "|"), x)))
greek_name <- hmdb_long[contains_greek, ]
sum(contains_greek == TRUE)
# 将希腊字母转换为英语字母
source("./Metabolome/ID Conversion/02_code/replace_greek_letters.R")
hmdb_long$synonyms <- replace_greek_letters(hmdb_long$synonyms)

# 大写字母转换为小写 
library(dplyr)

hmdb_tolower <- hmdb_long %>%
  mutate(
    name = tolower(name),
    synonyms = tolower(synonyms)
  )
write.csv(hmdb_tolower, file = "./Metabolome/ID Conversion/03_result/hmdb_tolower.csv")

# KEGG database ----
library(readr)
library(tidyr)
library(dplyr)

# 读入数据
kegg <- read_tsv("./Metabolome/ID Conversion/01_data/kegg_compound_list.txt",
                 col_names = c("kegg_id", "compound_names"))

# 拆分 compound_names 成多行
kegg_long <- kegg %>%
  separate_rows(compound_names, sep = ";\\s")  # 分隔符是 ;＋空格

# 去除可能的空格
kegg_long <- kegg_long %>%
  mutate(compound_names = trimws(compound_names))

# 查看结果
print(kegg_long)
write.csv(kegg_long, file = "./Metabolome/ID Conversion/03_result/kegg_long.csv")

# 检查是否有希腊文字
greek_letters <- c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", 
                   "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", 
                   "ρ", "σ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω")
contains_greek <- sapply(kegg_tolower$compound_names, function(x) any(grepl(paste(greek_letters, collapse = "|"), x)))
kegg_tolower[contains_greek, ]
sum(contains_greek == TRUE)

# 大写字母转换为小写 
library(dplyr)

kegg_tolower <- kegg_long %>%
  mutate(
    compound_names = tolower(compound_names)
  )
write.csv(kegg_tolower, file ="./Metabolome/ID Conversion/03_result/kegg_tolower.csv")

# Metabolites list input ----
library(openxlsx)
library(readr)
library(dplyr)
library(RSQLite)
library(DBI)
library(stringdist) # 模糊匹配

kegg_database <- read.csv("./Metabolome/ID Conversion/03_result/kegg_tolower.csv", row.names = 1)
hmdb_database <- read.csv("./Metabolome/ID Conversion/03_result/hmdb_tolower.csv", row.names = 1) 
metabolites <- read.xlsx("./Metabolome/ID Conversion/01_data/DEM_name.xlsx") # 测试数据来自MOLM13 High-Ctrl

kegg_database <- kegg_database[,c("compound_names","kegg_id")]
head(kegg_database)
metabolites <- metabolites[,-2]
head(metabolites)
# 大写换为小写
metabolites_tolower <- tolower(metabolites)
# 希腊字母换为英文
source("./Metabolome/ID Conversion/02_code/replace_greek_letters.R")
metabolites_tolower <- replace_greek_letters(metabolites_tolower)

# 写入 SQLite 数据库 ----
con <- dbConnect(RSQLite::SQLite(), "./Metabolome/ID Conversion/03_result/kegg_compound.db")
dbWriteTable(con, "kegg_map_table", kegg_database, overwrite = TRUE)

# 从数据库读取 compound_names 和 kegg_id
kegg_names <- dbGetQuery(con, "SELECT DISTINCT compound_names, kegg_id FROM kegg_map_table")

# 初始化结果向量
matched_kegg_ids_fuzzy <- vector("character", length(metabolites))
matched_names <- vector("character", length(metabolites))

# 设置匹配容差限制
max_dist <- 2

for (i in seq_along(metabolites_tolower)) {
  query <- metabolites_tolower[i]
  distances <- stringdist::stringdist(query, kegg_names$compound_names, method = "jw")
  min_dist <- min(distances)
  min_index <- which.min(distances)
  
  # 只接受距离小于等于 max_dist 的匹配，否则为 NA
  if (min_dist <= max_dist) {
    matched_kegg_ids_fuzzy[i] <- kegg_names$kegg_id[min_index]
    matched_names[i] <- kegg_names$compound_names[min_index]
  } else {
    matched_kegg_ids_fuzzy[i] <- NA
    matched_names[i] <- NA
  }
}

# 整合结果
conversion_result_fuzzy <- data.frame(
  original_name = metabolites,
  processed_name = metabolites_tolower,
  matched_kegg_name = matched_names,
  matched_kegg_id = matched_kegg_ids_fuzzy,
  stringsAsFactors = FALSE
)

# 查看未匹配的部分
unmatched_fuzzy <- subset(conversion_result_fuzzy, is.na(matched_kegg_id))
