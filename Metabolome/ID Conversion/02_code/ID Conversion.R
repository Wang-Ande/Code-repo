# HDMB database ----
library(readr)
hmdb <- read_tsv("./Metabolome/ID Conversion/01_data/hmdb_metabolites_parsed.tsv")
head(hmdb)

# 将 name 合并进 synonyms
hmdb <- hmdb %>%
  mutate(synonyms = ifelse(is.na(synonyms) | synonyms == "", 
                           name, 
                           paste(name, synonyms, sep = "; "))) %>%
  select(-name)

# 展开 synonyms
library(tidyr)
hmdb_long <- hmdb |> 
  separate_rows(synonyms, sep = ";") |> 
  dplyr::mutate(synonyms = trimws(synonyms))  # 去掉首尾空格
head(hmdb_long)

# 检查是否有希腊字母
greek_letters <- c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", 
                   "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", 
                   "ρ", "σ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω")
contains_greek <- sapply(hmdb_long$synonyms, 
                         function(x) any(grepl(paste(greek_letters, collapse = "|"), x)))
greek_name <- hmdb_long[contains_greek, ]

# 将希腊字母转换为英语字母
source("./Metabolome/ID Conversion/02_code/replace_greek_letters.R")
hmdb_long$synonyms <- replace_greek_letters(hmdb_long$synonyms)

# 大写字母转换为小写 
library(dplyr)
hmdb_tolower <- hmdb_long %>%
  mutate(
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
# 展开compound_names
kegg_long <- kegg |> 
  separate_rows(compound_names, sep = ";") |> 
  dplyr::mutate(compound_names = trimws(compound_names))  # 去掉首尾空格
head(kegg_long)

# 检查是否有希腊文字
greek_letters <- c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", 
                   "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", 
                   "ρ", "σ", "ς", "τ", "υ", "φ", "χ", "ψ", "ω")
contains_greek <- sapply(kegg_long$compound_names, 
                         function(x) any(grepl(paste(greek_letters, collapse = "|"), x)))
kegg_long[contains_greek, ]
# 将希腊字母转换为英语字母
source("./Metabolome/ID Conversion/02_code/replace_greek_letters.R")
kegg_long$compound_names <- replace_greek_letters(kegg_long$compound_names)

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
metabolites <- metabolites[,-2]
head(metabolites)

# 大写换为小写
metabolites_tolower <- tolower(metabolites)
# 希腊字母换为英文
source("./Metabolome/ID Conversion/02_code/replace_greek_letters.R")
metabolites_tolower <- replace_greek_letters(metabolites_tolower)

# 写入 KEGG 数据库 ----
con <- dbConnect(RSQLite::SQLite(), "./Metabolome/ID Conversion/03_result/kegg_compound.db")
dbWriteTable(con, "kegg_map_table", kegg_database, overwrite = TRUE)

# 从数据库读取 compound_names 和 kegg_id
kegg_names <- dbGetQuery(con, "SELECT DISTINCT compound_names, kegg_id FROM kegg_map_table")

# 初始化结果向量
matched_kegg_ids_fuzzy <- vector("character", length(metabolites))
matched_kegg_names <- vector("character", length(metabolites))
match_type_kegg <- vector("character", length(metabolites))
match_score_kegg <- numeric(length(metabolites))

# 设置匹配容差限制
# jw距离范围[0,1]
max_dist <- 0.6

for (i in seq_along(metabolites_tolower)) {
  query <- metabolites_tolower[i]
  distances <- stringdist::stringdist(query, kegg_names$compound_names, method = "jw")
  min_dist <- min(distances)
  min_index <- which.min(distances)
  match_score_kegg[i] <- ifelse(min_dist <= max_dist, min_dist, NA)
  
  if (min_dist <= max_dist) {
    matched_kegg_ids_fuzzy[i] <- kegg_names$kegg_id[min_index]
    matched_kegg_names[i] <- kegg_names$compound_names[min_index]
    
    if (min_dist == 0) {
      match_type_kegg[i] <- "exact"
    } else {
      match_type_kegg[i] <- "fuzzy"
    }
  } else {
    matched_kegg_ids_fuzzy[i] <- NA
    matched_kegg_names[i] <- NA
    match_type_kegg[i] <- NA
  }
}

# 写入 HMDB 数据库 ----
if (!dbIsValid(con)) {
  con <- dbConnect(RSQLite::SQLite(), "./Metabolome/ID Conversion/03_result/metabolite_map.db")
}

dbWriteTable(con, "hmdb_map_table", hmdb_database, overwrite = TRUE)

# 从数据库读取 compound_names 和 hmdb_id
hmdb_names <- dbGetQuery(con, "SELECT DISTINCT synonyms, accession FROM hmdb_map_table")

# 初始化 HMDB 匹配结果向量
matched_hmdb_ids_fuzzy <- vector("character", length(metabolites))
matched_hmdb_names <- vector("character", length(metabolites))
match_type_hmdb <- vector("character", length(metabolites))
match_score_hmdb <- numeric(length(metabolites))

# 设置匹配容差限制
# jw距离范围[0,1]
max_dist <- 0.6
# 逐个匹配
for (i in seq_along(metabolites_tolower)) {
  query <- metabolites_tolower[i]
  distances <- stringdist::stringdist(query, hmdb_names$synonyms, method = "jw")
  min_dist <- min(distances)
  min_index <- which.min(distances)
  match_score_hmdb[i] <- ifelse(min_dist <= max_dist, min_dist, NA)
  
  if (min_dist <= max_dist) {
    matched_hmdb_ids_fuzzy[i] <- hmdb_names$accession[min_index]
    matched_hmdb_names[i] <- hmdb_names$synonyms[min_index]
    
    if (min_dist == 0) {
      match_type_hmdb[i] <- "exact"
    } else {
      match_type_hmdb[i] <- "fuzzy"
    }
  } else {
    matched_hmdb_ids_fuzzy[i] <- NA
    matched_hmdb_names[i] <- NA
    match_type_hmdb[i] <- NA
  }
}

# 整合 KEGG + HMDB 匹配结果
conversion_result <- data.frame(
  original_name = metabolites,
  processed_name = metabolites_tolower,
  
  matched_kegg_name = matched_kegg_names,
  matched_kegg_id = matched_kegg_ids_fuzzy,
  match_type_kegg = match_type_kegg,
  match_score_kegg = match_score_kegg,
  
  matched_hmdb_name = matched_hmdb_names,
  matched_hmdb_id = matched_hmdb_ids_fuzzy,
  match_type_hmdb = match_type_hmdb,
  match_score_hmdb = match_score_hmdb,
  
  stringsAsFactors = FALSE
)

# 保存结果
write.xlsx(conversion_result, "./Metabolome/ID Conversion/03_result/conversion_result.xlsx")



