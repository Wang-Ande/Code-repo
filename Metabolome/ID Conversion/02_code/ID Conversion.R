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
metabolites <-
