# 使用示例
source("compounds_match.R")
# 示例：仅匹配 KEGG
result_kegg <- compounds_match(metabolites, id_types = "kegg")

# 示例：匹配 KEGG + HMDB（默认）
result_all <- compounds_match(metabolites)

# 添加新名称到现有ID
source("update_compound_db.R")
update_compound_db("hmdb_database.csv", "HMDB00001", "β-D-葡萄糖")


# 自定义路径
result_custom <- compounds_match(
  metabolites,
  kegg_db_path = "your/custom/kegg.db",
  hmdb_db_path = "your/custom/hmdb.db"
)
