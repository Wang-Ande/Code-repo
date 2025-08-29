#' 更新化合物数据库CSV文件（一对多结构）
#'
#' @param csv_path CSV文件路径
#' @param compound_id 化合物ID
#' @param new_name 要添加的新名称
#' @param check_duplicates 是否检查重复名称（默认TRUE）
#' 
#' @return 逻辑值，表示是否成功更新
#' @export
update_compound_db <- function(csv_path, compound_id, new_name, check_duplicates = TRUE) {
  
  # 读取数据库
  if (!file.exists(csv_path)) {
    stop("CSV文件不存在: ", csv_path)
  }
  db <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # 检查列名是否符合预期
  if (!all(c("id", "name") %in% colnames(db))) {
    stop("CSV文件必须包含'id'和'name'两列")
  }
  
  # 检查是否已存在完全相同的记录
  if (check_duplicates) {
    duplicate <- any(db$id == compound_id & db$name == new_name)
    if (duplicate) {
      message("记录已存在 - ID: ", compound_id, ", 名称: ", new_name)
      return(FALSE)
    }
  }
  
  # 创建新记录
  new_entry <- data.frame(id = compound_id, name = new_name, stringsAsFactors = FALSE)
  
  # 合并并保持原有排序（按ID排序）
  db_updated <- rbind(db, new_entry)
  db_updated <- db_updated[order(db_updated$id), ]
  
  # 写回文件
  tryCatch({
    write.csv(db_updated, csv_path, row.names = FALSE)
    message("成功添加新记录 - ID: ", compound_id, ", 名称: ", new_name)
    return(TRUE)
  }, error = function(e) {
    message("写入文件失败: ", e$message)
    return(FALSE)
  })
}

# 使用示例 ---------------------------------------------------

# 添加新名称到现有ID
update_compound_db("hmdb_database.csv", "HMDB00001", "β-D-葡萄糖")

# 添加全新ID和名称
update_compound_db("hmdb_database.csv", "HMDB99999", "测试化合物")

# 跳过重复检查（强制添加）
update_compound_db("hmdb_database.csv", "HMDB00001", "D-葡萄糖", check_duplicates = FALSE)
