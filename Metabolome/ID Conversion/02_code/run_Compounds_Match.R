#' Match compound names to various biological databases
#'
#' @param compound_list A vector of compound names to be matched
#' @param id_types Character vector specifying which databases to use for matching (e.g., c("kegg", "hmdb"))
#' @param db_path Path to the directory containing database files (default: "./")
#' @param max_dist Maximum allowed Jaro-Winkler distance for fuzzy matching (default: 0.6)
#' @param string_method String distance method to use (default: "jw" for Jaro-Winkler)
#' @param replace_greek Logical indicating whether to replace Greek letters (default: TRUE)
#' @param output_file Optional path to save results as Excel file (default: NULL)
#' 
#' @return A data frame with matching results where similarity scores range from 0-1 (higher = more similar)
#' @export
compounds_match <- function(compound_list,
                            id_types = c("kegg", "hmdb"),
                            db_path = "./",
                            max_dist = 0.6,
                            string_method = "jw",
                            replace_greek = TRUE,
                            output_file = NULL) {
  
  # Load required libraries
  suppressPackageStartupMessages({
    library(openxlsx)
    library(readr)
    library(dplyr)
    library(RSQLite)
    library(DBI)
    library(stringdist)
  })
  
  # Validate input
  if (!is.character(compound_list)) {
    stop("compound_list must be a character vector")
  }
  
  if (!all(id_types %in% c("kegg", "hmdb"))) {
    stop("Currently only 'kegg' and 'hmdb' are supported as id_types")
  }
  
  # Preprocess compound names
  compound_list <- as.character(compound_list)
  
  if (replace_greek) {
    source(paste0(db_path, "../02_code/replace_greek_letters.R"))
    processed_names <- replace_greek_letters(compound_list)
  }
  processed_names <- tolower(processed_names)
  
  # Initialize result data frame
  result_df <- data.frame(
    original_name = compound_list,
    processed_name = processed_names,
    stringsAsFactors = FALSE
  )
  
  # Connect to SQLite database
  db_file <- paste0(db_path, "metabolite_map.db")
  con <- dbConnect(RSQLite::SQLite(), db_file)
  
  # KEGG matching
  if ("kegg" %in% id_types) {
    # Read KEGG data
    kegg_data <- tryCatch({
      dbGetQuery(con, "SELECT DISTINCT compound_names, kegg_id FROM kegg_map_table")
    }, error = function(e) {
      stop("Failed to read KEGG data from database. Please ensure the KEGG table exists.")
    })
    
    # Initialize KEGG result columns
    result_df$matched_kegg_name <- NA_character_
    result_df$matched_kegg_id <- NA_character_
    result_df$match_type_kegg <- NA_character_
    result_df$similarity_score_kegg <- NA_real_  # Changed from match_score to similarity_score
    
    # Perform matching
    for (i in seq_along(processed_names)) {
      query <- processed_names[i]
      distances <- stringdist::stringdist(query, kegg_data$compound_names, method = string_method)
      min_dist <- min(distances)
      min_index <- which.min(distances)
      
      if (min_dist <= max_dist) {
        result_df$matched_kegg_name[i] <- kegg_data$compound_names[min_index]
        result_df$matched_kegg_id[i] <- kegg_data$kegg_id[min_index]
        # Convert distance to similarity score (1 - distance)
        result_df$similarity_score_kegg[i] <- 1 - min_dist
        result_df$match_type_kegg[i] <- ifelse(min_dist == 0, "exact", "fuzzy")
      }
    }
  }
  
  # HMDB matching
  if ("hmdb" %in% id_types) {
    # Read HMDB data
    hmdb_data <- tryCatch({
      dbGetQuery(con, "SELECT DISTINCT synonyms, accession FROM hmdb_map_table")
    }, error = function(e) {
      stop("Failed to read HMDB data from database. Please ensure the HMDB table exists.")
    })
    
    # Initialize HMDB result columns
    result_df$matched_hmdb_name <- NA_character_
    result_df$matched_hmdb_id <- NA_character_
    result_df$match_type_hmdb <- NA_character_
    result_df$similarity_score_hmdb <- NA_real_  # Changed from match_score to similarity_score
    
    # Perform matching
    for (i in seq_along(processed_names)) {
      query <- processed_names[i]
      distances <- stringdist::stringdist(query, hmdb_data$synonyms, method = string_method)
      min_dist <- min(distances)
      min_index <- which.min(distances)
      
      if (min_dist <= max_dist) {
        result_df$matched_hmdb_name[i] <- hmdb_data$synonyms[min_index]
        result_df$matched_hmdb_id[i] <- hmdb_data$accession[min_index]
        # Convert distance to similarity score (1 - distance)
        result_df$similarity_score_hmdb[i] <- 1 - min_dist
        result_df$match_type_hmdb[i] <- ifelse(min_dist == 0, "exact", "fuzzy")
      }
    }
  }
  
  # Disconnect from database
  dbDisconnect(con)
  
  # Save results if output file specified
  if (!is.null(output_file)) {
    write.xlsx(result_df, output_file)
    message("Results saved to: ", output_file)
  }
  
  return(result_df)
}

# 使用示例
# metabolites <- read.xlsx("./Metabolome/ID Conversion/01_data/DEM_name.xlsx")[,1]
# results <- compounds_match(metabolites, id_types = c("kegg", "hmdb"))
