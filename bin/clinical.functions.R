

`%ni%` <- Negate(`%in%`)

make_clean_names <- function(all.clinic){
  clean_names <- iconv(names(all.clinic), from = "", to = "UTF-8", sub = "")
  spaced_names <- gsub(" ", "_", clean_names)
  
  alpha_names <- gsub("[^A-Za-z0-9_]", "", spaced_names)
  names(all.clinic) <- alpha_names
  names(all.clinic) <- make.unique(names(all.clinic), sep = "_")
  return(all.clinic)
}

align_df2_to_df1 <- function(df1, df2) {
  # Get the column names from df1
  target_cols <- names(df1)
  
  # Create missing columns in df2 with NA
  for (col in setdiff(target_cols, names(df2))) {
    df2[[col]] <- NA
  }
  
  # Drop extra columns and reorder to match df1
  df2_aligned <- df2[, target_cols, drop = FALSE]
  
  combined <- rbind(df1, df2_aligned)
  return(combined)
}

fix_excel_dates <- function(df) {
  df_fixed <- df
  
  for (col_name in names(df_fixed)) {
    col <- df_fixed[[col_name]]
    
    # Check if column is character and all non-NA values are numeric strings
    if (is.character(col)) {
      non_na_values <- col[!is.na(col)]
      
      if (all(grepl("^\\d+$", non_na_values))) {
        col_numeric <- as.numeric(col)
        
        df_fixed[[col_name]] <- as.Date(col_numeric, origin = "1899-12-30")
        
      } else {
        df_fixed[[col_name]] <- as.Date(col)
      }
    }
  }
  
  return(df_fixed)
}
