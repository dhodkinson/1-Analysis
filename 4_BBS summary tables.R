#updates
# need to change non-native species to black text

rm(list=ls()) # clear all previous

# ============================================================================
# SPECIES x DATE / FIELD / SITE SUMMARY TABLE GENERATOR
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(flextable)
library(officer)
library(lubridate)
library(sf)

# ----------------------------------------------------------------------------
# 0. CONFIG - EDIT THESE
# ----------------------------------------------------------------------------

CONFIG <- list(
  input_path              = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg',
  layer_name              = '2026_Counts_final',
  bocc5_lookup_path       = "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv",
  taxonomic_lookup_path   = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BBS_species taxonomic order lookup.csv',
  output_dir              = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results/Summer',
  bocc5_order             = c("Red", "Amber", "Green"),
  
  # Set specific dates (e.g., c("14/04/2026", "30/05/2026")) to filter/split by date.
  # SET TO NULL to aggregate all dates into single columns (no date filters).
  dates_to_include        = c("2026-04-14","2026-05-30"),
  
  # Filter breakdown table to a specific group: "Woodland", "Farmland", or NULL (for all)
  filter_field_group      = NULL, 
  
  # Column split options: "site" or "none"
  split_columns_by        = "site",
  
  # Code lists used to categorize Field_Code into Field_Group for filtering
  woodland_codes          = as.character(1:50), 
  farmland_codes          = c(LETTERS, paste0("F", 1:50)),            
  
  include_breeding_status = TRUE,
  sort_by                 = "taxonomic",
  template_path           = NULL,
  
  bocc_text_colors        = c(
    "Red"   = "#C00000",
    "Amber" = "#D97706",
    "Green" = "#15803D"
  )
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1. HELPER FUNCTIONS & DATA CLEANING
# ----------------------------------------------------------------------------

apply_date_filter <- function(df, dates_to_include) {
  if (!"Date_formatted" %in% names(df) || is.null(dates_to_include) || length(dates_to_include) == 0) {
    return(df)
  }
  
  target_dates_parsed <- parse_date_time(dates_to_include, orders = c("dmy", "ymd", "mdy"))
  target_dates_formatted <- format(target_dates_parsed, "%d/%m/%Y")
  target_dates_formatted <- target_dates_formatted[!is.na(target_dates_formatted)]
  
  if (length(target_dates_formatted) > 0) {
    df <- df %>% filter(Date_formatted %in% target_dates_formatted)
  }
  return(df)
}

assign_field_group <- function(df, woodland_codes, farmland_codes) {
  df %>%
    mutate(
      Field_Code_clean = trimws(as.character(Field_Code)),
      Field_Group = case_when(
        Field_Code_clean %in% as.character(woodland_codes) | 
          str_detect(Field_Code_clean, "^[0-9]+$") | 
          str_detect(tolower(Field_Code_clean), "^w") ~ "Woodland",
        
        toupper(Field_Code_clean) %in% toupper(as.character(farmland_codes)) | 
          str_detect(Field_Code_clean, "^[A-Za-z]$") | 
          str_detect(tolower(Field_Code_clean), "^f") ~ "Farmland",
        
        TRUE ~ Field_Code_clean
      )
    ) %>%
    select(-Field_Code_clean)
}

# Read spatial layer from GeoPackage and drop geometry for tabular analysis
dt_spatial <- st_read(dsn = CONFIG$input_path, layer = CONFIG$layer_name, quiet = TRUE)
dt <- st_drop_geometry(dt_spatial)

has_date_col <- any(c("Date", "date") %in% names(dt))

# --- PRINT AVAILABLE DATES TO CONSOLE ---
if (has_date_col) {
  date_col_name <- if ("Date" %in% names(dt)) "Date" else "date"
  
  available_dates <- dt %>%
    mutate(
      Parsed = parse_date_time(.data[[date_col_name]], orders = c("dmy", "ymd", "mdy")),
      Formatted = format(Parsed, "%d/%m/%Y")
    ) %>%
    filter(!is.na(Parsed)) %>%
    group_by(Formatted, Parsed) %>%
    summarize(Record_Count = n(), .groups = "drop") %>%
    arrange(Parsed)
  
  message("\n--- AVAILABLE SURVEY DATES IN DATASET ---")
  print(as.data.frame(available_dates %>% select(Survey_Date = Formatted, Record_Count)))
  message("-----------------------------------------\n")
}

dt <- dt %>%
  mutate(
    Count           = coalesce(as.numeric(Count), 0),
    Species         = as.character(Species),
    Site            = if ("Site" %in% names(.)) {
      str_to_title(str_squish(coalesce(as.character(Site), "Unspecified")))
    } else "Unspecified",
    Field_Code      = if ("Field.Code" %in% names(.)) coalesce(as.character(`Field.Code`), "Unspecified")
    else if ("Field_Code" %in% names(.)) coalesce(as.character(Field_Code), "Unspecified")
    else if ("Field" %in% names(.)) coalesce(as.character(Field), "Unspecified")
    else "Unspecified",
    Breeding_status = if ("Breeding_status" %in% names(.)) coalesce(as.character(Breeding_status), "")
    else if ("Breeding_sta" %in% names(.)) coalesce(as.character(Breeding_sta), "")
    else if ("Status" %in% names(.)) coalesce(as.character(Status), "")
    else "",
    BTO_code        = if ("BTO_code" %in% names(.)) as.character(BTO_code) 
    else if ("BTO.code" %in% names(.)) as.character(BTO.code) 
    else NA_character_
  ) %>%
  assign_field_group(CONFIG$woodland_codes, CONFIG$farmland_codes)

if (has_date_col) {
  date_col_name <- if ("Date" %in% names(dt)) "Date" else "date"
  dt <- dt %>%
    mutate(
      Date_parsed    = parse_date_time(.data[[date_col_name]], orders = c("dmy", "ymd", "mdy")),
      Date_formatted = format(Date_parsed, "%d/%m/%Y")
    ) %>%
    filter(!is.na(Species), Species != "", !is.na(Date_parsed))
} else {
  dt <- dt %>%
    mutate(
      Date_parsed    = NA,
      Date_formatted = "Count"
    ) %>%
    filter(!is.na(Species), Species != "")
}

# ----------------------------------------------------------------------------
# 2. BOCC5 & TAXONOMIC LOOKUPS - SPECIES ORDERING
# ----------------------------------------------------------------------------

if (file.exists(CONFIG$taxonomic_lookup_path)) {
  tax_lookup <- read.csv(CONFIG$taxonomic_lookup_path, stringsAsFactors = FALSE) %>%
    mutate(
      species     = as.character(species),
      taxon_order = as.numeric(taxon_order)
    ) %>%
    distinct(species, .keep_all = TRUE)
} else {
  tax_lookup <- data.frame(species = character(), taxon_order = numeric())
}

# BoCC5 status is always taken from the species lookup (columns: species, bocc),
# so that changes to the lookup (e.g. the 2024 seabird update) are picked up
# without re-running the data cleaning scripts. The BOCC5 column stored in the
# layer is used only for species that are missing from the lookup.
if (!file.exists(CONFIG$bocc5_lookup_path)) stop("BoCC5 lookup not found: ", CONFIG$bocc5_lookup_path)

bocc5_lookup <- read.csv(CONFIG$bocc5_lookup_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM") %>%
  transmute(
    Species      = str_squish(as.character(species)),
    BOCC5_lookup = str_to_title(str_squish(as.character(bocc)))
  ) %>%
  distinct(Species, .keep_all = TRUE)

if ("BOCC5" %in% names(dt)) {
  layer_bocc <- dt %>%
    transmute(Species, BOCC5_layer = str_to_title(str_squish(as.character(BOCC5)))) %>%
    filter(!is.na(BOCC5_layer), BOCC5_layer != "") %>%
    distinct(Species, .keep_all = TRUE)
} else {
  layer_bocc <- data.frame(Species = character(), BOCC5_layer = character())
}

bocc_data <- dt %>%
  distinct(Species) %>%
  left_join(bocc5_lookup, by = "Species") %>%
  left_join(layer_bocc, by = "Species") %>%
  mutate(BOCC5 = coalesce(BOCC5_lookup, BOCC5_layer))

# Report species whose status in the lookup differs from the layer
changed <- bocc_data %>% filter(!is.na(BOCC5_lookup), !is.na(BOCC5_layer), BOCC5_lookup != BOCC5_layer)
if (nrow(changed) > 0) {
  message("BoCC5 status taken from the lookup where it differs from the layer:")
  print(as.data.frame(changed %>% select(Species, layer = BOCC5_layer, lookup = BOCC5_lookup)))
}
missing <- bocc_data %>% filter(is.na(BOCC5_lookup)) %>% pull(Species)
if (length(missing) > 0) message("Species not in the BoCC5 lookup (layer status used): ", paste(sort(missing), collapse = ", "))

bocc_data <- bocc_data %>% select(Species, BOCC5)

species_order <- dt %>%
  distinct(Species) %>%
  left_join(bocc_data, by = "Species") %>%
  left_join(tax_lookup %>% select(species, taxon_order), by = c("Species" = "species")) %>%
  mutate(
    BOCC5       = factor(BOCC5, levels = CONFIG$bocc5_order),
    taxon_order = coalesce(taxon_order, 9999)
  ) %>%
  distinct(Species, .keep_all = TRUE)

get_highest_breeding_status <- function(df) {
  df %>%
    mutate(
      status_clean = str_to_title(trimws(Breeding_status)),
      status_rank  = case_when(
        str_detect(status_clean, "^Conf") ~ 3,
        str_detect(status_clean, "^Prob") ~ 2,
        str_detect(status_clean, "^Poss") ~ 1,
        TRUE                              ~ 0
      )
    ) %>%
    group_by(Species) %>%
    arrange(desc(status_rank)) %>%
    slice(1) %>%
    ungroup() %>%
    select(Species, Breeding_status = status_clean)
}

# ----------------------------------------------------------------------------
# 3. BUILD SPECIES SUMMARY TABLES (BREAKDOWN TABLE)
# ----------------------------------------------------------------------------

build_species_date_table <- function(dt_all, species_order, config, has_date = TRUE) {
  
  if (nrow(dt_all) == 0) return(NULL)
  
  dt_work_data <- apply_date_filter(dt_all, config$dates_to_include)
  if (nrow(dt_work_data) == 0) return(NULL)
  
  if (!is.null(config$filter_field_group) && config$filter_field_group != "") {
    dt_work_data <- dt_work_data %>% 
      filter(tolower(Field_Group) == tolower(config$filter_field_group))
  }
  
  if (nrow(dt_work_data) == 0) return(NULL)
  
  split_mode   <- tolower(config$split_columns_by)
  use_breeding <- isTRUE(config$include_breeding_status)
  split_col    <- if (split_mode == "site") "Site" else "none"
  
  has_active_date_split <- has_date && !is.null(config$dates_to_include) && length(config$dates_to_include) > 0
  
  if (has_active_date_split) {
    if (split_col != "none" && split_col %in% names(dt_work_data)) {
      col_info <- dt_work_data %>%
        distinct(Date_formatted, Date_parsed, SubGroup = .data[[split_col]]) %>%
        arrange(Date_parsed, SubGroup) %>%
        mutate(Col_Name = paste(Date_formatted, SubGroup, sep = " | "))
      
      data_cols <- col_info$Col_Name
      dt_work   <- dt_work_data %>% mutate(Col_Name = paste(Date_formatted, .data[[split_col]], sep = " | "))
      two_tier_header <- TRUE
    } else {
      date_order <- dt_work_data %>%
        distinct(Date_formatted, Date_parsed) %>%
        arrange(Date_parsed) %>%
        pull(Date_formatted)
      
      data_cols <- date_order
      dt_work   <- dt_work_data %>% mutate(Col_Name = Date_formatted)
      two_tier_header <- FALSE
    }
  } else {
    if (split_col != "none" && split_col %in% names(dt_work_data)) {
      subgroups <- dt_work_data %>% 
        distinct(SubGroup = .data[[split_col]]) %>% 
        arrange(SubGroup) %>% 
        pull(SubGroup)
      data_cols <- subgroups
      dt_work   <- dt_work_data %>% mutate(Col_Name = .data[[split_col]])
      two_tier_header <- FALSE
    } else {
      data_cols <- "Count"
      dt_work   <- dt_work_data %>% mutate(Col_Name = "Count")
      two_tier_header <- FALSE
    }
  }
  
  table_wide <- dt_work %>%
    group_by(Species, Col_Name) %>%
    summarise(count_total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Col_Name, values_from = count_total, values_fill = 0)
  
  data_cols <- intersect(data_cols, names(table_wide))
  if (length(data_cols) == 0) return(NULL)
  
  table_wide <- table_wide %>% select(Species, all_of(data_cols)) %>% filter(if_any(all_of(data_cols), ~ . > 0))
  if (nrow(table_wide) == 0) return(NULL)
  
  if (use_breeding) {
    species_breeding <- get_highest_breeding_status(dt_work_data)
    table_wide <- table_wide %>% left_join(species_breeding, by = "Species") %>% relocate(Breeding_status, .after = Species)
  }
  
  table_wide <- table_wide %>% left_join(species_order, by = "Species")
  
  if (tolower(config$sort_by) == "taxonomic") {
    table_wide <- table_wide %>% arrange(taxon_order, Species)
  } else {
    table_wide <- table_wide %>% arrange(BOCC5, taxon_order, Species)
  }
  
  show_total_col <- length(data_cols) > 1
  if (show_total_col) {
    table_wide <- table_wide %>% mutate(Total = rowSums(pick(all_of(data_cols)), na.rm = TRUE))
    sum_cols <- c(data_cols, "Total")
  } else {
    sum_cols <- data_cols
  }
  
  grand_total <- table_wide %>%
    summarise(across(all_of(sum_cols), \ (x) sum(x, na.rm = TRUE))) %>%
    mutate(Species = "Total", BOCC5 = NA_character_, .before = 1)
  
  if (use_breeding) grand_total$Breeding_status <- ""
  table_wide <- bind_rows(table_wide, grand_total)
  
  lead_cols <- if (use_breeding) c("Species", "Breeding_status") else "Species"
  col_keys  <- if (show_total_col) c(lead_cols, data_cols, "Total") else c(lead_cols, data_cols)
  
  ft <- flextable(table_wide, col_keys = col_keys)
  
  if (two_tier_header) {
    header_df <- data.frame(
      col_keys = col_keys,
      Date_Header = c(if (use_breeding) c("Species", "Breeding Status") else "Species", sub(" \\| .*", "", data_cols), "Total"),
      Sub_Header  = c(if (use_breeding) c("Species", "Breeding Status") else "Species", sub(".* \\| ", "", data_cols), "Total"),
      stringsAsFactors = FALSE
    )
    ft <- set_header_df(ft, mapping = header_df, key = "col_keys") %>% merge_h(part = "header") %>% merge_v(part = "header")
  } else {
    labels <- c("Species", if (use_breeding) "Breeding Status", data_cols, if (show_total_col) "Total")
    ft <- ft %>% set_header_labels(values = setNames(as.list(labels), col_keys))
  }
  
  ft <- ft %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = setdiff(col_keys, "Species"), align = "center", part = "all")
  
  if (use_breeding) {
    ft <- ft %>% flextable::align(j = "Breeding_status", align = "left", part = "all")
  }
  
  ft <- ft %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Red",   color = config$bocc_text_colors["Red"],   part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Amber", color = config$bocc_text_colors["Amber"], part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Green", color = config$bocc_text_colors["Green"], part = "body") %>%
    bold(i = nrow(table_wide)) %>%
    bg(i = nrow(table_wide), bg = "grey85", part = "body") %>%
    color(i = nrow(table_wide), color = "black", part = "body")
  
  if (show_total_col) {
    ft <- ft %>% bold(j = "Total") %>% bg(j = "Total", bg = "grey85", part = "body") %>% color(j = "Total", color = "black", part = "body")
  }
  
  ft <- ft %>%
    hline(i = nrow(table_wide) - 1, border = fp_border(width = 1, style = "solid"), part = "body") %>%
    hline_top(border = fp_border(width = 2, style = "solid"), part = "header") %>%
    hline_bottom(border = fp_border(width = 1, style = "solid"), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = "Species Counts Summary") %>%
    set_table_properties(width = 1, layout = "autofit")
  
  return(ft)
}

# ----------------------------------------------------------------------------
# 4. BUILD MAXIMUM / OVERALL COUNTS TABLE (TABLE 1)
# ----------------------------------------------------------------------------

build_max_counts_table <- function(dt_all, species_order, config) {
  if (nrow(dt_all) == 0) return(NULL)
  
  dt_work_data <- apply_date_filter(dt_all, config$dates_to_include)
  if (nrow(dt_work_data) == 0) return(NULL)
  
  use_breeding <- isTRUE(config$include_breeding_status)
  
  table_max <- dt_work_data %>%
    group_by(Species, Date_formatted) %>%
    summarise(visit_count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(Max_Count = max(visit_count, na.rm = TRUE), .groups = "drop")
  
  if (use_breeding) {
    species_breeding <- get_highest_breeding_status(dt_work_data)
    table_max <- table_max %>% left_join(species_breeding, by = "Species") %>% relocate(Breeding_status, .after = Species)
  }
  
  table_max <- table_max %>% left_join(species_order, by = "Species")
  
  if (tolower(config$sort_by) == "taxonomic") {
    table_max <- table_max %>% arrange(taxon_order, Species)
  } else {
    table_max <- table_max %>% arrange(BOCC5, taxon_order, Species)
  }
  
  lead_cols <- if (use_breeding) c("Species", "Breeding_status") else "Species"
  col_keys  <- c(lead_cols, "Max_Count")
  
  ft <- flextable(table_max, col_keys = col_keys) %>%
    set_header_labels(values = setNames(as.list(c("Species", if (use_breeding) "Breeding Status", "Max Count")), col_keys)) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = "Max_Count", align = "center", part = "all")
  
  if (use_breeding) ft <- ft %>% flextable::align(j = "Breeding_status", align = "left", part = "all")
  
  ft <- ft %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Red",   color = config$bocc_text_colors["Red"],   part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Amber", color = config$bocc_text_colors["Amber"], part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Green", color = config$bocc_text_colors["Green"], part = "body") %>%
    hline_top(border = fp_border(width = 2, style = "solid"), part = "header") %>%
    hline_bottom(border = fp_border(width = 1, style = "solid"), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = "Maximum Species Counts Across Selected Visits") %>%
    set_table_properties(width = 1, layout = "autofit")
  
  return(ft)
}

# ----------------------------------------------------------------------------
# 4B. BUILD MAXIMUM COUNTS BREAKDOWN BY SITE (TABLE 1b)
# ----------------------------------------------------------------------------

build_max_counts_by_site_table <- function(dt_all, species_order, config) {
  if (nrow(dt_all) == 0) return(NULL)
  
  dt_work_data <- apply_date_filter(dt_all, config$dates_to_include)
  if (nrow(dt_work_data) == 0) return(NULL)
  
  use_breeding <- isTRUE(config$include_breeding_status)
  
  table_max_site <- dt_work_data %>%
    group_by(Species, Site, Date_formatted) %>%
    summarise(visit_count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(Species, Site) %>%
    summarise(Max_Count = max(visit_count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Site, values_from = Max_Count, values_fill = 0)
  
  site_cols <- setdiff(names(table_max_site), "Species")
  if (length(site_cols) == 0) return(NULL)
  
  table_max_site <- table_max_site %>% filter(if_any(all_of(site_cols), ~ . > 0))
  if (nrow(table_max_site) == 0) return(NULL)
  
  if (use_breeding) {
    species_breeding <- get_highest_breeding_status(dt_work_data)
    table_max_site <- table_max_site %>% left_join(species_breeding, by = "Species") %>% relocate(Breeding_status, .after = Species)
  }
  
  table_max_site <- table_max_site %>% left_join(species_order, by = "Species")
  
  if (tolower(config$sort_by) == "taxonomic") {
    table_max_site <- table_max_site %>% arrange(taxon_order, Species)
  } else {
    table_max_site <- table_max_site %>% arrange(BOCC5, taxon_order, Species)
  }
  
  lead_cols <- if (use_breeding) c("Species", "Breeding_status") else "Species"
  col_keys  <- c(lead_cols, site_cols)
  labels    <- c("Species", if (use_breeding) "Breeding Status", site_cols)
  
  ft <- flextable(table_max_site, col_keys = col_keys) %>%
    set_header_labels(values = setNames(as.list(labels), col_keys)) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = setdiff(col_keys, "Species"), align = "center", part = "all")
  
  if (use_breeding) ft <- ft %>% flextable::align(j = "Breeding_status", align = "left", part = "all")
  
  ft <- ft %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Red",   color = config$bocc_text_colors["Red"],   part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Amber", color = config$bocc_text_colors["Amber"], part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Green", color = config$bocc_text_colors["Green"], part = "body") %>%
    hline_top(border = fp_border(width = 2, style = "solid"), part = "header") %>%
    hline_bottom(border = fp_border(width = 1, style = "solid"), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = "Maximum Species Counts Breakdown by Site") %>%
    set_table_properties(width = 1, layout = "autofit")
  
  return(ft)
}

# ----------------------------------------------------------------------------
# 4C. BUILD MAXIMUM COUNTS BREAKDOWN BY SITE & FIELD GROUP (TABLE 1c)
# ----------------------------------------------------------------------------

build_max_counts_by_site_and_group_table <- function(dt_all, species_order, config) {
  if (nrow(dt_all) == 0) return(NULL)
  
  dt_work_data <- apply_date_filter(dt_all, config$dates_to_include)
  if (nrow(dt_work_data) == 0) return(NULL)
  
  if (!is.null(config$filter_field_group) && config$filter_field_group != "") {
    dt_work_data <- dt_work_data %>% 
      filter(tolower(Field_Group) == tolower(config$filter_field_group))
  }
  if (nrow(dt_work_data) == 0) return(NULL)
  
  use_breeding <- isTRUE(config$include_breeding_status)
  
  table_max_site_group <- dt_work_data %>%
    group_by(Species, Site, Date_formatted) %>%
    summarise(visit_count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(Species, Site) %>%
    summarise(Max_Count = max(visit_count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Site, values_from = Max_Count, values_fill = 0)
  
  site_cols <- setdiff(names(table_max_site_group), "Species")
  if (length(site_cols) == 0) return(NULL)
  
  table_max_site_group <- table_max_site_group %>% filter(if_any(all_of(site_cols), ~ . > 0))
  if (nrow(table_max_site_group) == 0) return(NULL)
  
  if (use_breeding) {
    species_breeding <- get_highest_breeding_status(dt_work_data)
    table_max_site_group <- table_max_site_group %>% left_join(species_breeding, by = "Species") %>% relocate(Breeding_status, .after = Species)
  }
  
  table_max_site_group <- table_max_site_group %>% left_join(species_order, by = "Species")
  
  if (tolower(config$sort_by) == "taxonomic") {
    table_max_site_group <- table_max_site_group %>% arrange(taxon_order, Species)
  } else {
    table_max_site_group <- table_max_site_group %>% arrange(BOCC5, taxon_order, Species)
  }
  
  lead_cols <- if (use_breeding) c("Species", "Breeding_status") else "Species"
  col_keys  <- c(lead_cols, site_cols)
  labels    <- c("Species", if (use_breeding) "Breeding Status", site_cols)
  
  ft <- flextable(table_max_site_group, col_keys = col_keys) %>%
    set_header_labels(values = setNames(as.list(labels), col_keys)) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = setdiff(col_keys, "Species"), align = "center", part = "all")
  
  if (use_breeding) ft <- ft %>% flextable::align(j = "Breeding_status", align = "left", part = "all")
  
  ft <- ft %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Red",   color = config$bocc_text_colors["Red"],   part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Amber", color = config$bocc_text_colors["Amber"], part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Green", color = config$bocc_text_colors["Green"], part = "body") %>%
    hline_top(border = fp_border(width = 2, style = "solid"), part = "header") %>%
    hline_bottom(border = fp_border(width = 1, style = "solid"), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = paste0("Maximum Species Counts Breakdown by Site (", config$filter_field_group, ")")) %>%
    set_table_properties(width = 1, layout = "autofit")
  
  return(ft)
}

# ----------------------------------------------------------------------------
# 5. GENERATE WORD DOCUMENT WITH TABLES
# ----------------------------------------------------------------------------

doc <- if (!is.null(CONFIG$template_path)) read_docx(CONFIG$template_path) else read_docx()

# Table 1: Maximum Counts per Species (site-wide)
ft_max <- build_max_counts_table(dt, species_order, CONFIG)
if (!is.null(ft_max)) {
  doc <- doc %>%
    body_add_par("1. Maximum Counts per Species Across Visits", style = "heading 1") %>%
    body_add_flextable(ft_max)
}

# Table 1b: Maximum Counts Breakdown by Site (site-wide, all field groups)
ft_max_site <- build_max_counts_by_site_table(dt, species_order, CONFIG)
if (!is.null(ft_max_site)) {
  doc <- doc %>%
    body_add_break() %>%
    body_add_par("1b. Maximum Counts per Species by Site", style = "heading 1") %>%
    body_add_flextable(ft_max_site)
}

# Table 2: Counts per Visit (Site-wide)
if (has_date_col) {
  cfg_visit <- CONFIG
  cfg_visit$filter_field_group <- NULL 
  ft_visit <- build_species_date_table(dt, species_order, cfg_visit, has_date = TRUE)
  if (!is.null(ft_visit)) {
    doc <- doc %>%
      body_add_break() %>%
      body_add_par("2. Species Counts per Visit", style = "heading 1") %>%
      body_add_flextable(ft_visit)
  }
}

# --- DYNAMIC LOOPS FOR EACH FIELD GROUP TYPE ---
all_field_groups <- sort(unique(dt$Field_Group))
all_field_groups <- all_field_groups[!is.na(all_field_groups) & all_field_groups != "Unspecified"]

section_num <- 3 

for (grp in all_field_groups) {
  
  cfg_loop <- CONFIG
  cfg_loop$filter_field_group <- grp
  
  # Max counts by site for this field group
  ft_max_site_group <- build_max_counts_by_site_and_group_table(dt, species_order, cfg_loop)
  if (!is.null(ft_max_site_group)) {
    doc <- doc %>%
      body_add_break() %>%
      body_add_par(paste0("Maximum Counts per Species by Site (", grp, ")"), style = "heading 1") %>%
      body_add_flextable(ft_max_site_group)
  }
  
  # Detailed Breakdown for this field group
  ft_breakdown <- build_species_date_table(dt, species_order, cfg_loop, has_date = has_date_col)
  if (!is.null(ft_breakdown)) {
    doc <- doc %>%
      body_add_break() %>%
      body_add_par(paste0(section_num, ". Species Counts Breakdown (", grp, ")"), style = "heading 1") %>%
      body_add_flextable(ft_breakdown)
    
    section_num <- section_num + 1
  }
}

# Dynamic Export Name including GeoPackage layer name
file_name <- paste0("survey_summary_tables_", CONFIG$layer_name, ".docx")
out_file  <- file.path(CONFIG$output_dir, file_name)

print(doc, target = out_file)
message("Saved summary tables to: ", out_file)