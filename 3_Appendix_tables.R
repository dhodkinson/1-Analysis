#updates
# need to change non-native species to black text

rm(list=ls()) #clear all previous

# ============================================================================
# SPECIES x FIELD.CODE APPENDIX TABLE GENERATOR (GEOPACKAGE & TAXONOMIC)
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(sf)
library(flextable)
library(officer)
library(lubridate)

# ----------------------------------------------------------------------------
# 0. CONFIG - EDIT THESE
# ----------------------------------------------------------------------------

CONFIG <- list(
  input_path              = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg',
  layer_name              = '2026_Count_final',
  bocc5_lookup_path       = "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv",
  taxonomic_lookup_path   = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BBS_species taxonomic order lookup.csv',
  output_dir              = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results',
  bocc5_order             = c("Red", "Amber", "Green"),
  
  dates_to_include        = NULL,
  filter_field_group      = "Woodland", 
  split_columns_by        = "site",
  
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
# 1. LOAD & STANDARDISE DATES FROM GEOPACKAGE
# ----------------------------------------------------------------------------

if (!is.null(CONFIG$layer_name)) {
  dt <- st_read(CONFIG$input_path, layer = CONFIG$layer_name, quiet = TRUE)
} else {
  layers <- st_layers(CONFIG$input_path)
  dt <- st_read(CONFIG$input_path, layer = layers$name[1], quiet = TRUE)
}

dt <- as.data.frame(dt) %>% select(-any_of("geom"))

has_date_col <- any(c("Date", "date") %in% names(dt))

dt <- dt %>%
  mutate(
    Count      = coalesce(as.numeric(Count), 0),
    Field.Code = as.character(Field.Code),
    Species    = as.character(Species),
    BTO_code   = if ("BTO_code" %in% names(.)) as.character(BTO_code) 
    else if ("BTO.code" %in% names(.)) as.character(BTO.code) 
    else NA_character_
  ) %>%
  filter(!is.na(Field.Code), Field.Code != "",
         !is.na(Species), Species != "")

if (has_date_col) {
  date_col_name <- if ("Date" %in% names(dt)) "Date" else "date"
  
  dt <- dt %>%
    mutate(
      raw_date       = .data[[date_col_name]],
      Parsed_Date    = parse_date_time(raw_date, orders = c("dmy", "ymd", "mdy", "dOmY")),
      Date_formatted = format(Parsed_Date, "%d/%m/%Y")
    ) %>%
    filter(!is.na(Parsed_Date)) %>%
    select(-raw_date, -Parsed_Date)
  
} else {
  message("Note: 'Date' column not found in data. Assigning default 'Summer_Survey'.")
  dt <- dt %>%
    mutate(
      Date_formatted = "Summer_Survey"
    )
}

apply_date_filter <- function(df, dates_to_include) {
  if (is.null(dates_to_include) || length(dates_to_include) == 0) {
    return(df)
  }
  parsed_config_dates <- parse_date_time(dates_to_include, orders = c("dmy", "ymd", "mdy", "dOmY"))
  std_config_dates <- format(parsed_config_dates, "%d/%m/%Y")
  
  df %>% filter(Date_formatted %in% std_config_dates)
}

dt <- apply_date_filter(dt, CONFIG$dates_to_include)

# ----------------------------------------------------------------------------
# 2. SPLIT FIELD.CODE INTO NUMERIC (1-8) / A-L / M-X GROUPS
# ----------------------------------------------------------------------------

assign_field_group <- function(field_code) {
  first_char <- toupper(substr(field_code, 1, 1))
  case_when(
    first_char %in% as.character(1:8)     ~ "Numeric_1_8",
    first_char >= "A" & first_char <= "L" ~ "A_L",
    TRUE                                  ~ "M_X"
  )
}

dt <- dt %>% mutate(field_group = assign_field_group(Field.Code))

# ----------------------------------------------------------------------------
# 3. LOOKUPS: BOCC5 & TAXONOMIC ORDERING
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
    BOCC5         = factor(BOCC5, levels = CONFIG$bocc5_order),
    taxon_order   = coalesce(taxon_order, 9999)
  ) %>%
  distinct(Species, .keep_all = TRUE)

# ----------------------------------------------------------------------------
# 4. TABLE-BUILDING FUNCTION
# ----------------------------------------------------------------------------

build_species_table <- function(dt_sub, group_label, species_order, config) {
  
  if (nrow(dt_sub) == 0) return(NULL)
  
  table_wide <- dt_sub %>%
    group_by(Species, Field.Code) %>%
    summarise(count_total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Field.Code,
                values_from = count_total,
                values_fill = 0)
  
  field_cols <- sort(setdiff(names(table_wide), "Species"))
  if (length(field_cols) == 0) return(NULL)
  
  table_wide <- table_wide %>% select(Species, all_of(field_cols))
  table_wide <- table_wide %>% filter(if_any(all_of(field_cols), ~ . > 0))
  
  if (nrow(table_wide) == 0) return(NULL)
  
  table_wide <- table_wide %>%
    left_join(species_order, by = "Species")
  
  # Dynamic Sorting depending on configuration flag
  if (isTRUE(config$sort_by_bocc)) {
    table_wide <- table_wide %>%
      mutate(bocc_sort_val = as.integer(BOCC5)) %>% # Red=1, Amber=2, Green=3, NA=4
      mutate(bocc_sort_val = coalesce(bocc_sort_val, 99)) %>%
      arrange(bocc_sort_val, taxon_order, Species) %>%
      select(-bocc_sort_val)
  } else {
    table_wide <- table_wide %>%
      arrange(taxon_order, Species)
  }
  
  table_wide <- table_wide %>%
    mutate(Total = rowSums(pick(all_of(field_cols)), na.rm = TRUE))
  
  grand_total <- table_wide %>%
    summarise(across(all_of(c(field_cols, "Total")), \ (x) sum(x, na.rm = TRUE))) %>%
    mutate(Species = "Total", taxon_order = 99999, BOCC5 = NA_character_, .before = 1)
  
  table_wide <- bind_rows(table_wide, grand_total)
  
  col_keys <- c("Species", field_cols, "Total")
  
  ft <- flextable(table_wide, col_keys = col_keys) %>%
    set_header_labels(
      values = setNames(as.list(c("Species", field_cols, "Total")), col_keys)
    ) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = c(field_cols, "Total"), align = "center", part = "all") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Red",   color = config$bocc_text_colors["Red"],   part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Amber", color = config$bocc_text_colors["Amber"], part = "body") %>%
    color(i = ~ !is.na(BOCC5) & BOCC5 == "Green", color = config$bocc_text_colors["Green"], part = "body") %>%
    bold(i = nrow(table_wide)) %>%
    bold(j = "Total") %>%
    bg(j = "Total", bg = "grey85", part = "body") %>%
    bg(i = nrow(table_wide), bg = "grey85", part = "body") %>%
    color(i = nrow(table_wide), color = "black", part = "body") %>%
    color(j = "Total", color = "black", part = "body") %>%
    hline(i = nrow(table_wide) - 1, 
          border = fp_border(width = 1, style = "solid"), part = "body") %>%
    hline_top(border = fp_border(width = 2, style = "solid"), part = "header") %>%
    hline_bottom(border = fp_border(width = 1, style = "solid"), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = group_label) %>%
    set_table_properties(width = 1, layout = "autofit")
  
  return(ft)
}

# ----------------------------------------------------------------------------
# 5. GENERATE WORD DOCUMENT FOR EACH DATE
# ----------------------------------------------------------------------------

generate_tables_for_date <- function(one_date, dt_all, config, species_order) {
  
  dt_date <- dt_all %>% filter(Date_formatted == one_date)
  
  ft_num <- build_species_table(dt_date %>% filter(field_group == "Numeric_1_8"),
                                "Field Codes 1-8", species_order, config)
  ft_al  <- build_species_table(dt_date %>% filter(field_group == "A_L"),
                                "Field Codes A-L", species_order, config)
  ft_mx  <- build_species_table(dt_date %>% filter(field_group == "M_X"),
                                "Field Codes M-X", species_order, config)
  
  if (is.null(ft_num) && is.null(ft_al) && is.null(ft_mx)) {
    message("No data for date: ", one_date, " - skipping")
    return(invisible(NULL))
  }
  
  doc <- if (!is.null(config$template_path)) {
    read_docx(config$template_path)
  } else {
    read_docx()
  }
  
  doc <- doc %>%
    body_add_par(paste("Species counts by Field Code -", one_date),
                 style = "heading 1")
  
  if (!is.null(ft_num)) {
    doc <- doc %>%
      body_add_par("Woodlands 1-8", style = "heading 2") %>%
      body_add_flextable(ft_num) %>%
      body_add_par("")
  }
  
  if (!is.null(ft_al)) {
    doc <- doc %>%
      body_add_par("Field Codes A-L", style = "heading 2") %>%
      body_add_flextable(ft_al) %>%
      body_add_par("")
  }
  
  if (!is.null(ft_mx)) {
    doc <- doc %>%
      body_add_par("Field Codes M-X", style = "heading 2") %>%
      body_add_flextable(ft_mx)
  }
  
  safe_date <- gsub("[^A-Za-z0-9_-]", "_", one_date)
  out_file <- file.path(config$output_dir,
                        paste0("appendix_table_counts_", safe_date, ".docx"))
  print(doc, target = out_file)
  message("Saved: ", out_file)
}

# ----------------------------------------------------------------------------
# 6. RUN SCRIPT
# ----------------------------------------------------------------------------

all_dates <- sort(unique(dt$Date_formatted))

walk(all_dates, generate_tables_for_date,
     dt_all = dt, config = CONFIG, species_order = species_order)