## =================================================================
## Prepare raw bird survey data for analysis in R (Single Layer Version)
## =================================================================

rm(list=ls()) #NEVER REMOVE

## -----------------------------------------------------------------
## Package check -- installs anything missing, then loads it
## -----------------------------------------------------------------
required_packages <- c("sf", "dplyr", "readr", "stringr", "stringdist")

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(sf)
library(dplyr)
library(readr)
library(stringr)
library(data.table)

## -----------------------------------------------------------------
## 0. File paths -- EDIT THESE to match your setup
## -----------------------------------------------------------------
gpkg_path          <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'
input_layer_name   <- "2026_Territories"
# Species lookup with columns: species, bto_code, taxon_order, sci_name, bocc,
# winter_visitor. These are renamed to Code, Species and BOCC5 when read in.
lookup_path        <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv"
status_lookup_path <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/Species_status_codes.csv"
output_dir         <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning'

output_layer_name  <- "2026_Territories_clean"
notes_layer_name   <- "Shimpling_2026_Notes_for_review"

## -----------------------------------------------------------------
## Project-specific output folder & file definitions
## -----------------------------------------------------------------
project_name <- "Shimpling_2026"

project_output_dir <- file.path(output_dir, project_name)
if (!dir.exists(project_output_dir)) {
  dir.create(project_output_dir, recursive = TRUE)
}

categorised_path <- file.path(
  project_output_dir,
  paste0(project_name, "_bto_code_categorisation.csv")
)

notes_csv_path <- file.path(
  project_output_dir,
  paste0(project_name, "_Notes_for_review.csv")
)

final_csv_path <- file.path(
  project_output_dir,
  paste0(output_layer_name, ".csv")
)

## -----------------------------------------------------------------
## 1. Read the raw layer and standardise column names/types
## -----------------------------------------------------------------
cat("Reading layer '", input_layer_name, "' from geopackage...\n", sep = "")
raw_data <- st_read(gpkg_path, layer = input_layer_name, quiet = TRUE)

combined_data <- raw_data %>%
  mutate(
    record_id   = paste(input_layer_name, row_number(), sep = "_"),
    BTO_code    = as.character(.[[match("BTO.code", names(.))]]),
    Status      = as.character(.[[match("BoCC.status", names(.))]]),
    Territories = coalesce(Territories, 1)  # change NA to 1, keep others
  )

## -----------------------------------------------------------------
## 2. Capitalise BTO codes
## -----------------------------------------------------------------
combined_data <- combined_data %>%
  mutate(BTO_code = str_to_upper(str_trim(BTO_code)))

## -----------------------------------------------------------------
## 3. Check BTO codes and status codes against lookups
## -----------------------------------------------------------------
species_lookup_raw <- read_csv(lookup_path, show_col_types = FALSE) %>%
  rename(any_of(c(Code = "bto_code", Species = "species", BOCC5 = "bocc")))
species_lookup <- species_lookup_raw %>%
  select(Code, Species, BOCC5) %>%
  mutate(Code = str_to_upper(str_trim(Code)))

status_lookup_raw <- read_csv(status_lookup_path, show_col_types = FALSE)
status_lookup <- status_lookup_raw %>%
  mutate(Status_code = str_to_upper(str_trim(Status_code)))

blank_like_status <- c("", "NA", "N/A", "NULL", "NONE", "-")

combined_data <- combined_data %>%
  mutate(
    code_valid     = BTO_code %in% species_lookup$Code,
    Status_compare = str_to_upper(str_trim(Status)),
    status_valid   = is.na(Status_compare) |
      Status_compare %in% blank_like_status |
      Status_compare %in% status_lookup$Status_code
  )

## -----------------------------------------------------------------
## Categorise invalid BTO codes / status codes
## -----------------------------------------------------------------
needs_review <- any(!combined_data$code_valid) || any(!combined_data$status_valid)

if (needs_review) {
  combined_data <- combined_data %>%
    mutate(needs_review = !code_valid | !status_valid)
  
  invalid_records <- combined_data %>%
    st_drop_geometry() %>%
    filter(needs_review) %>%
    select(record_id, BTO_code, code_valid, Status, status_valid, Site)
  
  if (!file.exists(categorised_path)) {
    invalid_records %>%
      mutate(category = NA_integer_, corrected_BTO_code = NA_character_, corrected_status = NA_character_) %>%
      write_csv(categorised_path)
    
    stop(
      "A categorisation template has been written to:\n", categorised_path,
      "\n\nOpen it, fill out the 'category' (1-4) and correction columns, save, and re-run."
    )
  }
  
  raw_categorised <- read_csv(categorised_path, show_col_types = FALSE)
  if ("...1" %in% names(raw_categorised) && !"record_id" %in% names(raw_categorised)) {
    raw_categorised <- raw_categorised %>% rename(record_id = ...1)
  }
  raw_categorised <- raw_categorised %>% mutate(record_id = as.character(record_id))
  
  record_categories <- raw_categorised %>%
    select(record_id, category, corrected_BTO_code, corrected_status)
  
  if (any(is.na(record_categories$category))) {
    stop("Some rows in categorisation CSV still have a blank 'category'. Fill in 1, 2, 3, or 4 for all rows.")
  }
  
  combined_data <- combined_data %>%
    left_join(record_categories, by = "record_id")
  
  export_data <- combined_data %>% filter(category == 1)
  if (nrow(export_data) > 0) {
    write_csv(st_drop_geometry(export_data), notes_csv_path)
    st_write(export_data, dsn = gpkg_path, layer = notes_layer_name, delete_layer = TRUE, quiet = TRUE)
  }
  
  combined_data <- combined_data %>%
    mutate(
      BTO_code = if_else(!is.na(category) & category == 2 & !code_valid, str_to_upper(str_trim(corrected_BTO_code)), BTO_code),
      Status   = if_else(!is.na(category) & category == 2 & !status_valid, str_to_upper(str_trim(corrected_status)), Status),
      code_valid     = BTO_code %in% species_lookup$Code,
      Status_compare = str_to_upper(str_trim(Status)),
      status_valid   = is.na(Status_compare) | Status_compare %in% blank_like_status | Status_compare %in% status_lookup$Status_code
    )
  
  combined_data <- combined_data %>%
    mutate(manually_accepted = !is.na(category) & category == 4) %>%
    filter(is.na(category) | category %in% c(2, 4)) %>%
    select(-any_of(c("category", "corrected_BTO_code", "corrected_status", "needs_review")))
  
} else {
  cat("All BTO codes and status codes matched their respective lookup files.\n")
  combined_data <- combined_data %>% mutate(manually_accepted = FALSE)
}

## -----------------------------------------------------------------
## Bring in human-readable lookup info
## -----------------------------------------------------------------
combined_data <- combined_data %>%
  left_join(species_lookup, by = c("BTO_code" = "Code")) %>%
  left_join(status_lookup, by = c("Status_compare" = "Status_code")) %>%
  select(-Status_compare)


# set name of Territories column to count to aid compatibility with other scripts
setnames(combined_data, "Territories", "Count")


## -----------------------------------------------------------------
## Save the final result into the geopackage
## -----------------------------------------------------------------
st_write(
  combined_data,
  dsn          = gpkg_path,
  layer        = output_layer_name,
  delete_layer = TRUE,
  quiet        = TRUE
)

write_csv(st_drop_geometry(combined_data), final_csv_path)

cat(
  "Final processed data successfully saved to:\n",
  "  - Geopackage layer: '", output_layer_name, "' in '", gpkg_path, "'\n",
  "  - CSV file: '", final_csv_path, "'\n",
  sep = ""
)