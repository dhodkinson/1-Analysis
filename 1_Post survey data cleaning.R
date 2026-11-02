## =================================================================
## Prepare raw bird survey data (geopackage) for analysis in R
## =================================================================
##
## WHAT THIS SCRIPT DOES
## ----------------------
## Takes several raw survey layers sitting inside a QGIS geopackage
## (.gpkg), stitches them together into one table, cleans up two
## coded columns (BTO species code and breeding-status code) against
## reference lookup files, standardises dates into ISO format, 
## and writes the cleaned result back into the geopackage as a brand 
## new layer -- leaving all of the original layers untouched.
##
## THE BIG PICTURE, STEP BY STEP
## ------------------------------
##   1. List every layer in the geopackage and combine a chosen subset
##      of them (e.g. two Winter Bird Survey visits + two Breeding
##      Bird Survey visits) into a single data frame.
##   2. Standardise the BTO species code column (uppercase, trimmed and add BoCC)
##   3. Standardise date formats into standard ISO (YYYY-MM-DD).
##   4. Validate BOTH the BTO species code and the status code for
##      every record against two reference CSVs:
##        - BTO_species_codes_lookup.csv   (valid species codes)
##        - Species_status_codes.csv (valid breeding-status codes)
##      Any record that fails either check gets routed through a
##      manual review process (see below) before the script will
##      let you save the final output.
##   5. Write the cleaned, validated data into the geopackage as a
##      new layer.
##
## THE MANUAL REVIEW PROCESS (the trickiest part of the script)
## ---------------------------------------------------------------
## Real survey data will always contain some typos or genuinely novel
## codes that aren't in the lookup files. Rather than silently
## dropping or guessing at these, the script pauses and asks YOU to
## make a decision for each problem record, via a project-specific CSV file
## (e.g., "[project_name]_bto_code_categorisation.csv") inside a project folder
## rather than interactive console prompts.
##
## The workflow is:
##   - RUN 1: no records have been checked yet. The script finds every
##     record with an invalid BTO code and/or an invalid status code,
##     writes them out to project_output_dir, and then deliberately
##     STOPS (via stop()) so you can go and fill it in.
##   - YOU fill in a 'category' column for every row:
##        1 = remove this record from the main dataset, but keep a
##            copy of it (exported to [project]_Notes_for_review.csv and to its
##            own geopackage layer) in case it's needed later
##        2 = the record is fine, it just needs correcting -- also
##            fill in corrected_BTO_code and/or corrected_status,
##            whichever field(s) were actually invalid for that row.
##            corrected_status may be left blank/NA on purpose, which
##            means "clear the status entirely" (a blank status is valid)
##        3 = genuinely bad/duplicate record -- discard it completely,
##            don't keep a copy anywhere
##        4 = keep the record exactly as it is -- no correction is applied
##            or required, it stays in the final output unchanged even
##            though its BTO code and/or status didn't match the lookup
##   - RUN 2 (re-running the script after editing the CSV): the script
##     reads your decisions back in, applies them (removes / corrects
##     / discards / keeps-as-is as instructed), re-checks any corrected
##     codes against the lookup files, and then proceeds to write the
##     final output.
##   - If you add more raw data later and re-run, any NEW problem
##     records are merged into the same categorisation CSV (previous
##     decisions you made are preserved), so you only ever have to
##     categorise each problem record once.
##
## FILES THIS SCRIPT READS
## -------------------------
##   - gpkg_path        : the geopackage containing your raw survey layers
##   - lookup_path      : CSV of valid BTO species codes (column: Code)
##   - status_lookup_path : CSV of valid status codes
##                        (columns: Status_code, Breeding_status)
##
## FILES THIS SCRIPT WRITES
## -----------------------------------------------------------------
##   Inside project_output_dir ( output_dir / project_name / ):
##   - [project_name]_bto_code_categorisation.csv : the manual review file
##   - [project_name]_Notes_for_review.csv        : category-1 records CSV
##
##   Inside geopackage (gpkg_path):
##   - layer = notes_layer_name  : category-1 records with geometry
##   - layer = output_layer_name : final cleaned dataset
## =================================================================

rm(list=ls()) #NEVER REMOVE

## -----------------------------------------------------------------
## Package check -- installs anything missing, then loads it
## -----------------------------------------------------------------
required_packages <- c("sf", "dplyr", "readr", "stringr", "lubridate", "stringdist")

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages)
}

library(sf)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

## -----------------------------------------------------------------
## 0. File paths -- EDIT THESE to match your setup
## -----------------------------------------------------------------
# The geopackage holding your raw digitised survey layers.
gpkg_path          <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/Broxted/Broxted Summer/Broxted.gpkg'

# Reference list of valid BTO species codes.
# Species lookup with columns: species, bto_code, taxon_order, sci_name, bocc,
# winter_visitor. These are renamed to Code, Species and BOCC5 when read in.
lookup_path        <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv"

# Reference list of valid breeding-status codes. Must contain columns
# 'Status_code' and 'Breeding_status'.
status_lookup_path <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/Species_status_codes.csv"

# Base folder for outputs
output_dir <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/Broxted/Species specific surveys/Turtle Dove'

# Name of the NEW layer the final processed data will be saved as inside gpkg_path
output_layer_name <- "Turtle_Dove_woodland_survey_2026_clean"

# Name of the layer the category 1 (removed/for review) records will be saved as
notes_layer_name <- "Turtle_Dove_2026_Notes_for_review"



## -----------------------------------------------------------------
## Project-specific output folder & file definitions
## -----------------------------------------------------------------
# Make a filesystem-safe version of the output layer name.
project_name <- "Results_2026"

# Store all CSV outputs for this project in their own folder.
project_output_dir <- file.path(output_dir, project_name)
if (!dir.exists(project_output_dir)) {
  dir.create(project_output_dir, recursive = TRUE)
}

# Define project-specific CSV paths:
categorised_path <- file.path(
  project_output_dir,
  paste0(project_name, "_bto_code_categorisation.csv")
)

site_categorised_path <- file.path(
  project_output_dir,
  paste0(project_name, "_site_categorisation.csv")
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
## 1. List layers in the geopackage and select which to combine
## -----------------------------------------------------------------
layers_available <- st_layers(gpkg_path)

cat("Layers found in geopackage:\n")
print(layers_available)

# Fixed set of layers to combine (edit this list if the layers change).
selected_layers <- c("Turtle_Dove_woodland_survey_2026")

# Sanity check: make sure these names actually exist in the geopackage.
missing_layers <- setdiff(selected_layers, layers_available$name)
if (length(missing_layers) > 0) {
  stop(
    "The following requested layers were not found in the geopackage:\n",
    paste(missing_layers, collapse = ", "),
    "\nLayers available are:\n",
    paste(layers_available$name, collapse = ", ")
  )
}

if (length(selected_layers) == 0) {
  stop("No layers selected -- please choose at least one layer to proceed.")
}

cat("Combining layers:", paste(selected_layers, collapse = ", "), "\n")

## -----------------------------------------------------------------
## Read and combine the selected layers
## -----------------------------------------------------------------
read_one_layer <- function(layer, path) {
  this_layer <- st_read(path, layer = layer, quiet = TRUE)
  
  # Find any column that looks like a BTO code column
  raw_bto_col <- names(this_layer)[
    str_detect(str_to_lower(names(this_layer)), "^bto[ ._]?code$")
  ]
  
  if (length(raw_bto_col) == 0) {
    stop("Could not find a BTO code column in layer '", layer, "'.")
  }
  if (length(raw_bto_col) > 1) {
    stop("Found more than one possible BTO code column in layer '", layer, "'.")
  }
  
  cat("Layer '", layer, "': using '", raw_bto_col, "' as the BTO code column.\n", sep = "")
  
  # Rename to canonical form temporarily
  names(this_layer)[names(this_layer) == raw_bto_col] <- "BTO_code"
  
  # Add metadata + keep Date as character initially for parsing
  this_layer <- this_layer %>%
    mutate(
      source_layer = layer,
      record_id    = paste(layer, row_number(), sep = "_"),
      Date         = as.character(Date)
    )
  
  return(this_layer)
}

# ACTUALLY COMBINE THE LAYERS INTO combined_data HERE:
combined_data <- bind_rows(lapply(selected_layers, function(l) {
  read_one_layer(l, gpkg_path)
}))

### -----------------------------------------------------------------
## Standardise Date Format
## -----------------------------------------------------------------
# Safely parse date-times or dates across various potential formats,
# stripping timestamps and converting everything to a standard YYYY-MM-DD string.
combined_data <- combined_data %>%
  mutate(
    Parsed_Date = suppressWarnings(
      parse_date_time(Date, orders = c("ymd HMS", "ymd HM", "ymd", "dmy HMS", "dmy", "mdy", "Ymd"))
    ),
    Date = if_else(!is.na(Parsed_Date), as.character(as.Date(Parsed_Date)), Date)
  ) %>%
  select(-Parsed_Date)

# Check for any unparseable dates
unparseable_dates <- combined_data %>% 
  filter(is.na(Date) | Date == "" | str_detect(Date, "[A-Za-z]"))

if (nrow(unparseable_dates) > 0) {
  warning(
    "Found ", nrow(unparseable_dates), 
    " record(s) with missing, blank, or unparseable dates. Check raw layers if needed."
  )
}

## -----------------------------------------------------------------
## 2. Capitalise BTO codes
## -----------------------------------------------------------------
combined_data <- combined_data %>%
  mutate(BTO_code = str_to_upper(str_trim(BTO_code)))

## -----------------------------------------------------------------
## Site Name Auto-Check & Manual Review Workflow
## -----------------------------------------------------------------
# Keep original Site casing intact, but create an uppercase tracking version
combined_data <- combined_data %>%
  mutate(
    Site_original = Site,
    Site_clean    = str_to_upper(str_trim(Site))
  )

# Get a frequency summary using the ORIGINAL casing for the report/CSV
site_summary <- combined_data %>%
  st_drop_geometry() %>%
  group_by(Site_original) %>%
  summarise(
    record_count = n(),
    Site_clean   = first(Site_clean),
    .groups = "drop"
  ) %>%
  select(Site_original, Site_clean, record_count) %>%
  arrange(desc(record_count))

# If a site categorisation CSV already exists from a previous run, load corrections
if (file.exists(site_categorised_path)) {
  saved_site_reviews <- read_csv(site_categorised_path, show_col_types = FALSE)
  
  # Join decisions back to main dataset if corrections were made
  if ("corrected_site" %in% names(saved_site_reviews)) {
    
    dupe_site_clean <- saved_site_reviews %>%
      count(Site_clean) %>%
      filter(n > 1)
    
    if (nrow(dupe_site_clean) > 0) {
      stop(
        "'", site_categorised_path, "' contains duplicate Site_clean entries, ",
        "which would duplicate matching records when joined: ",
        paste(dupe_site_clean$Site_clean, collapse = ", "),
        "\nPlease remove the duplicate row(s) (one row per site) and re-run."
      )
    }
    
    combined_data <- combined_data %>%
      left_join(
        saved_site_reviews %>% select(Site_clean, corrected_site),
        by = "Site_clean"
      ) %>%
      mutate(
        # Apply correction using whatever exact casing you typed in the CSV
        Site = if_else(!is.na(corrected_site) & corrected_site != "", str_trim(corrected_site), Site_original)
      ) %>%
      select(-corrected_site)
    
    # Recalculate tracking names after corrections
    combined_data <- combined_data %>%
      mutate(
        Site_original = Site,
        Site_clean    = str_to_upper(str_trim(Site))
      )
    
  }
}

# Smart Fuzzy-Matching to detect potential typos among the data's own sites
unique_clean_sites <- site_summary$Site_clean
if (length(unique_clean_sites) > 1) {
  match_res <- sapply(unique_clean_sites, function(s) {
    other_sites <- setdiff(unique_clean_sites, s)
    if(length(other_sites) == 0) return(NA)
    distances <- stringdist::stringdist(s, other_sites, method = "osa")
    closest_idx <- which.min(distances)
    if(distances[closest_idx] <= 2) { 
      return(site_summary$Site_original[match(other_sites[closest_idx], site_summary$Site_clean)])
    }
    return(NA)
  })
  
  site_summary$suggested_fix <- match_res
} else {
  site_summary$suggested_fix <- NA_character_
}

if (!file.exists(site_categorised_path)) {
  site_summary %>%
    mutate(corrected_site = NA_character_) %>%
    select(Site_original, record_count, suggested_fix, corrected_site, Site_clean) %>%
    write_csv(categorised_path <- site_categorised_path) # write out template
  
  stop(
    "---- SITE REVIEW REQUIRED ----\n",
    "A site review template has been written to:\n", site_categorised_path, "\n\n",
    "Please open this CSV file:\n",
    " - Look at 'Site_original' to see your actual casing.\n",
    " - If a name has a typo, type the corrected name into 'corrected_site'.\n",
    " - If the name is fine, leave 'corrected_site' blank.\n",
    "Save the file, then re-run the script."
  )
} else {
  current_csv_data <- read_csv(site_categorised_path, show_col_types = FALSE)
  new_sites <- setdiff(site_summary$Site_clean, current_csv_data$Site_clean)
  
  if (length(new_sites) > 0) {
    updated_csv_data <- site_summary %>%
      left_join(current_csv_data %>% select(Site_clean, corrected_site), by = "Site_clean")
    
    write_csv(updated_csv_data, site_categorised_path)
    
    stop(
      "New unreviewed site names detected and added to:\n", site_categorised_path, "\n\n",
      "Please fill in any necessary 'corrected_site' fields for the new entries and re-run."
    )
  }
}

# Clean up helper columns, leaving your original 'Site' text untouched
combined_data <- combined_data %>% select(-Site_original, -Site_clean)
cat("Site names successfully reviewed and preserved in original case.\n")


## -----------------------------------------------------------------
## 3. Check BTO codes and status codes against lookups
## -----------------------------------------------------------------

## ---- 3a. Load the BTO species code lookup ----
required_species_cols <- c("Code", "Species", "BOCC5")
species_lookup_raw <- read_csv(lookup_path, show_col_types = FALSE) %>%
  rename(any_of(c(Code = "bto_code", Species = "species", BOCC5 = "bocc")))

missing_species_cols <- setdiff(required_species_cols, names(species_lookup_raw))

if (length(missing_species_cols) > 0) {
  stop(
    "BTO_species_codes_lookup.csv is missing expected column(s): ",
    paste(missing_species_cols, collapse = ", "),
    ". Columns available are:\n",
    paste(names(species_lookup_raw), collapse = ", ")
  )
}

species_lookup <- species_lookup_raw %>%
  select(Code, Species, BOCC5) %>%
  mutate(Code = str_to_upper(str_trim(Code)))

cat("\nColumns in BTO_species_codes_lookup.csv:\n")
print(names(species_lookup))

dupe_codes <- species_lookup %>%
  count(Code) %>%
  filter(n > 1)

if (nrow(dupe_codes) > 0) {
  stop(
    "BTO_species_codes_lookup.csv contains duplicate Code entries, which would ",
    "duplicate matching records when joined: ",
    paste(dupe_codes$Code, collapse = ", "),
    "\nPlease remove the duplicate row(s) from BTO_species_codes_lookup.csv and re-run."
  )
}

## ---- 3b. Load the breeding-status code lookup ----
required_status_cols <- c("Status_code", "Breeding_status")
status_lookup_raw <- read_csv(status_lookup_path, show_col_types = FALSE)

missing_status_cols <- setdiff(required_status_cols, names(status_lookup_raw))
if (length(missing_status_cols) > 0) {
  stop(
    "Species_status_codes.csv is missing expected column(s): ",
    paste(missing_status_cols, collapse = ", "),
    ". Columns available are:\n", paste(names(status_lookup_raw), collapse = ", ")
  )
}

status_lookup <- status_lookup_raw %>%
  mutate(Status_code = str_to_upper(str_trim(Status_code)))

dupe_status_codes <- status_lookup %>%
  count(Status_code) %>%
  filter(n > 1)

if (nrow(dupe_status_codes) > 0) {
  stop(
    "Species_status_codes.csv contains duplicate Status_code entries, which ",
    "would duplicate matching records when joined: ",
    paste(dupe_status_codes$Status_code, collapse = ", "),
    "\nPlease remove the duplicate row(s) from Species_status_codes.csv and re-run."
  )
}

## ---- 3c. Validate record codes ----
blank_like_status <- c("", "NA", "N/A", "NULL", "NONE", "-")

combined_data <- combined_data %>%
  mutate(
    code_valid     = BTO_code %in% species_lookup$Code,
    Status_compare = str_to_upper(str_trim(Status)),
    status_valid   = is.na(Status_compare) |
      Status_compare %in% blank_like_status |
      Status_compare %in% status_lookup$Status_code
  )

invalid_codes <- combined_data %>%
  filter(!code_valid) %>%
  st_drop_geometry() %>%
  distinct(BTO_code) %>%
  arrange(BTO_code)

invalid_statuses <- combined_data %>%
  filter(!status_valid) %>%
  st_drop_geometry() %>%
  distinct(Status) %>%
  arrange(Status)

if (nrow(invalid_statuses) > 0) {
  cat("\n---- Status code diagnostic ----\n")
  cat("Distinct Status values in the data that did NOT match the lookup:\n")
  print(invalid_statuses$Status)
  cat("\nFull list of valid Status_code values from Species_status_codes.csv:\n")
  print(sort(unique(status_lookup$Status_code)))
  cat("Compare the two lists above character-by-character.\n")
}

## -----------------------------------------------------------------
## Categorise invalid BTO codes / status codes
## -----------------------------------------------------------------
needs_review <- (nrow(invalid_codes) > 0) || (nrow(invalid_statuses) > 0)

if (needs_review) {
  
  combined_data <- combined_data %>%
    mutate(needs_review = !code_valid | !status_valid)
  
  error_summary <- combined_data %>%
    st_drop_geometry() %>%
    filter(needs_review) %>%
    mutate(
      issue_type = case_when(
        !code_valid & !status_valid ~ "BTO_code & Status",
        !code_valid                 ~ "BTO_code",
        !status_valid               ~ "Status",
        TRUE ~ NA_character_
      )
    ) %>%
    group_by(issue_type, BTO_code, Status) %>%
    summarise(
      n_records     = n(),
      example_sites = paste(unique(Site)[1:min(3, length(unique(Site)))], collapse = "; "),
      date_range    = paste(min(Date, na.rm = TRUE), "to", max(Date, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(desc(n_records))
  
  cat("\n---- BTO codes / status codes not found in the lookup files ----\n")
  print(error_summary, n = Inf)
  
  invalid_records <- combined_data %>%
    st_drop_geometry() %>%
    filter(needs_review) %>%
    select(record_id, BTO_code, code_valid, Status, status_valid, Site, Date, source_layer)
  
  if (!file.exists(categorised_path)) {
    
    invalid_records %>%
      mutate(category = NA_integer_, corrected_BTO_code = NA_character_, corrected_status = NA_character_) %>%
      write_csv(categorised_path)
    
    stop(
      "A categorisation template has been written to:\n", categorised_path,
      "\n\nOpen it and fill in the 'category' column for every row with 1, 2, 3, or 4:\n",
      "  1 - remove record from main dataset and save in separate file\n",
      "  2 - needs correcting\n",
      "  3 - delete record from file and don't save\n",
      "  4 - keep the record exactly as it is (no correction needed)\n",
      "For any row marked 2:\n",
      "  - if 'code_valid' is FALSE, fill in 'corrected_BTO_code'\n",
      "  - if 'status_valid' is FALSE, you may fill in 'corrected_status' -- or\n",
      "    leave it blank/NA to clear the status entirely (a blank status is valid)\n",
      "Save the file, then re-run this script to apply your categorisation."
    )
  }
  
  raw_categorised <- read_csv(categorised_path, show_col_types = FALSE)
  
  cat("\nColumns found in categorisation CSV:\n")
  print(names(raw_categorised))
  
  if (!"record_id" %in% names(raw_categorised) && "...1" %in% names(raw_categorised)) {
    cat("Detected blank first column heading; renaming '...1' to 'record_id'.\n")
    raw_categorised <- raw_categorised %>% rename(record_id = ...1)
  }
  
  if (!"record_id" %in% names(raw_categorised)) {
    stop(
      "The categorisation CSV does not contain a 'record_id' column.\n\n",
      "Current columns are:\n",
      paste(names(raw_categorised), collapse = ", ")
    )
  }
  
  if ("record_id" %in% names(raw_categorised)) {
    raw_categorised <- raw_categorised %>% mutate(record_id = as.character(record_id))
  }
  
  required_cols <- c("record_id", "category", "corrected_BTO_code", "corrected_status")
  missing_cols  <- setdiff(required_cols, names(raw_categorised))
  
  id_format_stale <- !"record_id" %in% missing_cols &&
    length(intersect(raw_categorised$record_id, invalid_records$record_id)) == 0
  
  if (length(missing_cols) > 0 || id_format_stale) {
    updated_template <- invalid_records %>%
      mutate(category = NA_integer_, corrected_BTO_code = NA_character_, corrected_status = NA_character_) %>%
      left_join(
        raw_categorised %>% select(any_of(c("record_id", "category", "corrected_BTO_code", "corrected_status"))),
        by = "record_id", suffix = c("", "_old")
      )
    
    if ("category_old" %in% names(updated_template)) {
      updated_template <- updated_template %>%
        mutate(category = coalesce(category_old, category)) %>%
        select(-category_old)
    }
    if ("corrected_BTO_code_old" %in% names(updated_template)) {
      updated_template <- updated_template %>%
        mutate(corrected_BTO_code = coalesce(corrected_BTO_code_old, corrected_BTO_code)) %>%
        select(-corrected_BTO_code_old)
    }
    if ("corrected_status_old" %in% names(updated_template)) {
      updated_template <- updated_template %>%
        mutate(corrected_status = coalesce(corrected_status_old, corrected_status)) %>%
        select(-corrected_status_old)
    }
    
    write_csv(updated_template, categorised_path)
    
    stop(
      "A table of all code/status errors has been saved here: ", categorised_path, "\n\n",
      "What to do next:\n",
      "  1. Open the file above.\n",
      "  2. Fill in 'category' for every row (1, 2, 3, or 4).\n",
      "  3. For any row marked 2, fill in the relevant corrected_* column(s).\n",
      "  4. Save the file.\n",
      "  5. Re-run this script."
    )
  }
  
  record_categories <- raw_categorised %>%
    select(record_id, category, corrected_BTO_code, corrected_status)
  
  dupe_record_ids <- record_categories %>%
    count(record_id) %>%
    filter(n > 1)
  
  if (nrow(dupe_record_ids) > 0) {
    stop(
      "'", categorised_path, "' contains duplicate record_id entries, ",
      "which would duplicate matching records when joined: ",
      paste(dupe_record_ids$record_id, collapse = ", "),
      "\nPlease remove the duplicate row(s) (one row per record_id) and re-run."
    )
  }
  
  missing_from_file <- setdiff(invalid_records$record_id, record_categories$record_id)
  extra_in_file      <- setdiff(record_categories$record_id, invalid_records$record_id)
  
  if (length(missing_from_file) > 0) {
    stop(
      "Categorisation CSV is missing row(s) for record_id(s):\n",
      paste(missing_from_file, collapse = ", ")
    )
  }
  
  if (length(extra_in_file) > 0) {
    warning(
      "Categorisation CSV contains extra record_id(s) that no longer need review: ",
      paste(extra_in_file, collapse = ", "),
      ". These will be ignored."
    )
  }
  
  if (any(is.na(record_categories$category))) {
    stop("Some rows in categorisation CSV still have a blank 'category'. Fill in 1, 2, 3, or 4 for all rows.")
  }
  
  if (!all(record_categories$category %in% c(1, 2, 3, 4))) {
    stop("Categorisation CSV contains a category value other than 1, 2, 3, or 4.")
  }
  
  missing_corrections <- record_categories %>%
    left_join(
      combined_data %>% st_drop_geometry() %>% select(record_id, code_valid, status_valid),
      by = "record_id"
    ) %>%
    filter(
      category == 2,
      !code_valid & (is.na(corrected_BTO_code) | str_trim(corrected_BTO_code) == "")
    )
  
  if (nrow(missing_corrections) > 0) {
    stop(
      "The following record_id(s) are marked category 2 with an invalid BTO code, but missing 'corrected_BTO_code':\n",
      paste(missing_corrections$record_id, collapse = ", ")
    )
  }
  
  ## ---- Apply decisions ----
  if ("category" %in% names(combined_data)) {
    cat("Note: a leftover 'category' column was found in combined_data and dropped before rejoining.\n")
    combined_data <- combined_data %>% select(-category)
  }
  
  combined_data <- combined_data %>%
    left_join(record_categories, by = "record_id")
  
  if (!"category" %in% names(combined_data)) {
    stop("The join did not produce a 'category' column as expected.")
  }
  
  # ---- Category 1: export and save ----
  export_data <- combined_data %>%
    filter(category == 1)
  
  if (nrow(export_data) > 0) {
    write_csv(st_drop_geometry(export_data), notes_csv_path)
    cat(
      "\n",
      nrow(export_data),
      " record(s) removed and exported to '",
      notes_csv_path,
      "'.\n",
      sep = ""
    )
    
    st_write(
      export_data,
      dsn          = gpkg_path,
      layer        = notes_layer_name,
      delete_layer = TRUE,
      quiet        = TRUE
    )
    cat("Also written to layer '", notes_layer_name, "' in '", gpkg_path, "'.\n", sep = "")
  }
  
  # ---- Category 2: apply corrections ----
  n_flagged <- sum(combined_data$category == 2, na.rm = TRUE)
  if (n_flagged > 0) {
    combined_data <- combined_data %>%
      mutate(
        BTO_code = if_else(
          !is.na(category) & category == 2 & !code_valid,
          str_to_upper(str_trim(corrected_BTO_code)),
          BTO_code
        ),
        Status = if_else(
          !is.na(category) & category == 2 & !status_valid,
          str_to_upper(str_trim(corrected_status)),
          Status
        ),
        code_valid     = BTO_code %in% species_lookup$Code,
        Status_compare = str_to_upper(str_trim(Status)),
        status_valid   = is.na(Status_compare) |
          Status_compare %in% blank_like_status |
          Status_compare %in% status_lookup$Status_code
      )
    
    still_invalid_codes <- combined_data %>%
      st_drop_geometry() %>%
      filter(category == 2, !code_valid) %>%
      distinct(record_id, BTO_code)
    
    if (nrow(still_invalid_codes) > 0) {
      warning(
        "The following corrected BTO codes still don't match lookup:\n",
        paste0("record_id ", still_invalid_codes$record_id, ": ", still_invalid_codes$BTO_code, collapse = "; ")
      )
    }
    
    still_invalid_status <- combined_data %>%
      st_drop_geometry() %>%
      filter(category == 2, !status_valid) %>%
      distinct(record_id, Status)
    
    if (nrow(still_invalid_status) > 0) {
      warning(
        "The following corrected status codes still don't match lookup:\n",
        paste0("record_id ", still_invalid_status$record_id, ": ", still_invalid_status$Status, collapse = "; ")
      )
    }
    
    cat(n_flagged, "record(s) corrected (category 2) and merged back into the main dataset.\n")
  }
  
  # ---- Category 3: discard ----
  n_discarded <- sum(combined_data$category == 3, na.rm = TRUE)
  if (n_discarded > 0) {
    cat(n_discarded, "record(s) discarded as category 3 (not saved anywhere).\n")
  }
  
  # ---- Category 4: keep unchanged ----
  n_kept_as_is <- sum(combined_data$category == 4, na.rm = TRUE)
  if (n_kept_as_is > 0) {
    cat(n_kept_as_is, "record(s) kept exactly as-is (category 4), despite an unmatched code.\n")
  }
  
  combined_data <- combined_data %>%
    mutate(manually_accepted = !is.na(category) & category == 4)
  
  combined_data <- combined_data %>%
    filter(is.na(category) | category %in% c(2, 4)) %>%
    select(-category, -corrected_BTO_code, -corrected_status, -needs_review)
  
} else {
  cat("All BTO codes and status codes matched their respective lookup files.\n")
  combined_data <- combined_data %>%
    mutate(manually_accepted = FALSE)
}

## -----------------------------------------------------------------
## Bring in human-readable lookup info
## -----------------------------------------------------------------
combined_data <- combined_data %>%
  left_join(species_lookup, by = c("BTO_code" = "Code")) %>%
  left_join(status_lookup, by = c("Status_compare" = "Status_code"))

if (!"BOCC5" %in% names(combined_data)) {
  stop("The BOCC5 column was not added to the final dataset. Check BTO_species_codes_lookup.csv.")
}

## -----------------------------------------------------------------
## Final check
## -----------------------------------------------------------------
unresolved <- combined_data %>%
  st_drop_geometry() %>%
  filter((!code_valid | !status_valid) & !manually_accepted)

if (nrow(unresolved) > 0) {
  stop(
    "Cannot proceed to save output -- ", nrow(unresolved),
    " record(s) still have unresolved invalid codes. record_id(s): ",
    paste(unresolved$record_id, collapse = ", ")
  )
}

combined_data <- combined_data %>%
  select(-Status_compare)

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

cat("Final processed data written to layer '", output_layer_name,
    "' in '", gpkg_path, "'.\n", sep = "")

# Save flat data table to project directory CSV
write_csv(st_drop_geometry(combined_data), final_csv_path)

cat(
  "Final processed data successfully saved to:\n",
  "  - Geopackage layer: '", output_layer_name, "' in '", gpkg_path, "'\n",
  "  - CSV file: '", final_csv_path, "'\n",
  sep = ""
)

