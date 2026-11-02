#updates
# need to change non-native species to black text

rm(list=ls()) # clear all previous

# ============================================================================
# WINTER BIRD SURVEY - SUMMARY TABLES
# ============================================================================
# Produces a Word document with four tables:
#   Table 1  Whole site, mean count per visit, one column
#   Table 2  Whole site, mean count per visit, one column per Site + Total
#   Table 3  Farmland records only (Field.Code A-Z), maximum count, one column per Site
#   Table 4  Woodland records only (Field.Code 1-8), maximum count, one column per Field.Code
#
# Mean    = total count / number of visits, so a visit on which a species was
#           not recorded counts as zero.
# Maximum = the highest count of the species on any one visit, as used in the
#           2015/16 and 2020/21 reports. No Total column is given for these
#           tables, as the column maxima may come from different visits.
#
# All tables are in taxonomic order (taxon_order in the species lookup).
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(flextable)
library(officer)
library(lubridate)
library(sf)

# ----------------------------------------------------------------------------
# 0a. INPUT DATA - EDIT THESE, THEN CHECK THE DATES PRINTED IN THE CONSOLE
# ----------------------------------------------------------------------------

INPUT <- list(
  input_path = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg',
  layer_name = '2026_Counts_final'
)

dt_raw <- st_drop_geometry(st_read(dsn = INPUT$input_path, layer = INPUT$layer_name, quiet = TRUE))

# Find a column by any of its possible names
find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

date_col  <- find_col(dt_raw, c("Date", "date"))
field_col <- find_col(dt_raw, c("Field.Code", "Field_Code", "Field"))

if (is.na(date_col))  stop("No 'Date' column found in layer ", INPUT$layer_name)
if (is.na(field_col)) stop("No 'Field.Code' column found in layer ", INPUT$layer_name)
if (!"Site" %in% names(dt_raw)) stop("No 'Site' column found in layer ", INPUT$layer_name)

dt <- dt_raw %>%
  mutate(
    Count          = coalesce(as.numeric(Count), 0),
    Species        = str_squish(as.character(Species)),
    Site           = str_to_title(str_squish(coalesce(as.character(Site), "Unspecified"))),
    Field_Code     = toupper(str_squish(coalesce(as.character(.data[[field_col]]), ""))),
    Date_parsed    = parse_date_time(.data[[date_col]], orders = c("dmy", "ymd", "mdy")),
    Date_formatted = format(Date_parsed, "%d/%m/%Y")
  ) %>%
  filter(!is.na(Species), Species != "", !is.na(Date_parsed))

message("\n--- AVAILABLE SURVEY DATES IN DATASET ---")
print(as.data.frame(
  dt %>% count(Date_parsed, Survey_Date = Date_formatted, name = "Record_Count") %>%
    arrange(Date_parsed) %>% select(Survey_Date, Record_Count)
))
message("-----------------------------------------\n")

# ----------------------------------------------------------------------------
# 0b. CONFIG - EDIT THESE (use the dates printed above for dates_to_include)
# ----------------------------------------------------------------------------

CONFIG <- list(
  # Species lookup with columns: species, taxon_order, bocc (Red / Amber / Green)
  species_lookup_path   = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv',
  output_dir            = '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results/Winter',
  
  # Visits to include, e.g. c("14/01/2026", "20/02/2026"), or NULL for all dates
  dates_to_include      = c("2025-12-20", "2026-02-09"),
  
  # Field.Code values used for Table 3 (farmland) and Table 4 (woodland)
  farmland_codes        = LETTERS,
  woodland_codes        = as.character(1:8),
  
  bocc5_order           = c("Red", "Amber", "Green"),
  bocc_text_colors      = c("Red" = "#C00000", "Amber" = "#D97706", "Green" = "#15803D")
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1. DATE FILTER
# ----------------------------------------------------------------------------

if (!is.null(CONFIG$dates_to_include) && length(CONFIG$dates_to_include) > 0) {
  wanted <- format(parse_date_time(CONFIG$dates_to_include, orders = c("dmy", "ymd", "mdy")), "%d/%m/%Y")
  if (any(is.na(wanted))) stop("Could not read some of CONFIG$dates_to_include as dates")
  missing_dates <- setdiff(wanted, dt$Date_formatted)
  if (length(missing_dates) > 0) warning("No records for: ", paste(missing_dates, collapse = ", "))
  dt <- dt %>% filter(Date_formatted %in% wanted)
}
if (nrow(dt) == 0) stop("No records left after the date filter")

# ----------------------------------------------------------------------------
# 2. SPECIES ORDER (TAXONOMIC) AND BOCC5 STATUS
# ----------------------------------------------------------------------------

if (!file.exists(CONFIG$species_lookup_path)) stop("Species lookup not found: ", CONFIG$species_lookup_path)

species_lookup <- read.csv(CONFIG$species_lookup_path, stringsAsFactors = FALSE) %>%
  transmute(
    Species     = str_squish(as.character(species)),
    taxon_order = suppressWarnings(as.numeric(taxon_order)),
    BOCC5       = str_to_title(str_squish(as.character(bocc)))
  ) %>%
  distinct(Species, .keep_all = TRUE)

if (all(is.na(species_lookup$taxon_order))) {
  stop("The taxon_order column in the species lookup is empty or not numeric: ", CONFIG$species_lookup_path)
}

species_order <- dt %>%
  distinct(Species) %>%
  left_join(species_lookup, by = c("Species")) %>%
  mutate(
    BOCC5       = factor(BOCC5, levels = CONFIG$bocc5_order),
    taxon_order = coalesce(taxon_order, 9999)
  )

# Species in the data that are not in the lookup (or have no taxon_order):
# listed last, uncoloured
no_order <- species_order %>% filter(taxon_order == 9999) %>% pull(Species)
if (length(no_order) > 0) {
  warning("Species with no taxon_order in the lookup (listed at the end, no BOCC colour): ",
          paste(sort(no_order), collapse = ", "),
          "\nCheck the spelling against the 'species' column of the lookup.")
}

# ----------------------------------------------------------------------------
# 3. TABLE BUILDER
# ----------------------------------------------------------------------------

# At most 1 decimal place, no trailing ".0" (3 -> "3", 2.5 -> "2.5").
# A recorded species whose value rounds to zero is shown as "<0.1".
format_value <- function(x) {
  r   <- round(x, 1)
  out <- sub("\\.0$", "", formatC(r, format = "f", digits = 1))
  out[x > 0 & r == 0] <- "<0.1"
  out
}

# df        : the records to summarise (already filtered)
# split_col : NULL for a single column, or "Site" / "Field_Code"
#             for one column per value
# stat      : "mean" = total count / number of visits (with a Total column
#                      when split_col is given)
#             "max"  = highest count on any one visit (no Total column)
build_summary_table <- function(df, split_col = NULL, caption, stat = c("mean", "max")) {
  stat <- match.arg(stat)
  if (nrow(df) == 0) {
    warning("No records for: ", caption)
    return(NULL)
  }
  
  n_visits   <- n_distinct(df$Date_formatted)
  value_name <- if (stat == "max") "Max Count" else "Mean Count"
  summarise_visits <- function(x) if (stat == "max") max(x) else sum(x) / n_visits
  
  # Count per species per visit (and per column)
  per_visit <- df %>%
    mutate(Column = if (is.null(split_col)) value_name else as.character(.data[[split_col]])) %>%
    group_by(Species, Date_formatted, Column) %>%
    summarise(n = sum(Count), .groups = "drop")
  
  tab <- per_visit %>%
    group_by(Species, Column) %>%
    summarise(value = summarise_visits(n), .groups = "drop") %>%
    pivot_wider(names_from = Column, values_from = value, values_fill = 0)
  
  # Numeric columns in number order (1, 2, 10), then names A-Z
  data_cols <- setdiff(names(tab), "Species")
  data_cols <- data_cols[order(suppressWarnings(as.numeric(data_cols)), data_cols)]
  
  # Total column for mean tables only (sum of the unrounded means)
  if (!is.null(split_col) && stat == "mean") {
    tab$Total <- rowSums(tab[data_cols])
    data_cols <- c(data_cols, "Total")
  }
  
  tab <- tab %>%
    filter(if_any(all_of(data_cols), ~ . > 0)) %>%
    left_join(species_order, by = "Species") %>%
    arrange(taxon_order, Species) %>%
    mutate(across(all_of(data_cols), format_value))
  if (nrow(tab) == 0) return(NULL)
  
  ft <- flextable(tab, col_keys = c("Species", data_cols)) %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Calibri", part = "all") %>%
    bold(part = "header") %>%
    flextable::align(j = "Species", align = "left", part = "all") %>%
    flextable::align(j = data_cols, align = "center", part = "all") %>%
    hline_top(border = fp_border(width = 2), part = "header") %>%
    hline_bottom(border = fp_border(width = 1), part = "header") %>%
    padding(padding = 3, part = "all") %>%
    set_caption(caption = caption) %>%
    set_table_properties(width = 1, layout = "autofit")
  
  for (status in names(CONFIG$bocc_text_colors)) {
    rows <- which(tab$BOCC5 == status)
    if (length(rows) > 0) ft <- color(ft, i = rows, color = CONFIG$bocc_text_colors[[status]], part = "body")
  }
  
  if ("Total" %in% data_cols) {
    ft <- ft %>%
      bold(j = "Total") %>%
      bg(j = "Total", bg = "grey85", part = "body") %>%
      color(j = "Total", color = "black", part = "body")
  }
  
  ft
}

# ----------------------------------------------------------------------------
# 4. BUILD THE FOUR TABLES AND WRITE THE WORD DOCUMENT
# ----------------------------------------------------------------------------

dt_farmland <- dt %>% filter(Field_Code %in% toupper(CONFIG$farmland_codes))
dt_woodland <- dt %>% filter(Field_Code %in% toupper(CONFIG$woodland_codes))

# Records that are in neither Table 3 nor Table 4 (they are still in Tables 1 and 2)
unassigned <- dt %>% filter(!Field_Code %in% toupper(c(CONFIG$farmland_codes, CONFIG$woodland_codes)))
if (nrow(unassigned) > 0) {
  message("Records with a Field.Code outside the farmland and woodland lists (Tables 1 and 2 only):")
  print(as.data.frame(unassigned %>% count(Field_Code, name = "Record_Count")))
}

tables <- list(
  list(title = "Table 1. Mean count per species across the whole site",
       ft    = build_summary_table(dt, NULL, "Mean count per species across the whole site", stat = "mean")),
  list(title = "Table 2. Mean count per species by site",
       ft    = build_summary_table(dt, "Site", "Mean count per species by site", stat = "mean")),
  list(title = "Table 3. Maximum count per species by site (farmland)",
       ft    = build_summary_table(dt_farmland, "Site", "Maximum count per species by site (farmland)", stat = "max")),
  list(title = "Table 4. Maximum count per species by woodland field code",
       ft    = build_summary_table(dt_woodland, "Field_Code", "Maximum count per species by woodland field code", stat = "max"))
)

doc   <- read_docx()
first <- TRUE
for (tb in tables) {
  if (is.null(tb$ft)) next
  if (!first) doc <- body_add_break(doc)
  doc   <- doc %>% body_add_par(tb$title, style = "heading 1") %>% body_add_flextable(tb$ft)
  first <- FALSE
}

out_file <- file.path(CONFIG$output_dir, paste0("winter_survey_summary_tables_", INPUT$layer_name, ".docx"))
print(doc, target = out_file)
message("Saved summary tables to: ", out_file)