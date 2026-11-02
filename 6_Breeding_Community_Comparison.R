rm(list=ls())

# =============================================================================
# Bird Survey Community Analysis — Temporal effort-matching version (v4)
# Using BBS_bird_dataset.csv (Zenodo DOI: 10.5281/zenodo.18598453)
# Massimino et al. 2025, Global Ecology & Biogeography
#
# Outputs produced:
#   Figures: site RAI, species abundance vs BBS, species divergence (regional
#            and national), species overlap, rarefied richness, occupancy
#   Tables:  site RAI, abundance comparison, species overlap, rarefied
#            richness, absent species, site-only species, occupancy, summary
#   Report:  narrative Word summary document
# =============================================================================

library(tidyverse)
library(vegan)
library(scales)
library(data.table)
library(officer)
library(flextable)
library(sf)
library(ggrepel)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# --- File paths ---
BBS_BIRD_FILE  <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/BBS_bird_dataset.csv"
BBS_VISIT_FILE <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/BBS_visit_dataset.csv"
SPECIES_LOOKUP <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv"
GRID_COORDS    <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/grid_square_coordinates_lookup.csv"
# --- Survey data sources ---
# Use the GeoPackage, the CSV, or both (rows are combined). Set either to NA to skip it.
input_path     <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'
layer_name     <- '2026_Counts_final'
csv_path       <- NA   # e.g. '/Users/david/.../2026_Counts_extra.csv'
# Column names in the CSV (change the right-hand side to match your file).
# BOCC5 is optional.
CSV_COLS <- c(species_code = "BTO_code",
              transect     = "Site",
              count        = "Count",
              date_str     = "Date")
# --- Site details ---
SITE_NAME      <- "Shimpling Park Farm"   # appears in all figure titles and report
SITE_GRID_REF  <- "TL883512"              # BTO/OS grid reference

# --- Date selection & Season settings ---
# Set start/end dates (YYYY-MM-DD) or set to NA to include all rows in CSV
START_DATE        <- "2026-04-01"   
END_DATE          <- "2026-06-30"   
SEASON_NAME       <- "Spring 2026"
OVERRIDE_BBS_YEAR <- 2024           # Set target BBS reference year, or NA to auto-detect max year

# --- BBS cache ---
# Cache file saves pre-processed max(E, L) square counts to speed up future runs.
BBS_CACHE_FILE <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/bbs_cache_max.rds'

# --- BBS processing ---
INCLUDE_DISTANCE_BAND_3 <- TRUE   # TRUE = all bands (1, 2, 3 and F); cache rebuilds if changed

# --- Species lookups for the client-facing figures (sections 15a/15b) ---
# BOCC5 status comes from the BTO species codes lookup (columns bto_code, bocc).
SPECIES_CODES_FILE <- SPECIES_LOOKUP   # same BTO lookup supplies names and BOCC5

# --- BBS species code harmonisation ---
# The BBS dataset's species list differs from the BTO lookup in places.
# Map BBS codes to BTO codes (BBS code = BTO code), and drop codes that are
# not wild birds so they don't inflate BBS richness benchmarks.
BBS_CODE_MAP      <- c(RD = "LR")                 # BBS "Redpoll" (RD) = BTO code LR, named "Redpoll" in outputs
BBS_EXCLUDE_CODES <- c("PX", "ZF", "RV")          # Indian Peafowl, Domestic Mallard, Reeves's Pheasant
# Sensitive species were withheld from the published BBS dataset (see
# list_of_sensitive_species.csv in the BBS download). Site species with no BBS
# record anywhere are left out of the count-based comparisons (abundance,
# richness, percentiles, composition), because BBS has no counts for them.
# Codes listed here are additionally ASSUMED PRESENT in the regional and
# national pools, so they count as shared species in the pool overlap
# (Jaccard/Sorensen, Fig 7, Table 3) instead of appearing "site-only".
ASSUME_PRESENT_IN_BBS <- c("KT")   # Red Kite
SENSITIVE_FILE <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/list_of_sensitive_species.csv"  # NA to skip
# Guilds are YOURS to edit (columns species, bto_code, guild, farmland_index):
#   guild          - any group name you like (Farmland, Woodland, Wetland, ...).
#                    Each name becomes a section of the heat strip (Fig 11).
#                    Blank = "Other".
#   farmland_index - "Specialist" or "Generalist" puts the species in the
#                    farmland bird score (Fig 13); blank = not included.
#   diet, nest_site, migration - not currently used (kept for future figures).
GUILD_FILE  <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/species_guilds.csv"
GUILD_ORDER <- c("Farmland", "Woodland", "Wetland", "Generalist")  # display order; other names follow, "Other" last

# --- Effort control: TEMPORAL effort matching ---
# Each BBS square contributes one 2 km transect per year. Rather than pooling
# several squares (which mixes habitats and inflates richness), effort is
# matched by pooling the SAME square across successive years. A 9.2 km site
# transect is therefore compared against "pseudo-transects" made of 4.6 years
# of data from a single BBS square: the 4 most recent surveyed years at full
# weight plus the 5th most recent year at 0.6 weight.
SITE_TRANSECT_KM <- 9.2   # your transect length in km
BBS_TRANSECT_KM  <- 2.0   # BBS standard — do not change
N_YEARS_EQUIV    <- SITE_TRANSECT_KM / BBS_TRANSECT_KM   # BBS square-years per pseudo-transect
K_FULL           <- floor(N_YEARS_EQUIV)                 # years at full weight
K_FRAC           <- N_YEARS_EQUIV - K_FULL               # weight of the last, partial year
K_NEEDED         <- ceiling(N_YEARS_EQUIV)               # surveyed years a square must have
LOOKBACK_MAX_YEARS <- 10  # only look this far back from the BBS reference year
# (years a square was not surveyed are skipped, e.g. 2020)

# --- Adaptive regional radius ---
RADIUS_MIN_KM  <- 30
RADIUS_MAX_KM  <- 80
RADIUS_STEP_KM <- 5
MIN_SQUARES    <- 50      # minimum number of ELIGIBLE squares (>= K_NEEDED surveyed years)

# --- Output directory ---
out_dir <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results/Summer/Community"
dir.create(out_dir, showWarnings = FALSE)

# --- Helpers ---
# Figures are saved WITHOUT titles and subtitles when SHOW_FIGURE_TITLES is
# FALSE. The title and subtitle text of every figure is recorded instead and
# written to figure_captions.docx / figure_captions.csv in the output folder,
# ready to paste as captions.
SHOW_FIGURE_TITLES <- FALSE
figure_captions <- tibble(file = character(), title = character(), notes = character())

save_fig <- function(plot, filename, width = 18, height = 12, dpi = 300) {
  lb <- tryCatch(ggplot2::get_labs(plot), error = function(e) plot$labels)
  one_line <- function(x, sep = " ") {
    if (is.null(x) || !is.character(x)) return("")
    str_squish(gsub("\n", sep, x))
  }
  figure_captions <<- bind_rows(
    figure_captions[figure_captions$file != filename, ],
    tibble(file = filename, title = one_line(lb$title, " \u2014 "), notes = one_line(lb$subtitle)))
  if (!SHOW_FIGURE_TITLES) plot <- plot + labs(title = NULL, subtitle = NULL)
  ggsave(file.path(out_dir, filename), plot = plot,
         width = width, height = height, units = "cm", dpi = dpi, bg = "white")
}

# Borderless labels that work in old and new ggplot2 (label.size was replaced by linewidth)
.new_label <- "border.colour" %in% names(formals(ggplot2::geom_label))
geom_label_nb <- function(...) if (.new_label) geom_label(..., linewidth = 0) else geom_label(..., label.size = 0)

save_word_table <- function(ft, filename, title_text = NULL, note_text = NULL) {
  doc <- read_docx()
  if (!is.null(title_text))
    doc <- body_add_par(doc, title_text, style = "heading 2")
  doc <- body_add_flextable(doc, ft)
  if (!is.null(note_text))
    doc <- body_add_par(doc, note_text, style = "Normal")
  print(doc, target = file.path(out_dir, filename))
  message("Saved: ", filename)
}

effort_note <- if (K_FRAC == 0) {
  sprintf("BBS benchmarks built from %d-year pseudo-transects: %d surveyed years at a single BBS square (%d × 2 km = %.1f km, matching the site transect). No spatial pooling or scaling.",
          K_FULL, K_FULL, K_FULL, SITE_TRANSECT_KM)
} else {
  sprintf("BBS benchmarks built from %.1f-year pseudo-transects: the %d most recent surveyed years at a single BBS square plus the next year at %.0f%% weight (%.1f × 2 km = %.1f km, matching the site transect). No spatial pooling or scaling.",
          N_YEARS_EQUIV, K_FULL, K_FRAC * 100, N_YEARS_EQUIV, SITE_TRANSECT_KM)
}
message("Effort control: ", effort_note)

# =============================================================================
# 1. LOAD AND FILTER SURVEY DATA
# =============================================================================

# Dates may arrive as Date objects, "YYYY-MM-DD" or UK "DD/MM/YYYY" (Excel)
parse_survey_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(parse_date_time(as.character(x), orders = c("Ymd", "dmY", "dmy", "Ymd HMS", "dmY HM"), quiet = TRUE))
}

survey_sources <- list()

if (!is.na(input_path)) {
  survey_sources$gpkg <- st_read(input_path, layer = layer_name, quiet = TRUE) %>%
    st_drop_geometry() %>%
    as_tibble() %>%
    rename(species_code = BTO_code,
           transect     = Site,
           count        = Count,
           date_str     = Date) %>%
    mutate(source = "gpkg")
  message(sprintf("GeoPackage: %d rows from layer '%s'.", nrow(survey_sources$gpkg), layer_name))
}

if (!is.na(csv_path)) {
  csv_raw <- read_csv(csv_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  missing_cols <- setdiff(CSV_COLS, names(csv_raw))
  if (length(missing_cols) > 0)
    stop("CSV is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Columns found: ", paste(names(csv_raw), collapse = ", "),
         ". Update CSV_COLS in the configuration.")
  survey_sources$csv <- csv_raw %>%
    rename(all_of(CSV_COLS)) %>%
    mutate(count = as.numeric(count), source = "csv")
  if (any(is.na(survey_sources$csv$count)))
    warning(sum(is.na(survey_sources$csv$count)), " CSV rows have a non-numeric count and will be treated as 0.")
  message(sprintf("CSV: %d rows from %s.", nrow(survey_sources$csv), basename(csv_path)))
}

if (length(survey_sources) == 0) stop("No survey data: set input_path and/or csv_path.")

data_raw <- survey_sources %>%
  map(~ mutate(.x, date_parsed = parse_survey_date(date_str),
               across(any_of(c("species_code", "transect", "BOCC5")), as.character))) %>%
  bind_rows() %>%
  mutate(species_code = str_trim(species_code),
         count        = replace_na(as.numeric(count), 0),
         visit        = format(date_parsed, "%Y-%m-%d"),   # visit is defined by the date
         year         = year(date_parsed))

# Flag exact duplicates between sources (same species, transect, date, count)
if (length(survey_sources) > 1) {
  dups <- data_raw %>%
    count(species_code, transect, visit, count, name = "n_rows") %>%
    filter(n_rows > 1)
  if (nrow(dups) > 0)
    warning(nrow(dups), " species/transect/date/count combinations appear more than once across ",
            "the GeoPackage and CSV - check the CSV isn't repeating records already in the GeoPackage.")
}

bad_dates <- data_raw %>% filter(is.na(date_parsed))
if (nrow(bad_dates) > 0) {
  warning(nrow(bad_dates), " rows have dates that could not be read and are excluded (e.g. '",
          paste(head(unique(bad_dates$date_str), 3), collapse = "', '"), "' from ",
          paste(unique(bad_dates$source), collapse = "/"), ").")
  data_raw <- data_raw %>% filter(!is.na(date_parsed))
}

# Filter by custom date range if defined
if (!is.na(START_DATE) && !is.na(END_DATE)) {
  data <- data_raw %>%
    filter(date_parsed >= as.Date(START_DATE) & date_parsed <= as.Date(END_DATE))
  message(sprintf("Filtered survey data from %s to %s. Detections remaining: %d",
                  START_DATE, END_DATE, nrow(data)))
} else {
  data <- data_raw
}

# Determine reference year for BBS matching
if (!is.na(OVERRIDE_BBS_YEAR)) {
  year_site <- OVERRIDE_BBS_YEAR
} else {
  year_site <- max(data$year, na.rm = TRUE)
}

message(sprintf("Analysis Period: %s | BBS Benchmark Target Year: %d | Site: %s (%s)",
                SEASON_NAME, year_site, SITE_NAME, SITE_GRID_REF))

# Aggregate counts across sub-sites for each date/visit first
data_aggregated <- data %>%
  group_by(species_code, visit) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

# Zero-fill: all species x visit combinations
all_combos <- expand_grid(
  species_code = unique(data_aggregated$species_code),
  visit        = unique(data_aggregated$visit)
)

data_full <- all_combos %>%
  left_join(data_aggregated, by = c("species_code", "visit")) %>%
  replace_na(list(count = 0)) %>%
  mutate(
    date_parsed = as.Date(visit),
    year        = year(date_parsed)
  )

# =============================================================================
# 2. SITE RAI (MAX COUNT ACROSS VISITS)
# =============================================================================

RAI_site <- data_full %>%
  group_by(species_code) %>%
  summarise(
    n_visits    = n(),
    total_count = sum(count),
    RAI         = max(count),  # Updated: Max count across visits per species
    .groups     = "drop"
  )

site_wide <- RAI_site %>%
  mutate(rel_abund = RAI / sum(RAI)) %>%
  select(species_code, rel_abund)

site_present       <- RAI_site %>% filter(total_count > 0) %>% pull(species_code)
site_total_indivs  <- sum(data_full$count)
site_total_species <- length(site_present)

message(sprintf("Site: %d species, %d total detections across %d visits (using max count per species).",
                site_total_species, site_total_indivs, n_distinct(data_full$visit)))

# =============================================================================
# 3. BBS YEAR MAPPING
# =============================================================================

BBS_MAX_YEAR <- 2024
BBS_MIN_YEAR <- 1994

year_bbs <- if (year_site > BBS_MAX_YEAR) {
  message(sprintf("Survey reference year %d is beyond BBS dataset range. Using %d.",
                  year_site, BBS_MAX_YEAR))
  BBS_MAX_YEAR
} else if (year_site < BBS_MIN_YEAR) {
  message(sprintf("Survey reference year %d is before BBS dataset range. Using %d.",
                  year_site, BBS_MIN_YEAR))
  BBS_MIN_YEAR
} else {
  year_site
}

message(sprintf("BBS reference year selected: %d", year_bbs))

# =============================================================================
# 4. LOAD OR RESTORE BBS DATA (WITH MAX E/L VISIT LOGIC)
# =============================================================================

# The cache records the settings it was built with; if they differ from the
# current settings (e.g. band 3 switched on/off) it is rebuilt automatically.
cache_settings <- list(include_band_3 = INCLUDE_DISTANCE_BAND_3, version = 5)

use_cache <- FALSE
if (!is.na(BBS_CACHE_FILE) && file.exists(BBS_CACHE_FILE)) {
  bbs_cache <- readRDS(BBS_CACHE_FILE)
  use_cache <- identical(bbs_cache$settings, cache_settings)
  if (!use_cache) message("Cache was built with different settings \u2014 rebuilding.")
}

if (use_cache) {
  
  message("Loading processed BBS data from cache...")
  bbs_by_square  <- bbs_cache$bbs_by_square
  bbs_visits     <- bbs_cache$bbs_visits
  grid_coords    <- bbs_cache$grid_coords
  message("Cache loaded.")
  
} else {
  
  message("Loading BBS dataset and collapsing visits to max(E, L)...")
  bbs_raw     <- fread(BBS_BIRD_FILE)
  bbs_visits  <- fread(BBS_VISIT_FILE)
  grid_coords <- fread(GRID_COORDS) %>% rename_with(tolower)
  message("BBS bird file columns: ", paste(names(bbs_raw), collapse = ", "))
  
  # --- Distance bands ---
  if ("distance_band" %in% names(bbs_raw)) {
    message("Distance bands present: ", paste(sort(unique(bbs_raw$distance_band)), collapse = ", "))
    if (!INCLUDE_DISTANCE_BAND_3) {
      bbs_raw <- bbs_raw[distance_band != 3]
      message("Excluded distance_band == 3 from BBS dataset.")
    } else {
      message("All distance bands included (band 3 retained).")
    }
  }
  
  section_cols <- paste0("section_", 1:10)
  bbs_raw[, row_total := rowSums(.SD, na.rm = TRUE), .SDcols = section_cols]
  
  # --- Collapse to one count per square x year x species ---
  # Step 1: SUM across distance bands (and any other split rows) within each
  #         visit, so each visit gets its full count.
  # Step 2: take the MAX across the early and late visits.
  # (Taking max() directly over rows would pick the largest single distance
  #  band when the file has one row per visit x band.)
  visit_col <- intersect(c("visit", "visit_number", "visit_period", "early_late", "visit_id"),
                         names(bbs_raw))[1]
  if (!is.na(visit_col)) {
    message("Visit column detected: '", visit_col, "'. Summing bands within visits, then max(E, L).")
    per_visit <- bbs_raw[, .(visit_total = sum(row_total, na.rm = TRUE)),
                         by = c("year", "square", visit_col, "species_code")]
    bbs_by_square <- per_visit[, .(square_count = max(visit_total, na.rm = TRUE)),
                               by = .(year, square, species_code)]
  } else {
    warning("No visit column found in the BBS bird file \u2014 taking max over rows. ",
            "If rows are split by distance band this undercounts; check the column list above.")
    bbs_by_square <- bbs_raw[, .(square_count = max(row_total, na.rm = TRUE)),
                             by = .(year, square, species_code)]
  }
  
  if (!is.na(BBS_CACHE_FILE)) {
    message("Saving processed BBS data to cache for future runs...")
    saveRDS(list(bbs_by_square = bbs_by_square,
                 bbs_visits    = bbs_visits,
                 grid_coords   = grid_coords,
                 settings      = cache_settings),
            BBS_CACHE_FILE)
    message("Cache saved: ", BBS_CACHE_FILE)
  }
}

# Species names are read fresh each run (not cached) from the BTO lookup
sp_lookup <- read_csv(SPECIES_LOOKUP, col_types = cols(.default = "c")) %>%
  transmute(species_code = str_trim(bto_code), english_name = str_trim(species)) %>%
  filter(!is.na(species_code)) %>%
  distinct(species_code, .keep_all = TRUE)
sp_names <- sp_lookup %>% select(species_code, english_name)

# --- 4b. Harmonise BBS codes with the BTO lookup ---
bbs_by_square <- as.data.table(bbs_by_square)
if (length(BBS_CODE_MAP) > 0) {
  n_mapped <- bbs_by_square[species_code %in% names(BBS_CODE_MAP), .N]
  bbs_by_square[species_code %in% names(BBS_CODE_MAP), species_code := unname(BBS_CODE_MAP[species_code])]
  message(sprintf("BBS codes remapped (%s): %d rows.",
                  paste(names(BBS_CODE_MAP), BBS_CODE_MAP, sep = " -> ", collapse = ", "), n_mapped))
}
if (length(BBS_EXCLUDE_CODES) > 0) {
  n_excl <- bbs_by_square[species_code %in% BBS_EXCLUDE_CODES, .N]
  bbs_by_square <- bbs_by_square[!species_code %in% BBS_EXCLUDE_CODES]
  message(sprintf("BBS codes excluded (%s): %d rows.", paste(BBS_EXCLUDE_CODES, collapse = ", "), n_excl))
}
# Re-collapse in case a mapping merged two codes into one
bbs_by_square <- bbs_by_square[, .(square_count = sum(square_count)), by = .(year, square, species_code)]

# Any remaining BBS codes the BTO lookup doesn't know (candidates for BBS_CODE_MAP / BBS_EXCLUDE_CODES)
bbs_unknown <- bbs_by_square[!species_code %in% sp_lookup$species_code,
                             .(squares = uniqueN(square), birds = sum(square_count)), by = species_code][order(-squares)]
if (nrow(bbs_unknown) > 0)
  message("BBS codes not in the BTO lookup (add to BBS_CODE_MAP or BBS_EXCLUDE_CODES): ",
          paste(sprintf("%s (%d squares)", bbs_unknown$species_code, bbs_unknown$squares), collapse = ", "))

# --- 4c. Site species that the BBS dataset does not cover at all ---
# These can't be benchmarked, so they are kept for the site-only figures
# (Fig 1, Table 1) but removed from every site-vs-BBS comparison.
RAI_site_all    <- RAI_site
site_not_in_bbs <- setdiff(site_present, unique(bbs_by_square$species_code))
if (length(site_not_in_bbs) > 0) {
  message("Site species not in the BBS dataset (excluded from comparisons): ",
          paste(site_not_in_bbs, collapse = ", "))
  RAI_site           <- RAI_site %>% filter(!species_code %in% site_not_in_bbs)
  site_wide          <- RAI_site %>% mutate(rel_abund = RAI / sum(RAI)) %>% select(species_code, rel_abund)
  site_present       <- setdiff(site_present, site_not_in_bbs)
  site_total_species <- length(site_present)
  site_total_indivs  <- sum(data_full$count[!data_full$species_code %in% site_not_in_bbs])
}
assumed_present <- intersect(site_not_in_bbs, ASSUME_PRESENT_IN_BBS)
if (length(assumed_present) > 0)
  message("Assumed present regionally and nationally (withheld from BBS): ", paste(assumed_present, collapse = ", "))

# Sensitive species recorded on the farm (matched to the BTO lookup by
# scientific name, or by English name for recent taxonomic changes)
site_sensitive <- character(0)
if (!is.na(SENSITIVE_FILE) && file.exists(SENSITIVE_FILE)) {
  sens <- fread(file = SENSITIVE_FILE)
  bto_full <- fread(file = SPECIES_LOOKUP)
  strip_q <- function(x) str_trim(gsub("^(western|eurasian|european|common|northern|pied|red-billed) ", "", tolower(x)))
  sens_codes <- bto_full[sci_name %in% sens$scientific_name |
                           strip_q(species) %in% strip_q(sens$English_name),
                         unique(bto_code)]
  site_sensitive <- intersect(RAI_site_all$species_code[RAI_site_all$total_count > 0], sens_codes)
  if (length(site_sensitive) > 0) {
    in_bbs <- site_sensitive %in% unique(bbs_by_square$species_code)
    message("BBS sensitive species recorded on the farm: ",
            paste0(site_sensitive, ifelse(in_bbs, " (in BBS data)", " (withheld from BBS)"), collapse = ", "))
  }
}

# =============================================================================
# 5. GRID REFERENCE PARSER AND ADAPTIVE RADIUS
# =============================================================================

parse_gridref <- function(gr) {
  gr <- toupper(trimws(gr))
  letters_500 <- matrix(c(
    "SV","SW","SX","SY","SZ","TV","TW",
    "SQ","SR","SS","ST","SU","TQ","TR",
    "SL","SM","SN","SO","SP","TL","TM",
    "SF","SG","SH","SJ","SK","TF","TG",
    "SA","SB","SC","SD","SE","TA","TB",
    "NV","NW","NX","NY","NZ","OV","OW",
    "NQ","NR","NS","NT","NU","OQ","OR",
    "NL","NM","NN","NO","NP","OL","OM",
    "NF","NG","NH","NJ","NK","OF","OG",
    "NA","NB","NC","ND","NE","OA","OB",
    "HV","HW","HX","HY","HZ","JV","JW",
    "HQ","HR","HS","HT","HU","JQ","JR",
    "HL","HM","HN","HO","HP","JL","JM"
  ), nrow = 13, ncol = 7, byrow = TRUE)
  prefix  <- substr(gr, 1, 2)
  digits  <- substr(gr, 3, nchar(gr))
  n_dig   <- nchar(digits)
  if (n_dig %% 2 != 0 || n_dig == 0) stop("Grid reference must have even digits: ", gr)
  half    <- n_dig / 2
  pos     <- which(letters_500 == prefix, arr.ind = TRUE)
  if (nrow(pos) == 0) stop("Unrecognised prefix: ", prefix)
  e100    <- (pos[1,2] - 1) * 100000
  n100    <- (pos[1,1] - 1) * 100000
  scale   <- 10^(5 - half)
  list(easting  = e100 + as.numeric(substr(digits, 1, half)) * scale,
       northing = n100 + as.numeric(substr(digits, half+1, n_dig)) * scale,
       precision_m = scale)
}

osgb_to_latlon <- function(E, N) {
  a <- 6377563.396; b <- 6356256.909; F0 <- 0.9996012717
  lat0 <- 49*pi/180; lon0 <- -2*pi/180; N0 <- -100000; E0 <- 400000
  e2 <- 1-(b^2/a^2); n <- (a-b)/(a+b)
  lat <- lat0; M <- 0
  repeat {
    lat <- (N-N0-M)/(a*F0)+lat
    Ma  <- (1+n+(5/4)*n^2+(5/4)*n^3)*(lat-lat0)
    Mb  <- (3*n+3*n^2+(21/8)*n^3)*sin(lat-lat0)*cos(lat+lat0)
    Mc  <- ((15/8)*n^2+(15/8)*n^3)*sin(2*(lat-lat0))*cos(2*(lat+lat0))
    Md  <- (35/24)*n^3*sin(3*(lat-lat0))*cos(3*(lat+lat0))
    M   <- b*F0*(Ma-Mb+Mc-Md)
    if (abs(N-N0-M) < 0.00001) break
  }
  nu   <- a*F0/sqrt(1-e2*sin(lat)^2)
  rho  <- a*F0*(1-e2)/(1-e2*sin(lat)^2)^1.5
  tL   <- tan(lat); sL <- 1/cos(lat); dE <- E-E0
  lat_out <- lat-(tL/(2*rho*nu))*dE^2+(tL/(24*rho*nu^3))*(5+3*tL^2+nu/rho-1-9*tL^2*(nu/rho-1))*dE^4
  lon_out <- lon0+sL/nu*dE-sL/(6*nu^3)*(nu/rho+2*tL^2)*dE^3+sL/(120*nu^5)*(5+28*tL^2+24*tL^4)*dE^5
  list(lat = lat_out*180/pi, lon = lon_out*180/pi)
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  dlat <- (lat2-lat1)*pi/180; dlon <- (lon2-lon1)*pi/180
  a <- sin(dlat/2)^2 + cos(lat1*pi/180)*cos(lat2*pi/180)*sin(dlon/2)^2
  2*R*asin(sqrt(a))
}

site_loc <- parse_gridref(SITE_GRID_REF)
site_ll  <- osgb_to_latlon(site_loc$easting, site_loc$northing)
message(sprintf("Site lat/lon: %.5f, %.5f", site_ll$lat, site_ll$lon))

grid_coords <- grid_coords %>%
  mutate(dist_km = haversine_km(site_ll$lat, site_ll$lon, etrs89lat, etrs89long))

all_bbs_squares <- unique(bbs_visits$square)

# --- Temporal effort windows ---------------------------------------------------
# For every BBS square, rank its surveyed years (from bbs_visits, so years with
# zero detections of a species still count as effort) backwards from year_bbs,
# and keep the K_NEEDED most recent within the lookback window. Each kept year
# gets a weight: 1 for the first K_FULL years, K_FRAC for the last one.
# Squares with fewer than K_NEEDED surveyed years in the window are ineligible.
square_years <- bbs_visits %>%
  as_tibble() %>%
  distinct(square, year) %>%
  filter(year <= year_bbs, year > year_bbs - LOOKBACK_MAX_YEARS) %>%
  group_by(square) %>%
  arrange(desc(year), .by_group = TRUE) %>%
  mutate(yr_rank = row_number()) %>%
  filter(yr_rank <= K_NEEDED) %>%
  mutate(n_years_avail = n(),
         weight        = if_else(yr_rank <= K_FULL, 1, K_FRAC)) %>%
  ungroup() %>%
  filter(n_years_avail >= K_NEEDED)

eligible_squares <- unique(square_years$square)
window_years     <- range(square_years$year)
message(sprintf(
  "Temporal window: %d of %d BBS squares have >= %d surveyed years between %d and %d.",
  length(eligible_squares), length(all_bbs_squares), K_NEEDED,
  year_bbs - LOOKBACK_MAX_YEARS + 1, year_bbs))

find_regional_squares <- function(radius_km) {
  grid_coords %>%
    filter(square %in% eligible_squares, dist_km <= radius_km) %>%
    pull(square)
}

radius_used <- RADIUS_MIN_KM
repeat {
  candidate_squares <- find_regional_squares(radius_used)
  sq_per_year <- length(candidate_squares)
  if (sq_per_year >= MIN_SQUARES || radius_used >= RADIUS_MAX_KM) break
  radius_used <- radius_used + RADIUS_STEP_KM
}

regional_squares <- candidate_squares

if (sq_per_year < MIN_SQUARES) {
  warning(sprintf(
    "At maximum radius %.0f km, only %d eligible squares found (< %d target). Consider reducing MIN_SQUARES or K years.",
    RADIUS_MAX_KM, sq_per_year, MIN_SQUARES
  ))
}

n_squares_used <- length(regional_squares)

regional_years <- square_years %>% filter(square %in% regional_squares) %>% pull(year) %>% range()
radius_label <- sprintf("Local BBS squares within %.0f km of %s (n = %d squares, %d\u2013%d)",
                        radius_used, SITE_NAME,
                        n_squares_used,
                        regional_years[1], regional_years[2])
message("Region: ", radius_label)

# =============================================================================
# 6. BBS PSEUDO-TRANSECTS AND ABUNDANCE BENCHMARKS
# =============================================================================
# Each eligible square becomes one pseudo-transect with the same effort as the
# site transect. For each species at each square:
#   pt_count    = sum over kept years of max(E, L) count x year weight
#                 (directly comparable to the site's peak count on 9.2 km)
#   pt_presence = 1 if recorded in a full-weight year, K_FRAC if recorded only
#                 in the partial year, else 0 (expected detection in 9.2 km)
# Benchmarks are then the DISTRIBUTION of these across squares, so the site
# can be placed within it (percentiles), not just compared to a mean.

build_pseudo_transects <- function(squares) {
  sy <- as.data.table(square_years)[square %in% squares, .(square, year, weight)]
  bb <- as.data.table(bbs_by_square)[square %in% squares]
  pt <- bb[sy, on = .(square, year), nomatch = 0L]
  pt <- pt[, .(pt_count    = sum(square_count * weight),
               pt_presence = if (any(square_count > 0)) max(weight[square_count > 0]) else 0),
           by = .(square, species_code)]
  # Zero-fill so every square has a row for every species in the pool
  as_tibble(pt) %>%
    complete(square = squares, species_code,
             fill = list(pt_count = 0, pt_presence = 0))
}

pt_regional <- build_pseudo_transects(regional_squares)
pt_national <- build_pseudo_transects(eligible_squares)

summarise_pt <- function(pt) {
  pt %>%
    group_by(species_code) %>%
    summarise(total_count = sum(pt_count),
              RAI_mean    = mean(pt_count),
              RAI_q25     = quantile(pt_count, 0.25),
              RAI_median  = median(pt_count),
              RAI_q75     = quantile(pt_count, 0.75),
              RAI_q90     = quantile(pt_count, 0.90),
              occupancy   = mean(pt_presence),
              .groups = "drop")
}

n_squares_regional <- n_squares_used
n_squares_national <- length(eligible_squares)

bbs_rai <- summarise_pt(pt_regional) %>%
  rename(bbs_RAI = RAI_mean)

bbs_wide <- bbs_rai %>%
  mutate(rel_abund = bbs_RAI / sum(bbs_RAI)) %>%
  select(species_code, rel_abund)

# National BBS
bbs_rai_national <- summarise_pt(pt_national) %>%
  rename(bbs_RAI_national = RAI_mean)

bbs_national_wide <- bbs_rai_national %>%
  mutate(rel_abund = bbs_RAI_national / sum(bbs_RAI_national)) %>%
  select(species_code, rel_abund)

# Site percentile for each species within the regional pseudo-transects
# (mid-rank for ties: 50% = typical; >90% = among the highest in the region)
site_percentiles <- pt_regional %>%
  inner_join(RAI_site %>% select(species_code, site_RAI = RAI), by = "species_code") %>%
  group_by(species_code) %>%
  summarise(pct_regional = mean(pt_count < site_RAI) + 0.5 * mean(pt_count == site_RAI),
            .groups = "drop")

# Composition typicality: Bray-Curtis similarity of each regional pseudo-transect
# to the regional mean composition, versus the site's similarity to it.
bc_to_mean <- function(pt, mean_comp) {
  pt %>%
    group_by(square) %>%
    mutate(p = if (sum(pt_count) > 0) pt_count / sum(pt_count) else 0) %>%
    ungroup() %>%
    left_join(mean_comp, by = "species_code") %>%
    replace_na(list(rel_abund = 0)) %>%
    group_by(square) %>%
    summarise(bc_sim = sum(pmin(p, rel_abund)), .groups = "drop")
}
bc_regional_dist <- bc_to_mean(pt_regional, bbs_wide)

# =============================================================================
# 7. SPECIES NAME MATCHING
# =============================================================================

unmatched <- setdiff(unique(data$species_code), unique(sp_lookup$species_code))
if (length(unmatched) > 0) {
  message("Species codes not in BBS lookup: ", paste(unmatched, collapse = ", "))
} else {
  message("All species codes matched to BBS lookup.")
}

add_names <- function(df) {
  df %>%
    left_join(sp_names, by = "species_code") %>%
    mutate(english_name = coalesce(english_name, species_code))
}

RAI_site         <- add_names(RAI_site)
RAI_site_all     <- add_names(RAI_site_all)
bbs_rai          <- add_names(bbs_rai)
bbs_rai_national <- add_names(bbs_rai_national)

# =============================================================================
# 8. SPECIES PRESENCE POOLS
# =============================================================================

regional_present <- bbs_rai %>%
  filter(total_count > 0) %>%
  pull(species_code)

national_present <- bbs_rai_national %>%
  filter(total_count > 0) %>%
  pull(species_code)

# =============================================================================
# 9. BRAY-CURTIS SIMILARITY
# =============================================================================

combined_bc <- full_join(
  site_wide  %>% rename(site = rel_abund),
  bbs_wide   %>% rename(bbs  = rel_abund),
  by = "species_code"
) %>% replace_na(list(site = 0, bbs = 0))

mat_bc <- t(select(combined_bc, site, bbs))
colnames(mat_bc) <- combined_bc$species_code

bc_similarity <- suppressWarnings(
  1 - as.numeric(vegdist(mat_bc, method = "bray"))
)
message(sprintf("Bray-Curtis similarity to regional BBS: %.3f", bc_similarity))

# Where does the site sit among regional pseudo-transects? (low % = unusual composition)
bc_site_pct <- mean(bc_regional_dist$bc_sim < bc_similarity)
message(sprintf("Site composition is more typical than %.0f%% of regional BBS pseudo-transects.",
                bc_site_pct * 100))

# =============================================================================
# 10. SPECIES DIVERGENCE
# =============================================================================

divergence <- full_join(
  site_wide         %>% rename(site         = rel_abund),
  bbs_wide          %>% rename(bbs_regional = rel_abund),
  by = "species_code"
) %>%
  full_join(bbs_national_wide %>% rename(bbs_national = rel_abund),
            by = "species_code") %>%
  replace_na(list(site = 0, bbs_regional = 0, bbs_national = 0)) %>%
  add_names() %>%
  mutate(
    diff_regional      = site - bbs_regional,
    diff_national      = site - bbs_national,
    direction_regional = case_when(
      site > bbs_regional ~ "Over-represented",
      site < bbs_regional ~ "Under-represented",
      TRUE                ~ "Equal"
    ),
    direction_national = case_when(
      site > bbs_national ~ "Over-represented",
      site < bbs_national ~ "Under-represented",
      TRUE                ~ "Equal"
    )
  )

# =============================================================================
# 11. SPECIES-LEVEL ABUNDANCE COMPARISON
# =============================================================================

species_comparison <- RAI_site %>%
  left_join(bbs_rai          %>% select(species_code, bbs_RAI, RAI_q25, RAI_median, RAI_q75, RAI_q90),
            by = "species_code") %>%
  left_join(bbs_rai_national %>% select(species_code, bbs_RAI_national),
            by = "species_code") %>%
  left_join(site_percentiles, by = "species_code") %>%
  rename(regional_mean_RAI = bbs_RAI,
         national_mean_RAI = bbs_RAI_national) %>%
  mutate(
    ratio_vs_regional = RAI / regional_mean_RAI,
    ratio_vs_national = RAI / national_mean_RAI
  ) %>%
  filter(!is.na(regional_mean_RAI)) %>%
  arrange(desc(abs(log(ratio_vs_regional))))

# =============================================================================
# 12. SPECIES RICHNESS — SAME-EFFORT, SAME-PLACE PSEUDO-TRANSECTS
# =============================================================================
# Each regional pseudo-transect is one BBS square surveyed over N_YEARS_EQUIV
# years, i.e. the same length of transect as the site, but from ONE location.
# Its richness is sum(pt_presence): species recorded in full-weight years count
# 1, species recorded only in the partial year count K_FRAC. The benchmark is
# the distribution of this across squares (no spatial pooling), and the site is
# placed as a percentile within it. Note that year-to-year turnover still adds
# some species (as repeat visits do at the site), but far fewer than mixing
# habitats across different squares.

pt_richness <- function(pt) {
  pt %>% group_by(square) %>% summarise(S = sum(pt_presence), .groups = "drop")
}
richness_regional_dist <- pt_richness(pt_regional)
richness_national_dist <- pt_richness(pt_national)

summarise_richness <- function(d) {
  list(S_expected = mean(d$S),
       CI_low     = unname(quantile(d$S, 0.025)),
       CI_high    = unname(quantile(d$S, 0.975)),
       pct_site   = mean(d$S < site_total_species) + 0.5 * mean(d$S == site_total_species))
}
regional_rare <- summarise_richness(richness_regional_dist)
national_rare <- summarise_richness(richness_national_dist)
n_units_equiv <- N_YEARS_EQUIV

richness_comparison <- tibble(
  Dataset          = c("Your site", "BBS regional pseudo-transects", "BBS national pseudo-transects"),
  Samples          = c(1, n_squares_regional, n_squares_national),
  Species_pool     = c(site_total_species,
                       sum(bbs_rai$total_count > 0),
                       sum(bbs_rai_national$total_count > 0)),
  Rarefied_S       = c(site_total_species,
                       round(regional_rare$S_expected, 1),
                       round(national_rare$S_expected, 1)),
  CI_low           = c(NA, regional_rare$CI_low,  national_rare$CI_low),
  CI_high          = c(NA, regional_rare$CI_high, national_rare$CI_high),
  Site_percentile  = c(NA, round(regional_rare$pct_site * 100), round(national_rare$pct_site * 100))
)
print(richness_comparison)
# =============================================================================
# 13. PRESENCE/ABSENCE OVERLAP
# =============================================================================

calc_overlap <- function(set_a, set_b, label_b) {
  shared <- length(intersect(set_a, set_b))
  only_a <- length(setdiff(set_a, set_b))
  only_b <- length(setdiff(set_b, set_a))
  tibble(
    Comparison        = paste("Site vs", label_b),
    Site_species      = length(set_a),
    Reference_species = length(set_b),
    Shared            = shared,
    Site_only         = only_a,
    Reference_only    = only_b,
    Jaccard           = round(shared / (shared + only_a + only_b), 3),
    Sorensen          = round(2 * shared / (length(set_a) + length(set_b)), 3)
  )
}

# Species withheld from BBS but assumed present (e.g. Red Kite) count as shared
site_present_pool <- union(site_present, assumed_present)
overlap_results <- bind_rows(
  calc_overlap(site_present_pool, union(regional_present, assumed_present), "BBS regional"),
  calc_overlap(site_present_pool, union(national_present, assumed_present), "BBS national")
)
print(overlap_results)

# =============================================================================
# 14. SPECIES GAP LISTS & OCCUPANCY
# =============================================================================

regional_absent <- bbs_rai %>%
  filter(species_code %in% regional_present,
         !species_code %in% site_present) %>%
  arrange(desc(bbs_RAI))

site_only_species <- sp_names %>%
  filter(species_code %in% site_present,
         !species_code %in% regional_present) %>%
  left_join(RAI_site %>% select(species_code, RAI), by = "species_code") %>%
  mutate(english_name = coalesce(english_name, species_code)) %>%
  arrange(desc(RAI))

# Regional occupancy: share of regional pseudo-transects (one square, same
# effort as the site) in which the species was recorded
regional_occupancy <- pt_regional %>%
  group_by(species_code) %>%
  summarise(
    n_squares_with = sum(pt_presence > 0),
    occupancy      = mean(pt_presence),
    .groups = "drop"
  ) %>%
  filter(n_squares_with > 0) %>%
  add_names() %>%
  mutate(recorded_at_site = species_code %in% site_present) %>%
  arrange(desc(occupancy))

absent_with_occupancy <- regional_occupancy %>%
  filter(!recorded_at_site) %>%
  left_join(regional_absent %>% select(species_code, bbs_RAI), by = "species_code") %>%
  arrange(desc(occupancy))

message(sprintf("%d regional species absent from site.", nrow(regional_absent)))
message(sprintf("%d site species absent from regional BBS.", nrow(site_only_species)))

# =============================================================================
# 15a. SPECIES TRAITS: GUILD, FARMLAND INDEX MEMBERSHIP, BOCC5 STATUS
# =============================================================================

bocc_levels <- c("Red", "Amber", "Green", "Non-native", "Not listed")

species_codes <- read_csv(SPECIES_CODES_FILE, col_types = cols(.default = "c")) %>%
  transmute(species_code = str_trim(bto_code),
            BOCC5 = bocc_levels[match(tolower(str_trim(bocc)), tolower(bocc_levels))]) %>%
  filter(!is.na(species_code)) %>%
  distinct(species_code, .keep_all = TRUE)

if (!is.na(GUILD_FILE) && file.exists(GUILD_FILE)) {
  guild_table <- read_csv(GUILD_FILE, col_types = cols(.default = "c")) %>%
    transmute(species_code   = str_trim(bto_code),
              guild          = na_if(str_to_title(str_trim(guild)), ""),
              farmland_index = na_if(str_to_title(str_trim(farmland_index)), ""),
              across(any_of(c("diet", "nest_site", "migration")), ~ na_if(str_trim(.x), ""))) %>%
    filter(!is.na(species_code)) %>%
    distinct(species_code, .keep_all = TRUE)
} else {
  warning("Guild file not found: ", GUILD_FILE,
          "\nContinuing without guilds - all species will be grouped as 'Other' in Figure 11 ",
          "and the farmland bird score (Figure 13) will be skipped. Save species_guilds.csv at that path to enable them.")
  guild_table <- tibble(species_code = character(), guild = character(), farmland_index = character())
}

bad_index <- setdiff(na.omit(guild_table$farmland_index), c("Specialist", "Generalist"))
if (length(bad_index) > 0)
  warning("Unrecognised farmland_index values in guild file (ignored): ", paste(bad_index, collapse = ", "))

species_traits <- full_join(species_codes, guild_table, by = "species_code") %>%
  mutate(guild = coalesce(guild, "Other"),
         fbi   = case_when(farmland_index == "Specialist" ~ "S",
                           farmland_index == "Generalist" ~ "G"),
         BOCC5 = coalesce(BOCC5, "Not listed"))

guild_levels <- c(intersect(str_to_title(GUILD_ORDER), unique(species_traits$guild)),
                  sort(setdiff(unique(species_traits$guild), c(str_to_title(GUILD_ORDER), "Other"))),
                  "Other")

# Report species recorded at the site or regionally that the lookups don't cover
uncovered <- setdiff(union(site_present, unique(pt_regional$species_code)), species_codes$species_code)
if (length(uncovered) > 0)
  message("Species codes not in the BTO lookup (BOCC5 = 'Not listed'): ", paste(uncovered, collapse = ", "))

# Species-wide lookup (used by Figure 8c and 15b), covering absent species too
bocc_lookup <- species_traits %>% select(species_code, BOCC5)

# Percentile of every regional species (site count 0 if not recorded) and of
# every regional square for every species (mid-rank, so ties split evenly)
site_counts_all <- tibble(species_code = unique(pt_regional$species_code)) %>%
  left_join(RAI_site %>% select(species_code, site_RAI = RAI), by = "species_code") %>%
  replace_na(list(site_RAI = 0))

species_bench <- pt_regional %>%
  inner_join(site_counts_all, by = "species_code") %>%
  group_by(species_code) %>%
  summarise(site_RAI   = first(site_RAI),
            pct        = mean(pt_count < site_RAI) + 0.5 * mean(pt_count == site_RAI),
            reg_median = median(pt_count),
            reg_mean   = mean(pt_count),
            occupancy  = mean(pt_presence),
            .groups = "drop") %>%
  filter(reg_mean > 0 | site_RAI > 0) %>%
  add_names() %>%
  left_join(species_traits, by = "species_code") %>%
  mutate(guild = coalesce(guild, "Other"), BOCC5 = coalesce(BOCC5, "Not listed"),
         recorded = species_code %in% site_present)

square_pct <- pt_regional %>%
  group_by(species_code) %>%
  mutate(pct = (rank(pt_count, ties.method = "average") - 0.5) / n()) %>%
  ungroup()

# =============================================================================
# 15. FIGURES
# =============================================================================

fig_title <- function(txt) txt
occ_axis_label <- "Regional occupancy (% of BBS squares where species is present)"

# --- Figure 1: Site RAI ---
p1 <- RAI_site_all %>%
  filter(total_count > 0) %>%
  mutate(english_name = fct_reorder(english_name, RAI)) %>%
  ggplot(aes(x = RAI, y = english_name)) +
  geom_col(fill = "steelblue") +
  labs(
    title    = fig_title("Site Relative Abundance Index"),
    subtitle = sprintf("Maximum count across visits along %.1f km transect", SITE_TRANSECT_KM),
    x = "Max count across visits", y = NULL
  )
save_fig(p1, "fig1_site_RAI.png", width = 18, height = 22)


# --- Figures 3 and 4: Difference in count, site vs BBS mean ---
# Bars: site count (highest visit, over the transect) minus the mean count on BBS
# pseudo-transects of the same length, so each bar depends on that species alone.
# Labels: the same difference as a percentage of the BBS mean.
count_diff_plot <- function(d, title, subtitle, where, units) {
  # d needs english_name, site (site count) and bench (benchmark mean count).
  # Bars = site - bench, for the 20 species with the largest difference in birds;
  # the label on each bar is the same difference as a percentage of the benchmark.
  d <- d %>% mutate(diff = site - bench) %>% filter(!is.na(diff)) %>%
    slice_max(abs(diff), n = 20, with_ties = FALSE) %>%
    mutate(english_name = fct_reorder(english_name, diff),
           pct = if_else(bench > 0, (site / bench - 1) * 100, NA_real_),
           lab = if_else(is.na(pct), "none in benchmark",
                         paste0(ifelse(pct > 0, "+", ""), comma(pct, accuracy = 1), "%")),
           direction = if_else(diff >= 0, paste("More than", where), paste("Fewer than", where)))
  ggplot(d, aes(x = diff, y = english_name, fill = direction)) +
    geom_col() +
    geom_text(aes(label = lab, hjust = ifelse(diff >= 0, -0.1, 1.1)), size = 3, colour = "grey20") +
    scale_fill_manual(values = setNames(c("steelblue", "tomato"), c(paste("More than", where), paste("Fewer than", where)))) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_x_continuous(labels = function(x) ifelse(x > 0, paste0("+", comma(x)), comma(x)),
                       expand = expansion(mult = c(0.14, 0.14))) +
    labs(title = title, subtitle = subtitle,
         x = units, y = NULL, fill = NULL) +
    theme(legend.position = "bottom")
}
site_vs_bbs <- species_comparison %>% filter(species_code %in% site_present)   # only species recorded on the farm
diff_note  <- "Labels: the difference as a percentage of the BBS mean (+100% = twice the BBS mean; -50% = half)."
diff_units <- "Difference in abundance"   # x-axis title

# --- Figure 3: Difference from regional ---
p3 <- count_diff_plot(
  site_vs_bbs %>% transmute(english_name, site = RAI, bench = regional_mean_RAI),
  fig_title("Differences in abundance compared to the regional bird community"),
  paste0(str_wrap(paste0("The 20 recorded species with the largest differences in number of birds. Bars: site count minus the regional BBS mean count. ",
                         diff_note), 95), "\n", radius_label),
  "the region", diff_units)
save_fig(p3, "fig3_divergence_regional.png", height = 16)

# --- Figure 4: Difference from national ---
p4 <- count_diff_plot(
  site_vs_bbs %>% transmute(english_name, site = RAI, bench = national_mean_RAI),
  fig_title("Differences in abundance compared to the national bird community"),
  str_wrap(paste0("The 20 recorded species with the largest differences in number of birds. Bars: site count minus the national BBS mean count. ",
                  diff_note, " National = full UK BBS."), 95),
  "the UK", diff_units)
save_fig(p4, "fig4_divergence_national.png", height = 16)

# --- Figure 5: Abundance ratio (log2) ---
p5 <- species_comparison %>%
  filter(total_count > 0) %>%
  mutate(
    log2_regional = log2(ratio_vs_regional),
    log2_national = log2(ratio_vs_national),
    english_name  = fct_reorder(english_name, log2_regional)
  ) %>%
  select(english_name, log2_regional, log2_national) %>%
  pivot_longer(cols = c(log2_regional, log2_national),
               names_to = "reference", values_to = "log2_ratio") %>%
  mutate(reference = recode(reference,
                            "log2_regional" = "vs BBS regional",
                            "log2_national" = "vs BBS national"
  )) %>%
  ggplot(aes(x = log2_ratio, y = english_name, colour = reference, shape = reference)) +
  geom_vline(xintercept = 0,      linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-1,1), linetype = "dotted", colour = "grey80") +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  labs(
    title    = fig_title("Species abundance ratio: site vs BBS"),
    subtitle = paste0("0 = equal; +1 = twice as abundant at site; \u22121 = half.\n", str_wrap(effort_note, 110)),
    x = "Log2(site RAI / BBS RAI)", y = NULL, colour = NULL, shape = NULL
  ) +
  theme(legend.position = "bottom")
save_fig(p5, "fig5_abundance_ratio.png", width = 20, height = 18)

# --- Figure 6a: Species richness — site vs expected at the same effort ---
# Bars = farm richness and the mean richness of same-effort BBS samples
# (one square over N_YEARS_EQUIV years = the farm transect length).
# Error bars = the range containing 95% of those samples.
p6a <- richness_comparison %>%
  mutate(Dataset = recode(Dataset,
                          "Your site" = SITE_NAME,
                          "BBS regional pseudo-transects" = "BBS regional\n(expected)",
                          "BBS national pseudo-transects" = "BBS national\n(expected)"),
         Dataset = fct_inorder(Dataset),
         lab = if_else(is.na(CI_low), sprintf("%d", round(Rarefied_S)),
                       sprintf("%.1f\n(%.0f\u2013%.0f)", Rarefied_S, CI_low, CI_high))) %>%
  ggplot(aes(x = Dataset, y = Rarefied_S, fill = Dataset)) +
  geom_col(width = 0.55) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.6, na.rm = TRUE) +
  geom_text(aes(y = pmax(Rarefied_S, coalesce(CI_high, 0)), label = lab),
            vjust = -0.4, size = 3.4, lineheight = 0.9) +
  scale_fill_manual(values = c("#2166AC", "tomato", "grey50"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = fig_title(paste0("Species richness: ", SITE_NAME, " vs BTO BBS")),
       subtitle = str_wrap(paste0(
         "Expected = mean species recorded by BBS squares surveyed with the same effort as the farm (",
         SITE_TRANSECT_KM, " km: one square over ", round(N_YEARS_EQUIV, 1), " years). ",
         "Error bars = range containing 95% of squares. Regional n = ", n_squares_regional,
         " squares; national n = ", n_squares_national, " squares."), 90),
       x = NULL, y = "Number of species") +
  theme_minimal(base_size = 10) + theme(plot.title.position = "plot", panel.grid.major.x = element_blank())
save_fig(p6a, "fig6a_richness_comparison.png", width = 16, height = 13)

# --- Figure 6d: Richness — raw same-effort samples as box plots + points ---
rich_cols <- c(SITE_NAME, "BBS regional", "BBS national")
rich_raw <- bind_rows(
  richness_regional_dist %>% transmute(group = "BBS regional", S),
  richness_national_dist %>% transmute(group = "BBS national", S)
) %>% mutate(group = factor(group, levels = rich_cols))
rich_farm <- tibble(group = factor(SITE_NAME, levels = rich_cols), S = site_total_species)
# Box = 25th-75th percentile + median; whiskers = 2.5th-97.5th percentile
# (the range containing 95% of squares, matching the error bars in Fig 6a)
rich_box <- rich_raw %>%
  group_by(group) %>%
  summarise(ymin = quantile(S, 0.025), lower = quantile(S, 0.25), middle = median(S),
            upper = quantile(S, 0.75), ymax = quantile(S, 0.975), .groups = "drop")
rich_n <- tibble(group = factor(c("BBS regional", "BBS national"), levels = rich_cols),
                 lab = sprintf("n = %d squares", c(n_squares_regional, n_squares_national)))

p6d <- ggplot(rich_raw, aes(x = group, y = S)) +
  geom_hline(yintercept = site_total_species, colour = "#2166AC", linetype = "dashed", linewidth = 0.5) +
  geom_boxplot(data = rich_box, inherit.aes = FALSE, stat = "identity",
               aes(x = group, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = group),
               width = 0.5, alpha = 0.35, colour = "grey30") +
  geom_point(aes(colour = group), position = position_jitter(width = 0.13, height = 0, seed = 1),
             size = 1.8, alpha = 0.75) +
  geom_point(data = rich_farm, shape = 23, size = 6, fill = "#2166AC", colour = "white", stroke = 1.2) +
  geom_text(data = rich_farm, aes(label = paste(S, "species")), vjust = -1.6,
            colour = "#2166AC", fontface = "bold", size = 3.6) +
  geom_text(data = rich_n, aes(label = lab, y = -Inf), vjust = -0.6, size = 3, colour = "grey40") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("BBS regional" = "tomato", "BBS national" = "grey50"), guide = "none") +
  scale_colour_manual(values = c("BBS regional" = "tomato3", "BBS national" = "grey35"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = fig_title("Species richness: site vs same-effort BBS samples"),
       subtitle = str_wrap(paste0(
         "Each dot is one BBS square surveyed with the same effort as the farm (", SITE_TRANSECT_KM,
         " km: one square over ", round(N_YEARS_EQUIV, 1), " years). Boxes show the middle 50% of squares and the median; ",
         "whiskers show the range containing 95% of squares (2.5th\u201397.5th percentile). Dashed line = the farm's richness. The farm has more species than ",
         round(regional_rare$pct_site * 100), "% of regional and ", round(national_rare$pct_site * 100),
         "% of national squares."), 78),
       x = NULL, y = "Number of species") +
  theme_minimal(base_size = 10) + theme(plot.title.position = "plot", panel.grid.major.x = element_blank(),
                                        axis.text.x = element_text(face = "bold", size = 10))
save_fig(p6d, "fig6d_richness_boxplot.png", width = 16, height = 13)

# --- Figure 6: Richness — where the site sits among same-effort samples ---
richness_dists <- bind_rows(
  richness_regional_dist %>% mutate(scale = "Regional"),
  richness_national_dist %>% mutate(scale = "National")
) %>% mutate(scale = factor(scale, levels = c("Regional", "National")))

richness_labels <- tibble(
  scale = factor(c("Regional", "National"), levels = c("Regional", "National")),
  lab   = sprintf("Site: %d species\nhigher than %.0f%% of %s BBS squares\nn = %d squares, %d\u2013%d",
                  site_total_species,
                  c(regional_rare$pct_site, national_rare$pct_site) * 100,
                  c("regional", "national"),
                  c(n_squares_regional, n_squares_national),
                  c(regional_years[1], window_years[1]),
                  c(regional_years[2], window_years[2]))
)

p6 <- richness_dists %>%
  ggplot(aes(x = S)) +
  geom_histogram(binwidth = 2, fill = "grey75", colour = "white", boundary = 0) +
  geom_vline(xintercept = site_total_species, colour = "steelblue", linewidth = 1.1) +
  geom_label_nb(data = richness_labels, aes(x = -Inf, y = Inf, label = lab),
                hjust = 0, vjust = 1, size = 3.2, colour = "steelblue", inherit.aes = FALSE,
                fill = alpha("white", 0.75), label.padding = unit(0.4, "lines")) +
  facet_wrap(~ scale, ncol = 1, scales = "free_y") +
  labs(
    title    = fig_title(paste0("Species richness at ", SITE_NAME, " vs BTO BBS")),
    subtitle = paste0("Each grey unit is one BBS square surveyed over ", round(N_YEARS_EQUIV, 1),
                      " years (= ", SITE_TRANSECT_KM, " km of transect).\n",
                      sprintf("Regional = local BBS squares within %.0f km of %s.", radius_used, SITE_NAME)),
    x = "Species recorded per BBS square", y = "Number of BBS squares"
  )
save_fig(p6, "fig6_richness_distribution.png", width = 18, height = 16)

# --- Figure 6b: Species percentile strip — where each species' count sits ---
p6b <- species_comparison %>%
  filter(!is.na(pct_regional)) %>%
  mutate(english_name = fct_reorder(english_name, pct_regional)) %>%
  ggplot(aes(y = english_name)) +
  annotate("rect", xmin = 0.25, xmax = 0.75, ymin = -Inf, ymax = Inf, fill = "grey75", alpha = 0.35) +
  geom_vline(xintercept = 0.5, colour = "grey60", linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = pct_regional, yend = english_name), colour = "grey70") +
  geom_point(aes(x = pct_regional, colour = pct_regional), size = 3) +
  scale_colour_gradient2(low = "tomato", mid = "grey50", high = "steelblue",
                         midpoint = 0.5, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = fig_title("How the site's counts rank against the region"),
    subtitle = paste0("Percentile of the site's peak count among regional ", SITE_TRANSECT_KM,
                      " km pseudo-transects.\nShaded band = middle half of the region.\n",
                      str_wrap(radius_label, 80)),
    x = "Site percentile (0% = lowest in region, 100% = highest)", y = NULL
  )
save_fig(p6b, "fig6b_species_percentiles.png", width = 18,
         height = max(10, n_distinct(p6b$data$english_name) * 0.45))

# --- Figure 6c: Compositional typicality ---
p6c <- bc_regional_dist %>%
  ggplot(aes(x = bc_sim)) +
  geom_histogram(binwidth = 0.025, fill = "grey75", colour = "white", boundary = 0) +
  geom_vline(xintercept = bc_similarity, colour = "steelblue", linewidth = 1.1) +
  annotate("text", x = bc_similarity, y = Inf, vjust = 1.5, size = 3.2,
           hjust = if (bc_similarity > median(bc_regional_dist$bc_sim)) 1.05 else -0.05,
           colour = "steelblue",
           label = sprintf("Site: %.2f\nmore typical than %.0f%% of squares",
                           bc_similarity, bc_site_pct * 100)) +
  labs(
    title    = fig_title("How typical is the site's bird community?"),
    subtitle = paste0("Bray-Curtis similarity of each regional pseudo-transect to the\nregional average composition.\n",
                      str_wrap(radius_label, 80)),
    x = "Similarity to regional average community (1 = identical)", y = "Number of BBS squares"
  )
save_fig(p6c, "fig6c_composition_typicality.png", width = 18, height = 12)


# --- Figure 7: Species overlap stacked bar ---
overlap_long <- overlap_results %>%
  select(Comparison, Shared, Site_only, Reference_only) %>%
  pivot_longer(cols = c(Shared, Site_only, Reference_only),
               names_to = "category", values_to = "n_species") %>%
  mutate(category = recode(category,
                           "Shared"         = "Recorded at both",
                           "Site_only"      = "Site only",
                           "Reference_only" = "BBS only"
  ))

p7 <- overlap_long %>%
  ggplot(aes(x = Comparison, y = n_species, fill = category)) +
  geom_col(position = "stack", width = 0.5) +
  scale_fill_manual(values = c("Recorded at both" = "steelblue",
                               "Site only"        = "goldenrod",
                               "BBS only"         = "tomato")) +
  geom_text(aes(label = n_species),
            position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  labs(
    title    = fig_title("Species pool overlap: site vs BBS"),
    subtitle = paste0(
      "Jaccard (regional): ", overlap_results$Jaccard[1],
      " | S\u00f8rensen: ", overlap_results$Sorensen[1], "\n",
      "Jaccard (national): ", overlap_results$Jaccard[2],
      " | S\u00f8rensen: ", overlap_results$Sorensen[2]
    ),
    x = NULL, y = "Number of species", fill = NULL
  ) +
  theme(legend.position = "bottom")
save_fig(p7, "fig7_species_overlap.png", width = 16, height = 12)

# --- Figure 8: Absent species ranked by regional RAI ---
p8 <- regional_absent %>%
  slice_max(bbs_RAI, n = 30) %>%
  mutate(english_name = fct_reorder(english_name, bbs_RAI)) %>%
  ggplot(aes(x = bbs_RAI, y = english_name)) +
  geom_col(fill = "tomato") +
  labs(
    title    = fig_title("Regional species not recorded at site"),
    subtitle = paste0("Top 30 of ", nrow(regional_absent),
                      " regional species absent from site. Ranked by regional abundance.\n",
                      radius_label),
    x = sprintf("Mean BBS peak count per %.1f km pseudo-transect", SITE_TRANSECT_KM),
    y = NULL
  )
save_fig(p8, "fig8_absent_species_ranked.png", width = 18, height = 20)

# --- Figure 8b: Site Species Breakdown by Regional Occupancy Tier ---
occupancy_tier_data <- regional_occupancy %>%
  filter(species_code %in% site_present) %>%
  mutate(
    occupancy_tier = case_when(
      occupancy >= 0.70 ~ "Ubiquitous (>70% regional pseudo-transects)",
      occupancy >= 0.40 ~ "Widespread (40-70%)",
      occupancy >= 0.15 ~ "Localized (15-40%)",
      TRUE              ~ "Scarce / Rare (<15%)"
    ),
    occupancy_tier = factor(occupancy_tier, levels = c(
      "Ubiquitous (>70% regional pseudo-transects)",
      "Widespread (40-70%)",
      "Localized (15-40%)",
      "Scarce / Rare (<15%)"
    ))
  ) %>%
  count(occupancy_tier, name = "number_of_species")

# =============================================================================
# FIGURE 8b: Top Regional Species Absent from Site (Sorted by Occupancy)
# =============================================================================

# Identify species present in the regional pool but absent from your site
absent_species_data <- regional_occupancy %>%
  filter(!species_code %in% site_present) %>%
  # Sort descending by regional occupancy so the commonest missing species are at the top
  arrange(desc(occupancy)) %>%
  slice_head(n = 40) %>%
  # Arrange ascending for correct ggplot ordering (commonest at top)
  arrange(occupancy)

p8b <- absent_species_data %>%
  ggplot(aes(x = occupancy, y = fct_inorder(english_name))) +
  geom_col(fill = "tomato3", width = 0.6) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = occ_axis_label, y = NULL
  )

save_fig(p8b, "fig8b_top_regional_absences.png", width = 18, height = 13)

# =============================================================================
# FIGURE 8c: Red & Amber List Regional Species Absent from Site (BOCC5)
# =============================================================================

# bocc_lookup comes from section 15a and covers species absent from the site
# (previously it came from site data only, so absent species had no status)

# Join regional occupancy with our cleaned BOCC5 lookup
red_amber_absences <- regional_occupancy %>%
  left_join(bocc_lookup, by = "species_code") %>%
  filter(!species_code %in% site_present) %>%
  # Filter strictly for Red and Amber listed species using BOCC5
  filter(BOCC5 %in% c("Red", "Amber")) %>%
  # Sort descending by regional occupancy so the commonest missing priority species are at the top
  arrange(desc(occupancy)) %>%
  slice_head(n = 30) %>%
  # Arrange ascending for correct ggplot ordering (commonest at top)
  arrange(occupancy)

p8c <- red_amber_absences %>%
  ggplot(aes(x = occupancy, y = fct_inorder(english_name), fill = BOCC5)) +
  geom_col(width = 0.6) +
  # Use standard conservation colors: deep red for Red list, warm orange/amber for Amber list
  scale_fill_manual(values = c("Red" = "#D95F02", "Amber" = "#E6AB02")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = occ_axis_label, y = NULL, fill = "BOCC5 Status"
  ) +
  theme(legend.position = "bottom")

save_fig(p8c, "fig8c_priority_species_absences.png", width = 18, height = 12)

# --- Figure 9: Occupancy overview ---
p9 <- regional_occupancy %>%
  slice_max(occupancy, n = 40) %>%
  mutate(
    english_name     = fct_reorder(english_name, occupancy),
    recorded_at_site = factor(recorded_at_site,
                              levels = c(TRUE, FALSE),
                              labels = c("Recorded at site", "Absent from site"))
  ) %>%
  ggplot(aes(x = occupancy, y = english_name, fill = recorded_at_site)) +
  geom_col() +
  scale_fill_manual(values = c("Recorded at site" = "steelblue",
                               "Absent from site" = "tomato")) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = fig_title("Regional occupancy: top 40 species"),
    subtitle = paste0("Blue = recorded at site; Red = absent. ", radius_label),
    x = occ_axis_label, y = NULL, fill = NULL
  ) +
  theme(legend.position = "bottom")
save_fig(p9, "fig9_occupancy_overview.png", width = 18, height = 22)

# --- Figure 10: Regional occupancy of site species ---
p10 <- regional_occupancy %>%
  filter(species_code %in% site_present) %>%
  # Sort explicitly: arrange ascending so rarest is first, then convert to factor in that order
  arrange(occupancy) %>%
  mutate(
    english_name = factor(english_name, levels = english_name)
  ) %>%
  ggplot(aes(x = occupancy, y = english_name)) +
  geom_col(fill = "steelblue", width = 0.6) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = fig_title("Regional Occupancy of Site Species"),
    subtitle = paste0("Showing how common or rare each recorded species is across the wider countryside. Rarest regionally at the top.\n", radius_label),
    x = occ_axis_label, y = NULL
  )
save_fig(p10, "fig10_site_species_occupancy.png", width = 18, height = max(10, nrow(filter(regional_occupancy, species_code %in% site_present)) * 0.4))


# =============================================================================
# 15b. CLIENT-FACING FIGURES (11 heat strip, 13 farmland bird score)
# =============================================================================

pct_fill <- scale_fill_gradient2(
  low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0.5,
  limits = c(0, 1), breaks = c(0.05, 0.5, 0.95),
  labels = c("Among the lowest\nin the region", "Typical", "Among the highest\nin the region"),
  name = NULL)

# --- Figure 11: Percentile heat strip by guild --------------------------------
# All species recorded at the site, plus regional species a same-effort survey
# would usually find (occupancy >= 25%) even if the site missed them.
strip_data <- species_bench %>%
  filter(recorded | occupancy >= 0.25) %>%
  mutate(guild = factor(guild, levels = guild_levels),
         label = if_else(recorded,
                         sprintf("%s  (region median %s)", comma(site_RAI, 1), comma(reg_median, 0.1)),
                         sprintf("not recorded  (region median %s)", comma(reg_median, 0.1))),
         english_name = fct_reorder(english_name, pct))

p11 <- strip_data %>%
  ggplot(aes(x = 1, y = english_name, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label, colour = abs(pct - 0.5) > 0.3), size = 2.8, hjust = 0.5) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  pct_fill +
  facet_grid(guild ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = fig_title("How each species on the farm compares with the region"),
       subtitle = str_wrap(paste0(
         "Colour = where the farm's peak count ranks among regional BBS squares surveyed with the same effort (",
         SITE_TRANSECT_KM, " km). Text = farm count (regional median). Includes regionally frequent species the farm did not record."), 95),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) + theme(plot.title.position = "plot") +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
        legend.position = "bottom", legend.key.width = unit(2.2, "cm"))
save_fig(p11, "fig11_percentile_heat_strip.png", width = 18,
         height = max(12, nrow(strip_data) * 0.5 + 5))

ordinal <- function(n) paste0(n, ifelse(n %% 100 %in% 11:13, "th",
                                        c("th","st","nd","rd","th","th","th","th","th","th")[n %% 10 + 1]))

# --- Figure 13: Farmland Bird Index thermometer -------------------------------
# Index = geometric mean over FBI species of (count + 1) / (regional mean + 1).
# 1 = matches the regional average; 2 = twice as many farmland birds. The same
# index is computed for every regional square so the farm can be ranked.
fbi_sp <- species_bench %>% filter(!is.na(fbi), fbi %in% c("S", "G")) %>%
  select(species_code, fbi, reg_mean)

if (nrow(fbi_sp) >= 3) {
  fbi_index <- function(counts, sp) {
    d <- counts %>% inner_join(sp, by = "species_code")
    d %>% group_by(unit) %>%
      summarise(index = exp(mean(log((count + 1) / (reg_mean + 1)))), .groups = "drop")
  }
  fbi_sets <- list("All FBI species" = fbi_sp,
                   "FBI specialists" = filter(fbi_sp, fbi == "S"))
  fbi_sets <- fbi_sets[map_int(fbi_sets, nrow) >= 3]
  
  fbi_res <- imap_dfr(fbi_sets, function(sp, nm) {
    sq <- fbi_index(pt_regional %>% transmute(unit = square, species_code, count = pt_count), sp)
    fm <- fbi_index(site_counts_all %>% transmute(unit = "farm", species_code, count = site_RAI), sp)
    bind_rows(sq %>% mutate(type = "square"), fm %>% mutate(type = "farm")) %>%
      mutate(set = sprintf("%s (%d)", nm, nrow(sp)))
  }) %>% mutate(set = fct_rev(fct_inorder(set)))
  
  fbi_farm <- fbi_res %>% filter(type == "farm")
  fbi_farm$rank <- map2_int(as.character(fbi_farm$set), fbi_farm$index, function(st, v)
    1L + sum(fbi_res$index[fbi_res$type == "square" & fbi_res$set == st] > v))
  fbi_bands <- fbi_res %>% filter(type == "square") %>% group_by(set) %>%
    summarise(q10 = quantile(index, .1), q25 = quantile(index, .25),
              q75 = quantile(index, .75), q90 = quantile(index, .9), .groups = "drop")
  
  p13 <- ggplot() +
    geom_vline(xintercept = 1, colour = "grey50", linetype = "dashed") +
    geom_linerange(data = fbi_bands, aes(y = set, xmin = q10, xmax = q90),
                   linewidth = 9, colour = "grey88") +
    geom_linerange(data = fbi_bands, aes(y = set, xmin = q25, xmax = q75),
                   linewidth = 9, colour = "grey72") +
    geom_point(data = filter(fbi_res, type == "square"), aes(x = index, y = set),
               position = position_jitter(height = 0.12, width = 0, seed = 1),
               colour = "grey35", size = 1.2, alpha = 0.7) +
    geom_point(data = fbi_farm, aes(x = index, y = set), shape = 21, fill = "#1B7837",
               colour = "white", size = 7, stroke = 1.2) +
    geom_text(data = fbi_farm, aes(x = index, y = set,
                                   label = sprintf("Farm %.2fx\n%s of %d", index, map_chr(rank, ordinal),
                                                   n_squares_regional + 1)),
              vjust = -0.9, size = 3.2, fontface = "bold", colour = "#1B7837", lineheight = 0.9) +
    scale_x_log10(labels = function(x) paste0(x, "x")) +
    labs(title = fig_title("Farmland bird score vs the region"),
         subtitle = str_wrap(paste0(
           "Geometric mean across Farmland Bird Index species of (count + 1) / (regional mean + 1). ",
           "1x = regional average. Dots = regional BBS squares at the same effort; ",
           "dark band = middle 50%, light band = middle 80%."), 85),
         x = "Farmland bird score (log scale)", y = NULL) +
    theme_minimal(base_size = 10) + theme(plot.title.position = "plot") +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold"))
  save_fig(p13, "fig13_fbi_thermometer.png", width = 18, height = 9)
} else {
  message("Too few FBI species in the regional pool — Figure 13 skipped.")
}

# --- Figure captions: the titles and subtitles removed from the images ---
write_csv(figure_captions, file.path(out_dir, "figure_captions.csv"))
cap_doc <- read_docx() %>%
  body_add_par(paste0("Figure captions \u2014 ", SITE_NAME, ", ", SEASON_NAME), style = "heading 1")
for (i in seq_len(nrow(figure_captions))) {
  cap_doc <- cap_doc %>%
    body_add_par(figure_captions$file[i], style = "heading 3") %>%
    body_add_fpar(fpar(ftext(figure_captions$title[i], fp_text(bold = TRUE, font.size = 11)))) %>%
    body_add_par(figure_captions$notes[i], style = "Normal")
}
print(cap_doc, target = file.path(out_dir, "figure_captions.docx"))
message("Figure captions saved: figure_captions.docx and figure_captions.csv")

message("All figures saved.")

# =============================================================================
# 16. TABLES
# =============================================================================

# Table 1: Site RAI
ft1 <- RAI_site_all %>%
  filter(total_count > 0) %>%
  arrange(english_name) %>%
  select(english_name, species_code, n_visits, total_count, RAI) %>%
  rename("Species" = english_name, "Code" = species_code,
         "Visits" = n_visits, "Total count" = total_count,
         "RAI (max count)" = RAI) %>%
  flextable() %>%
  colformat_double(j = "RAI (max count)", digits = 0) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption(sprintf("Site RAI: %s (%s). Maximum observed count across visits along %.1f km transect.",
                      SITE_NAME, SEASON_NAME, SITE_TRANSECT_KM))
save_word_table(ft1, "table1_site_RAI.docx", "Table 1: Site Relative Abundance Index")

# Table 2: Abundance comparison
ft2 <- species_comparison %>%
  arrange(desc(abs(log2(ratio_vs_regional)))) %>%
  select(english_name, species_code, RAI, regional_mean_RAI,
         national_mean_RAI, ratio_vs_regional, ratio_vs_national, pct_regional) %>%
  mutate(pct_regional = round(pct_regional * 100)) %>%
  rename("Species" = english_name, "Code" = species_code,
         "Site RAI (max)" = RAI,
         "BBS regional (pseudo-transect mean)" = regional_mean_RAI,
         "BBS national (pseudo-transect mean)" = national_mean_RAI,
         "Ratio site/regional" = ratio_vs_regional,
         "Ratio site/national" = ratio_vs_national,
         "Site percentile in region (%)" = pct_regional) %>%
  flextable() %>%
  colformat_double(j = 3:7, digits = 2) %>%
  colformat_double(j = 8, digits = 0) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption(paste0("Species abundance comparison: ", SITE_NAME, " (", SEASON_NAME,
                     "). BBS window ending ", year_bbs, ". ", effort_note)) %>%
  footnote(i = 1, j = 1, ref_symbol = "*",
           value = as_paragraph(paste("BBS source: Massimino et al. (2025).",
                                      "DOI: 10.5281/zenodo.18598453.",
                                      "Regional BBS:", radius_label)),
           part = "header")
save_word_table(ft2, "table2_abundance_comparison.docx", "Table 2: Abundance comparison")

# Table 3: Overlap
ft3 <- overlap_results %>%
  rename("Comparison" = Comparison, "Site species" = Site_species,
         "BBS species" = Reference_species, "Shared" = Shared,
         "Site only" = Site_only, "BBS only" = Reference_only,
         "Jaccard" = Jaccard, "Sorensen" = Sorensen) %>%
  flextable() %>%
  colformat_double(j = c("Jaccard","Sorensen"), digits = 3) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption("Species pool overlap: site vs BBS regional and national.") %>%
  footnote(i = 1, j = 1, ref_symbol = "*",
           value = as_paragraph(paste0(
             "Jaccard = shared/(shared+site-only+BBS-only). ",
             "S\u00f8rensen = 2\u00d7shared/(site+BBS totals). ",
             "Presence/absence only. Regional BBS: ", radius_label)),
           part = "header")
save_word_table(ft3, "table3_species_overlap.docx", "Table 3: Species pool overlap")

# Table 4: Richness vs same-effort pseudo-transects
ft4 <- richness_comparison %>%
  mutate(`95% range` = ifelse(is.na(CI_low), "\u2014",
                              sprintf("%.1f\u2013%.1f", CI_low, CI_high)),
         Site_percentile = ifelse(is.na(Site_percentile), "\u2014",
                                  paste0(Site_percentile, "%"))) %>%
  select(Dataset, Samples, Species_pool, Rarefied_S, `95% range`, Site_percentile) %>%
  rename("Samples (squares)" = Samples,
         "Species pool" = Species_pool,
         "Mean richness per sample" = Rarefied_S,
         "Site percentile" = Site_percentile) %>%
  flextable() %>%
  colformat_double(j = "Mean richness per sample", digits = 1, na_str = "\u2014") %>%
  colformat_double(j = "Samples (squares)", digits = 0) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption(paste0("Species richness per same-effort sample. ", effort_note,
                     " 95% range = middle 95% of squares."))
save_word_table(ft4, "table4_richness_pseudo_transects.docx", "Table 4: Species richness vs same-effort BBS samples")

# Table 5: Absent regional species
ft5 <- regional_absent %>%
  slice_max(bbs_RAI, n = 30) %>%
  mutate(occupancy_pct = round(occupancy * 100, 1)) %>%
  select(english_name, species_code, bbs_RAI, occupancy_pct) %>%
  rename("Species" = english_name, "Code" = species_code,
         "BBS regional RAI (pseudo-transect mean)" = bbs_RAI,
         "Regional occupancy (% pseudo-transects)" = occupancy_pct) %>%
  flextable() %>%
  colformat_double(j = "BBS regional RAI (pseudo-transect mean)", digits = 2) %>%
  colformat_double(j = "Regional occupancy (% pseudo-transects)", digits = 1) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption(paste0("Top 30 regional species absent from site, ranked by regional abundance. ",
                     radius_label))
save_word_table(ft5, "table5_absent_regional_species.docx",
                "Table 5: Regional species absent from site")

# Table 6: Site-only species
if (nrow(site_only_species) > 0) {
  ft6 <- site_only_species %>%
    select(english_name, species_code, RAI) %>%
    rename("Species" = english_name, "Code" = species_code, "Site RAI (max count)" = RAI) %>%
    flextable() %>%
    colformat_double(j = "Site RAI (max count)", digits = 0) %>%
    theme_booktabs() %>% autofit() %>%
    set_caption(paste0("Species recorded at site but absent from regional BBS (",
                       radius_used, " km radius). Likely winter visitors or species ",
                       "not covered by breeding season BBS."))
  save_word_table(ft6, "table6_site_only_species.docx",
                  "Table 6: Site-only species")
} else {
  message("No site-only species — table 6 skipped.")
}

# Table 7: Regional occupancy
ft7 <- regional_occupancy %>%
  left_join(bbs_rai %>% select(species_code, bbs_RAI), by = "species_code") %>%
  arrange(desc(occupancy)) %>%
  mutate(occupancy_pct    = round(occupancy * 100, 1),
         recorded_at_site = if_else(recorded_at_site, "Yes", "No")) %>%
  select(english_name, species_code, recorded_at_site, occupancy_pct, bbs_RAI) %>%
  rename("Species" = english_name, "Code" = species_code,
         "Recorded at site" = recorded_at_site,
         "Regional occupancy (% pseudo-transects)" = occupancy_pct,
         "BBS regional RAI (pseudo-transect mean)" = bbs_RAI) %>%
  flextable() %>%
  colformat_double(j = "Regional occupancy (% pseudo-transects)", digits = 1) %>%
  colformat_double(j = "BBS regional RAI (pseudo-transect mean)", digits = 2) %>%
  theme_booktabs() %>% autofit() %>%
  set_caption(paste0("Regional occupancy for all BBS species. ", radius_label))
save_word_table(ft7, "table7_regional_occupancy.docx",
                "Table 7: Regional species occupancy")

message("All tables saved.")

# =============================================================================
# 17. SUMMARY REPORT
# =============================================================================

top_over  <- divergence %>% filter(diff_regional > 0) %>%
  slice_max(diff_regional, n = 5) %>% pull(english_name)
top_under <- divergence %>% filter(diff_regional < 0) %>%
  slice_min(diff_regional, n = 5) %>% pull(english_name)
top_absent <- regional_absent %>% slice_max(bbs_RAI, n = 5) %>% pull(english_name)
top_site_only <- if (nrow(site_only_species) > 0) site_only_species$english_name else "None"

make_heading <- function(text, level = 1) {
  fpar(ftext(text, fp_text(
    bold = TRUE,
    font.size    = if (level == 1) 14 else 12,
    font.family  = "Arial",
    color        = if (level == 1) "#2E75B6" else "#404040"
  )))
}
make_body <- function(text) {
  fpar(ftext(text, fp_text(font.size = 11, font.family = "Arial")))
}

doc <- read_docx() %>%
  
  body_add_fpar(fpar(ftext(
    paste0("Bird Community Analysis \u2014 ", SITE_NAME),
    fp_text(bold = TRUE, font.size = 16, font.family = "Arial", color = "#2E75B6")
  ))) %>%
  body_add_fpar(make_body(paste0(
    "Grid reference: ", SITE_GRID_REF,
    "  |  Season: ", SEASON_NAME,
    "  |  BBS window ending: ", year_bbs,
    "  |  Regional reference: ", radius_label
  ))) %>%
  body_add_fpar(make_body(effort_note)) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_fpar(make_heading("1. Species Richness")) %>%
  body_add_fpar(make_body(paste0(
    SITE_NAME, " recorded ", site_total_species, " species from ",
    site_total_indivs, " individual detections along a ", SITE_TRANSECT_KM,
    " km transect. To compare like with like, each regional BBS square was treated as a ",
    "pseudo-transect of the same length, built from ", round(N_YEARS_EQUIV, 1),
    " years of surveys at that one square rather than by pooling several squares. ",
    "Across the ", n_squares_regional, " regional squares (", radius_label, "), ",
    "a pseudo-transect recorded on average ", round(regional_rare$S_expected, 1),
    " species (95% of squares: ", round(regional_rare$CI_low, 1),
    "\u2013", round(regional_rare$CI_high, 1), "); the site recorded more species than ",
    round(regional_rare$pct_site * 100), "% of them. ",
    "Nationally (", n_squares_national, " squares) the average was ",
    round(national_rare$S_expected, 1), " species (95% of squares: ",
    round(national_rare$CI_low, 1), "\u2013", round(national_rare$CI_high, 1),
    "), and the site ranked above ", round(national_rare$pct_site * 100), "% of squares."
  ))) %>%
  body_add_fpar(make_body(paste0(
    "Interpretation note: transect length is matched, but visit structure differs. ",
    "The site was surveyed on ", n_distinct(data_full$visit), " visits in one season, whereas each ",
    "BBS pseudo-transect comprises two visits per year over ", round(N_YEARS_EQUIV, 1), " years. ",
    "Year-to-year turnover can add occasional species to BBS pseudo-transects, so their richness ",
    "is a slightly generous benchmark; it is, however, far less inflated than pooling several ",
    "squares with different habitats."
  ))) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_fpar(make_heading("2. Community Composition")) %>%
  body_add_fpar(make_body(paste0(
    "Bray-Curtis similarity between the site assemblage and the regional BBS assemblage ",
    "was ", round(bc_similarity, 3),
    " (1.0 = identical proportional abundance structure; 0.0 = no shared abundance), ",
    "more typical than ", round(bc_site_pct * 100), "% of individual regional BBS squares. ",
    "Jaccard overlap with the regional BBS was ", overlap_results$Jaccard[1],
    " and with the national BBS was ", overlap_results$Jaccard[2], ". ",
    "Sørensen indices were ", overlap_results$Sorensen[1],
    " (regional) and ", overlap_results$Sorensen[2], " (national)."
  ))) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_fpar(make_heading("3. Species Divergence from Regional Assemblage")) %>%
  body_add_fpar(make_body(paste0(
    "The five species most over-represented at the site relative to their regional ",
    "proportional abundance were: ", paste(top_over, collapse = ", "), ". ",
    "These species characterise the site assemblage and likely reflect its dominant habitat types."
  ))) %>%
  body_add_fpar(make_body(paste0(
    "The five species most under-represented relative to the regional assemblage were: ",
    paste(top_under, collapse = ", "), ". ",
    "These absences or low abundances may reflect habitat unsuitability, ",
    "seasonal absence, or local population factors."
  ))) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_fpar(make_heading("4. Species Pool Gaps")) %>%
  body_add_fpar(make_body(paste0(
    overlap_results$Shared[1], " of the site's ", site_total_species,
    " species were also recorded in the regional BBS pool. ",
    nrow(regional_absent), " regional species were not recorded at the site. ",
    "The five most regionally abundant species absent from the site were: ",
    paste(top_absent, collapse = ", "), "."
  ))) %>%
  body_add_fpar(make_body(
    if (nrow(site_only_species) > 0) {
      paste0(nrow(site_only_species),
             " species were recorded at the site but not in the regional BBS pool within ",
             radius_used, " km: ", paste(top_site_only, collapse = ", "), ". ",
             "These may be scarce locally, poorly detected by BBS (e.g. nocturnal species) or passage migrants.")
    } else {
      "All site species were present in the regional BBS pool."
    }
  )) %>%
  body_add_fpar(make_body(
    if (length(site_not_in_bbs) > 0) {
      nm_of <- function(x) RAI_site_all$english_name[match(x, RAI_site_all$species_code)]
      other <- setdiff(site_not_in_bbs, assumed_present)
      paste0(
        if (length(assumed_present) > 0)
          paste0(paste(nm_of(assumed_present), collapse = ", "),
                 " is withheld from the published BBS dataset as a sensitive species. It is assumed to be ",
                 "present regionally and nationally, and counted as shared in the species-pool comparison, ",
                 "but excluded from abundance, richness and composition comparisons because BBS holds no counts for it. ")
        else "",
        if (length(other) > 0)
          paste0("Also recorded but absent from the published BBS dataset, so excluded from all comparisons: ",
                 paste(nm_of(other), collapse = ", "), ".")
        else "")
    } else ""
  )) %>%
  body_add_par("", style = "Normal") %>%
  
  body_add_fpar(make_heading("5. Methods and Data Sources")) %>%
  body_add_fpar(make_body(paste0(
    "Site survey: fixed-effort transect (", SITE_TRANSECT_KM, " km), ",
    "BTO two-letter species codes, ", SEASON_NAME, ". ",
    "Site RAI is calculated as the maximum observed count across survey visits per species. ",
    "BBS reference: Massimino et al. (2025), Zenodo DOI: 10.5281/zenodo.18598453. ",
    "BBS counts processed using max(E, L) visits per square and year",
    if (INCLUDE_DISTANCE_BAND_3) ", all distance bands included. " else ", excluding distance band 3. ",
    "Effort matching: each BBS square was treated as a pseudo-transect of ", SITE_TRANSECT_KM,
    " km by summing its ", K_FULL, " most recent surveyed years",
    if (K_FRAC > 0) sprintf(" plus the next most recent year at %.0f%% weight", K_FRAC * 100) else "",
    " (within ", LOOKBACK_MAX_YEARS, " years of ", year_bbs, "; unsurveyed years skipped). ",
    "Squares with fewer than ", K_NEEDED, " surveyed years were excluded. ",
    "BBS abundance benchmarks are the mean pseudo-transect count; site percentiles give the ",
    "share of regional pseudo-transects with a lower count (ties split). ",
    "Richness benchmarks are the distribution of species per pseudo-transect. ",
    "Community similarity: Bray-Curtis (vegan::vegdist()). ",
    "Overlap: Jaccard and S\u00f8rensen indices on presence/absence data."
  ))) %>%
  body_add_fpar(make_body(paste0(
    "Citation: Massimino D et al. (2025). Global Ecology and Biogeography 34(1): e13943. ",
    "https://doi.org/10.1111/geb.13943"
  )))

print(doc, target = file.path(out_dir, "summary_report.docx"))
message("Summary report saved.")

# =============================================================================
# END
# =============================================================================
message("\nAnalysis complete. All outputs saved to: ", out_dir)