rm(list=ls())

# =============================================================================
# Winter Bird Community Analysis — farm vs regional benchmark from eBird data
#
# Farm data  = winter survey counts from a GeoPackage and/or CSV (BTO codes),
#              as in the breeding-season tool.
# Benchmark  = your own winter farmland surveys entered in eBird ("Download My
#              Data" export).
#
# Effort varies between eBird checklists, so every benchmark visit is adjusted
# to the length of the farm's transect (STANDARD_KM = SITE_TRANSECT_KM):
#   counts     -> birds per STANDARD_KM (count x STANDARD_KM / km walked)
#   richness   -> species per STANDARD_KM, using a species-distance curve
#                 S = c x km^z fitted to the checklists
#   encounter probability -> chance a STANDARD_KM visit records the species, from an
#                 encounter-rate model fitted to presence/absence
# Counts entered as "X" are treated as present (richness, detection) and
# excluded from anything that needs a number.
#
# Unit of comparison: the winter survey = the average per visit, at one site, in
# one winter. The farm's focal winter is ranked against all benchmark
# winter surveys.
#
# Where the birds were on the farm (by crop, habitat or woodland type) is a
# separate script: 8_Winter_Crop_Analysis.R
# =============================================================================

library(tidyverse)
library(scales)
library(data.table)
library(officer)
library(flextable)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# --- Files ---
EBIRD_FILE     <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/ebird_farm_data.csv'
SPECIES_LOOKUP <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv"
GUILD_FILE     <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/BBS Community Analysis/species_guilds.csv"
out_dir        <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results/Winter/Community"

# --- The farm being assessed ---
SITE_NAME        <- "Shimpling Park Farm"
SITE_TRANSECT_KM <- 9.2            # length of the farm transect walked on each visit
FOCAL_WINTER     <- NA             # winter to assess, e.g. "2025/26"; NA = most recent in the farm data
# Farm survey data: GeoPackage and/or CSV (rows are combined). NA to skip either.
input_path <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'
if (!is.na(input_path)) { message("Layers in the GeoPackage:"); print(sf::st_layers(input_path)$name) }   # pick one of these
layer_name <- '2026_Counts_final'
csv_path   <- NA
CSV_COLS   <- c(species_code = "BTO_code",   # column names in the CSV
                transect     = "Site",
                count        = "Count",
                date_str     = "Date")
# If the farm's surveys are ALSO in the eBird export, list those eBird Location
# names / IDs here so the farm is not compared against itself.
FARM_EBIRD_LOCATIONS <- c("Shimpling Park Farm","USSAF Lavenham Airfield")

# --- Which checklists to use ---
WINTER_MONTHS   <- c(11, 12, 1, 2) # months of the year to include (1 = January); applies to farm and benchmark
COUNTIES        <- c("Suffolk","Norfolk", "Essex", "Cambridgeshire")    # benchmark counties; NA = all counties
# Benchmark locations: NA = every location in COUNTIES, or a specific list of
# Location names / IDs to use only those.
BENCHMARK_LOCATIONS <- NA
EXCLUDE_LOCATIONS   <- NA          # Location names / IDs to leave out of the benchmark
BENCHMARK_FIRST_WINTER <- NA       # e.g. "2015/16" to ignore older winters; NA = all
PROTOCOL_PATTERN <- "Travel"       # keep protocols containing this text (travelling counts)
COMPLETE_ONLY    <- TRUE           # only checklists where all species were reported
# Night-time checklists: sunrise and sunset are worked out for each checklist's
# date and coordinates. A checklist is kept only if it starts no earlier than
# TWILIGHT_BUFFER_MIN before sunrise and ends no later than that after sunset.
EXCLUDE_NOCTURNAL   <- TRUE
TWILIGHT_BUFFER_MIN <- 30               # minutes of twilight allowed either side
LOCAL_TZ            <- "Europe/London"  # time zone of the times in the eBird export
# Walking speed (distance / duration). Checklists slower than MIN or faster than
# MAX are left out of the benchmark; NA = no limit. An implausible speed usually
# means a typo in the distance or duration. Checklists with no duration are kept.
# Figure C2 (RUN_CHECKS) shows the spread of speeds to help choose the limits.
MIN_KM_PER_HOUR <- 1
MAX_KM_PER_HOUR <- 3
EXCLUDE_CHECKLISTS <- NA  # individual eBird checklist IDs to leave out, e.g. c("S123456789")

# --- Data checks ---
# TRUE = after the analysis, print tables and draw figures for tracking down odd
# points on figure W2 (species against distance), and save winter_w2_points_check.csv.
# Then use find_point() and inspect_visit() in the console (see the end of the script).
RUN_CHECKS     <- TRUE
CHECK_N_LABEL  <- 10      # how many of the most extreme visits to label on figure C1

# --- Effort standardisation ---
STANDARD_KM <- NA     # NA = the farm's transect length (SITE_TRANSECT_KM), so benchmark
# visits are adjusted to the same distance as the farm survey
MIN_KM      <- 0.5    # ignore checklists shorter than this ...
MAX_KM      <- 15     # ... or longer than this
MIN_VISITS  <- 1     # visits a site needs in a winter to count as a winter survey
# Shape of the species-distance curve used to adjust richness (see figure W2d,
# which compares all four on your data):
#   "power"     S = c x km^z            keeps rising; standard species-area form
#   "semilog"   S = a + b x log(km)     rises more slowly at long distances
#   "michaelis" S = Smax x km/(K + km)  levels off towards a maximum
#   "encounter" sum over species of the chance each is met in that distance
#               (built from the same encounter rates as the encounter probabilities)
RICHNESS_CURVE <- "power"
RICHNESS_Z  <- NA     # power curve only. NA = estimate the exponent from the data;
# or fix it (typical values 0.2-0.4)

# --- Species name matching ---
# eBird scientific name = name used in the BTO lookup, for the cases where they differ
SCI_SYNONYMS <- c("Curruca communis"     = "Sylvia communis",
                  "Curruca curruca"      = "Sylvia curruca",
                  "Curruca undata"       = "Sylvia undata",
                  "Acanthis flammea"     = "Acanthis cabaret",
                  "Astur gentilis"       = "Accipiter gentilis",
                  "Ardea ibis"           = "Bubulcus ibis",
                  "Botaurus minutus"     = "Ixobrychus minutus",
                  "Eudromias morinellus" = "Charadrius morinellus",
                  "Corvus monedula"      = "Coloeus monedula",
                  "Gulosus aristotelis"  = "Phalacrocorax aristotelis",
                  "Hydrobates leucorhous" = "Oceanodroma leucorhoa")

GUILD_ORDER <- c("Farmland", "Woodland", "Wetland", "Generalist")

# --- Figures ---
# Titles and subtitles are left off the images and written to
# figure_captions.docx / .csv instead. Set TRUE to print them on the figures.
SHOW_FIGURE_TITLES <- FALSE

if (is.na(STANDARD_KM)) STANDARD_KM <- SITE_TRANSECT_KM
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helpers ---
figure_captions <- tibble(file = character(), title = character(), notes = character())

save_fig <- function(plot, filename, width = 18, height = 12, dpi = 300) {
  lb <- tryCatch(ggplot2::get_labs(plot), error = function(e) plot$labels)
  one_line <- function(x, sep = " ") {
    if (is.null(x) || !is.character(x)) return("")
    str_squish(gsub("\n", sep, x))
  }
  figure_captions <<- bind_rows(
    figure_captions[figure_captions$file != filename, ],
    tibble(file = filename, title = one_line(lb$title, " — "), notes = one_line(lb$subtitle)))
  if (!SHOW_FIGURE_TITLES) plot <- plot + labs(title = NULL, subtitle = NULL)
  ggsave(file.path(out_dir, filename), plot = plot,
         width = width, height = height, units = "cm", dpi = dpi, bg = "white")
}

# Borderless labels that work in old and new ggplot2 (label.size was replaced by linewidth)
.new_label <- "border.colour" %in% names(formals(ggplot2::geom_label))
geom_label_nb <- function(...) if (.new_label) geom_label(..., linewidth = 0) else geom_label(..., label.size = 0)
annotate_label_nb <- function(...) if (.new_label) annotate("label", ..., linewidth = 0) else annotate("label", ..., label.size = 0)

# Sunrise and sunset (UTC) for a date and position; standard sunrise equation
sun_times <- function(date, lat, lon) {
  rad <- pi / 180
  n   <- as.numeric(date) + 2440588 - 2451545 + 0.0008
  Js  <- n - lon / 360
  M   <- (357.5291 + 0.98560028 * Js) %% 360
  C   <- 1.9148 * sin(M * rad) + 0.02 * sin(2 * M * rad) + 0.0003 * sin(3 * M * rad)
  L   <- (M + C + 180 + 102.9372) %% 360
  Jt  <- 2451545 + Js + 0.0053 * sin(M * rad) - 0.0069 * sin(2 * L * rad)
  dec <- asin(sin(L * rad) * sin(23.44 * rad))
  cw  <- (sin(-0.833 * rad) - sin(lat * rad) * sin(dec)) / (cos(lat * rad) * cos(dec))
  w   <- acos(pmin(pmax(cw, -1), 1)) / rad
  toP <- function(J) as.POSIXct((J - 2440587.5) * 86400, origin = "1970-01-01", tz = "UTC")
  list(sunrise = toP(Jt - w / 360), sunset = toP(Jt + w / 360))
}

# Encounter rate for one species: P(recorded on a visit of d km) = 1 - exp(-rate x d)
fit_rate <- function(present, km) {
  if (!any(present)) return(0)
  if (all(present)) return(Inf)
  nll <- function(lr) { r <- exp(lr); -sum(ifelse(present, log1p(-exp(-r * km)), -r * km)) }
  exp(optimize(nll, c(-12, 5))$minimum)
}

winter_label <- function(d) {
  y <- if_else(month(d) >= 7, year(d), year(d) - 1)
  sprintf("%d/%02d", y, (y + 1) %% 100)
}
winter_start <- function(w) as.integer(substr(w, 1, 4))
in_list <- function(location, loc_id, lst) location %in% lst | loc_id %in% lst
mid_rank_pct <- function(x, v) mean(x < v, na.rm = TRUE) + 0.5 * mean(x == v, na.rm = TRUE)
km_lab <- function() sprintf("%g km", STANDARD_KM)

# =============================================================================
# 1. LOAD THE EBIRD EXPORT
# =============================================================================

raw <- fread(file = EBIRD_FILE, colClasses = "character", encoding = "UTF-8", na.strings = "")
needed <- c("Submission ID", "Common Name", "Scientific Name", "Count", "County", "Location ID",
            "Location", "Date", "Protocol", "All Obs Reported", "Distance Traveled (km)")
if (length(setdiff(needed, names(raw))) > 0)
  stop("eBird file is missing column(s): ", paste(setdiff(needed, names(raw)), collapse = ", "))

obs <- as_tibble(raw) %>%
  transmute(checklist = `Submission ID`,
            common    = `Common Name`,
            sci       = `Scientific Name`,
            count_raw = str_trim(Count),
            county    = County,
            loc_id    = `Location ID`,
            location  = Location,
            date      = as.Date(parse_date_time(Date, orders = c("Ymd", "dmY", "mdY"), quiet = TRUE)),
            protocol  = Protocol,
            complete  = `All Obs Reported` %in% c("1", "TRUE", "true"),
            km        = suppressWarnings(as.numeric(`Distance Traveled (km)`)),
            time_str  = if ("Time" %in% names(raw)) Time else NA_character_,
            duration  = if ("Duration (Min)" %in% names(raw)) suppressWarnings(as.numeric(`Duration (Min)`)) else NA_real_,
            lat       = if ("Latitude" %in% names(raw)) suppressWarnings(as.numeric(Latitude)) else NA_real_,
            lon       = if ("Longitude" %in% names(raw)) suppressWarnings(as.numeric(Longitude)) else NA_real_)

# =============================================================================
# 1b. LOAD THE FARM SURVEY DATA (GeoPackage and/or CSV, BTO codes)
# =============================================================================

parse_survey_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(parse_date_time(as.character(x), orders = c("Ymd", "dmY", "dmy", "Ymd HMS", "dmY HM"), quiet = TRUE))
}
farm_sources <- list()
if (!is.na(input_path)) {
  gpkg_layers <- sf::st_layers(input_path)$name
  if (!layer_name %in% gpkg_layers)
    stop("Layer '", layer_name, "' is not in the GeoPackage. Layers available: ", paste(gpkg_layers, collapse = ", "))
  farm_sources$gpkg <- sf::st_read(input_path, layer = layer_name, quiet = TRUE) %>%
    sf::st_drop_geometry() %>% as_tibble() %>%
    rename(species_code = BTO_code, transect = Site, count = Count, date_str = Date)
  message(sprintf("Farm GeoPackage: %d rows from layer '%s'.", nrow(farm_sources$gpkg), layer_name))
}
if (!is.na(csv_path)) {
  csv_raw <- read_csv(csv_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  if (length(setdiff(CSV_COLS, names(csv_raw))) > 0)
    stop("Farm CSV is missing column(s): ", paste(setdiff(CSV_COLS, names(csv_raw)), collapse = ", "),
         ". Columns found: ", paste(names(csv_raw), collapse = ", "), ". Update CSV_COLS.")
  farm_sources$csv <- csv_raw %>% rename(all_of(CSV_COLS))
  message(sprintf("Farm CSV: %d rows from %s.", nrow(farm_sources$csv), basename(csv_path)))
}
if (length(farm_sources) == 0) stop("No farm data: set input_path and/or csv_path.")

farm_obs <- farm_sources %>%
  map(~ transmute(.x, species_code = str_trim(as.character(species_code)),
                  date = parse_survey_date(date_str),
                  count = replace_na(suppressWarnings(as.numeric(count)), 0))) %>%
  bind_rows()
if (any(is.na(farm_obs$date)))
  warning(sum(is.na(farm_obs$date)), " farm rows have unreadable dates and are excluded.")
farm_obs <- farm_obs %>%
  filter(!is.na(date), month(date) %in% WINTER_MONTHS, count > 0) %>%
  mutate(winter = winter_label(date))
if (nrow(farm_obs) == 0)
  stop("No farm records fall in WINTER_MONTHS (", paste(month.abb[WINTER_MONTHS], collapse = ", "), ").")
if (is.na(FOCAL_WINTER)) FOCAL_WINTER <- farm_obs$winter[which.max(winter_start(farm_obs$winter))]
message(sprintf("Farm data: %d records on %d visit dates; winters: %s. Assessing %s.",
                nrow(farm_obs), n_distinct(farm_obs$date), paste(sort(unique(farm_obs$winter)), collapse = ", "), FOCAL_WINTER))

# =============================================================================
# 2. FILTER CHECKLISTS
# =============================================================================

checklists <- obs %>%
  distinct(checklist, .keep_all = TRUE) %>%
  select(checklist, county, loc_id, location, date, protocol, complete, km, time_str, duration, lat, lon) %>%
  mutate(is_focal = FALSE)   # the farm's own data comes from the GeoPackage/CSV

# --- Daylight check: start and end of each checklist against sunrise and sunset ---
checklists <- checklists %>%
  mutate(start_local = parse_date_time(paste(format(date), time_str), orders = c("Ymd IMp", "Ymd HM", "Ymd HMS"),
                                       tz = LOCAL_TZ, quiet = TRUE),
         end_local   = start_local + minutes(round(coalesce(duration, 0))),
         sunrise     = with_tz(sun_times(date, lat, lon)$sunrise, LOCAL_TZ),
         sunset      = with_tz(sun_times(date, lat, lon)$sunset, LOCAL_TZ),
         light_known = !is.na(start_local) & !is.na(lat) & !is.na(lon),
         nocturnal   = light_known & (start_local < sunrise - minutes(TWILIGHT_BUFFER_MIN) |
                                        end_local   > sunset  + minutes(TWILIGHT_BUFFER_MIN)),
         km_per_hour = if_else(!is.na(duration) & duration > 0 & !is.na(km), km / (duration / 60), NA_real_),
         speed_out   = !is.na(km_per_hour) &
           ((!is.na(MIN_KM_PER_HOUR) & km_per_hour < MIN_KM_PER_HOUR) |
              (!is.na(MAX_KM_PER_HOUR) & km_per_hour > MAX_KM_PER_HOUR)))
speed_filter_on <- !is.na(MIN_KM_PER_HOUR) | !is.na(MAX_KM_PER_HOUR)
speed_label <- if (!speed_filter_on) "Walking speed filter off" else
  paste0("Walking speed ", if (is.na(MIN_KM_PER_HOUR)) "up to " else paste0(MIN_KM_PER_HOUR, if (is.na(MAX_KM_PER_HOUR)) " km/h or more" else "-"),
         if (is.na(MAX_KM_PER_HOUR)) "" else paste0(MAX_KM_PER_HOUR, " km/h"))

filter_log <- tibble(step = "All checklists in export", checklists = nrow(checklists))
keep_step <- function(df, cond, label) {
  out <- df[cond & !is.na(cond), ]
  filter_log <<- bind_rows(filter_log, tibble(step = label, checklists = nrow(out)))
  out
}

checklists <- checklists %>%
  keep_step(month(.$date) %in% WINTER_MONTHS,
            paste0("Months: ", paste(month.abb[WINTER_MONTHS], collapse = ", "))) %>%
  keep_step(all(is.na(COUNTIES)) | .$county %in% COUNTIES,
            paste0("Counties: ", if (all(is.na(COUNTIES))) "all" else paste(COUNTIES, collapse = ", "))) %>%
  keep_step(all(is.na(FARM_EBIRD_LOCATIONS)) | !in_list(.$location, .$loc_id, FARM_EBIRD_LOCATIONS),
            "Farm's own eBird locations removed") %>%
  keep_step(.$is_focal | all(is.na(BENCHMARK_LOCATIONS)) | in_list(.$location, .$loc_id, BENCHMARK_LOCATIONS),
            "Benchmark location list") %>%
  keep_step(.$is_focal | all(is.na(EXCLUDE_LOCATIONS)) | !in_list(.$location, .$loc_id, EXCLUDE_LOCATIONS),
            "Excluded locations removed") %>%
  keep_step(grepl(PROTOCOL_PATTERN, .$protocol, ignore.case = TRUE),
            paste0("Protocol contains '", PROTOCOL_PATTERN, "'")) %>%
  keep_step(!COMPLETE_ONLY | .$complete, "Complete checklists only") %>%
  keep_step(!is.na(.$km) & .$km >= MIN_KM & .$km <= MAX_KM,
            sprintf("Distance recorded, %g-%g km", MIN_KM, MAX_KM)) %>%
  { nocturnal_removed <<- if (EXCLUDE_NOCTURNAL) filter(., nocturnal) else filter(., FALSE)
  n_no_time <<- sum(!.$light_known); . } %>%
  keep_step(!EXCLUDE_NOCTURNAL | !.$nocturnal,
            if (EXCLUDE_NOCTURNAL) sprintf("Daytime only (within %g min of sunrise and sunset)", TWILIGHT_BUFFER_MIN)
            else "Daytime filter off") %>%
  { speed_before <<- .; . } %>%
  keep_step(!.$speed_out, speed_label) %>%
  keep_step(all(is.na(EXCLUDE_CHECKLISTS)) | !.$checklist %in% EXCLUDE_CHECKLISTS, "Individually excluded checklists removed") %>%
  mutate(winter = winter_label(date),
         site   = location)

if (!is.na(BENCHMARK_FIRST_WINTER))
  checklists <- checklists %>%
  keep_step(.$is_focal | winter_start(.$winter) >= winter_start(BENCHMARK_FIRST_WINTER),
            paste0("Benchmark winters from ", BENCHMARK_FIRST_WINTER))

if (EXCLUDE_NOCTURNAL) {
  if (nrow(nocturnal_removed) > 0) {
    write_csv(nocturnal_removed %>%
                transmute(checklist, location, date, start = format(start_local, "%H:%M"), end = format(end_local, "%H:%M"),
                          sunrise = format(sunrise, "%H:%M"), sunset = format(sunset, "%H:%M"), km) %>% arrange(date, start),
              file.path(out_dir, "winter_nocturnal_checklists_removed.csv"))
    message(nrow(nocturnal_removed), " night-time checklists removed (listed in winter_nocturnal_checklists_removed.csv).")
  }
  if (n_no_time > 0)
    message(n_no_time, " checklists have no start time or coordinates, so could not be checked for daylight; they are kept.")
}
speed_removed <- speed_before %>% filter(speed_out)
if (nrow(speed_removed) > 0) {
  write_csv(speed_removed %>% transmute(checklist, location, date, start = time_str, minutes = duration, km,
                                        km_per_hour = round(km_per_hour, 2)) %>% arrange(km_per_hour),
            file.path(out_dir, "winter_speed_checklists_removed.csv"))
  message(nrow(speed_removed), " checklists removed for walking speed (listed in winter_speed_checklists_removed.csv).")
}
message("Checklist filtering:")
print(as.data.frame(filter_log), row.names = FALSE)

# =============================================================================
# 3. SPECIES CLEANING AND NAME MATCHING
# =============================================================================
# Drops entries that are not a single wild species ("gull sp.", "A/B" pairs,
# hybrids, domestic types). Subspecies groups are merged into their species
# using the first two words of the scientific name.

sp_obs <- obs %>%
  semi_join(checklists, by = "checklist") %>%
  filter(!is.na(sci)) %>%
  # Judged on the scientific name: "Genus sp.", "A/B" pairs and "A x B" hybrids.
  # Common names are not used for this, because some true species have a slash
  # in their eBird name (e.g. "Pied Wagtail/White Wagtail" = Motacilla alba).
  mutate(not_species = grepl(" sp\\.|/| x |\\(hybrid\\)|\\(Domestic", sci))
dropped_taxa <- sp_obs %>% filter(not_species) %>% count(common, name = "records") %>% arrange(desc(records))
if (nrow(dropped_taxa) > 0)
  message("Entries dropped (not a single wild species): ",
          paste(sprintf("%s (%d)", dropped_taxa$common, dropped_taxa$records), collapse = "; "))

lookup <- read_csv(SPECIES_LOOKUP, col_types = cols(.default = "c")) %>%
  transmute(sci_name = str_trim(sci_name), species_code = str_trim(bto_code),
            english_name = str_trim(species), BOCC5 = str_to_title(str_trim(bocc))) %>%
  arrange(species_code != "FP") %>%                 # Columba livia on farmland = Feral Pigeon
  distinct(sci_name, .keep_all = TRUE)

sp_obs <- sp_obs %>%
  filter(!not_species) %>%
  mutate(binomial = word(sci, 1, 2),
         sci_key  = coalesce(unname(SCI_SYNONYMS[binomial]), binomial),
         ebird_name = str_trim(gsub("\\s*\\(.*\\)", "", common))) %>%
  left_join(lookup, by = c("sci_key" = "sci_name")) %>%
  mutate(english_name = coalesce(english_name, ebird_name),
         species_code = coalesce(species_code, sci_key),
         BOCC5        = if_else(BOCC5 %in% c("Red", "Amber", "Green", "Non-Native"), BOCC5, "Not listed"),
         BOCC5        = recode(BOCC5, "Non-Native" = "Non-native"),
         is_x         = toupper(count_raw) == "X",
         count        = suppressWarnings(as.numeric(count_raw)))

unmatched <- sp_obs %>% filter(species_code == sci_key) %>% distinct(english_name, sci_key)
if (nrow(unmatched) > 0)
  message("Species not in the BTO lookup (eBird names used, BOCC5 'Not listed'; add to SCI_SYNONYMS if a name differs): ",
          paste(sprintf("%s [%s]", unmatched$english_name, unmatched$sci_key), collapse = "; "))

# Farm species (BTO codes) get their names and status from the same lookup
lookup_by_code <- read_csv(SPECIES_LOOKUP, col_types = cols(.default = "c")) %>%
  transmute(species_code = str_trim(bto_code), english_name = str_trim(species),
            BOCC5 = str_to_title(str_trim(bocc))) %>%
  mutate(BOCC5 = recode(BOCC5, "Non-Native" = "Non-native"),
         BOCC5 = if_else(BOCC5 %in% c("Red", "Amber", "Green", "Non-native"), BOCC5, "Not listed")) %>%
  distinct(species_code, .keep_all = TRUE)
farm_unknown <- setdiff(unique(farm_obs$species_code), lookup_by_code$species_code)
if (length(farm_unknown) > 0)
  message("Farm species codes not in the BTO lookup (shown by code): ", paste(farm_unknown, collapse = ", "))

species_info <- bind_rows(
  sp_obs %>% distinct(species_code, english_name, BOCC5),
  farm_obs %>% distinct(species_code) %>% left_join(lookup_by_code, by = "species_code") %>%
    mutate(english_name = coalesce(english_name, species_code), BOCC5 = coalesce(BOCC5, "Not listed"))) %>%
  distinct(species_code, .keep_all = TRUE)

if (!is.na(GUILD_FILE) && file.exists(GUILD_FILE)) {
  guilds <- read_csv(GUILD_FILE, col_types = cols(.default = "c")) %>%
    transmute(species_code = str_trim(bto_code), guild = na_if(str_to_title(str_trim(guild)), "")) %>%
    distinct(species_code, .keep_all = TRUE)
} else {
  warning("Guild file not found - all species grouped as 'Other'.")
  guilds <- tibble(species_code = character(), guild = character())
}
species_info <- species_info %>% left_join(guilds, by = "species_code") %>%
  mutate(guild = coalesce(guild, "Other"))
guild_levels <- c(intersect(GUILD_ORDER, species_info$guild),
                  sort(setdiff(species_info$guild, c(GUILD_ORDER, "Other"))), "Other")

# =============================================================================
# 4. VISITS (site x date) AND STANDARDISATION
# =============================================================================

visits <- bind_rows(
  checklists %>%
    group_by(site, is_focal, winter, date) %>%
    summarise(km = sum(km), n_checklists = n(), .groups = "drop"),
  farm_obs %>% distinct(date, winter) %>%
    mutate(site = SITE_NAME, is_focal = TRUE, km = SITE_TRANSECT_KM, n_checklists = 1L)) %>%
  mutate(visit_id = row_number())

visit_species <- bind_rows(
  sp_obs %>%
    inner_join(checklists %>% select(checklist, site), by = "checklist") %>%
    inner_join(visits %>% filter(!is_focal) %>% select(visit_id, site, date), by = c("site", "date")) %>%
    group_by(visit_id, species_code) %>%
    summarise(count = sum(count), has_x = any(is_x), .groups = "drop"),   # NA if any X
  farm_obs %>%
    inner_join(visits %>% filter(is_focal) %>% select(visit_id, date), by = "date") %>%
    group_by(visit_id, species_code) %>%
    summarise(count = sum(count), has_x = FALSE, .groups = "drop"))

# Full grid with zeros: every visit x every species in the pool
grid <- expand_grid(visit_id = visits$visit_id, species_code = unique(visit_species$species_code)) %>%
  left_join(visit_species, by = c("visit_id", "species_code")) %>%
  mutate(present = !is.na(has_x),
         count   = if_else(present, count, 0),          # NA stays NA where count was X
         has_x   = coalesce(has_x, FALSE)) %>%
  left_join(visits %>% select(visit_id, km), by = "visit_id") %>%
  mutate(count_std = count * STANDARD_KM / km)

visits <- visits %>%
  left_join(grid %>% group_by(visit_id) %>%
              summarise(S = sum(present), any_x = any(has_x),
                        total_std = if (any(has_x)) NA_real_ else sum(count_std), .groups = "drop"),
            by = "visit_id")

# --- Species-distance curve: S = c x km^z ---
fit_data <- visits %>% filter(!is_focal, S > 0)   # curve is fitted to the benchmark only
fit_lm <- lm(log(S) ~ log(km), data = fit_data)    # straight line on log-log axes
z_ci   <- confint(fit_lm)[2, ]
fit_r2 <- summary(fit_lm)$r.squared
z_fit <- if (nrow(fit_data) >= 20 && max(fit_data$km) / min(fit_data$km) >= 1.5)
  unname(coef(fit_lm)[2]) else NA_real_
z_used <- if (!is.na(RICHNESS_Z)) RICHNESS_Z else if (is.na(z_fit)) 0.3 else min(max(z_fit, 0), 1)
z_source <- if (!is.na(RICHNESS_Z)) "fixed in the configuration" else
  if (is.na(z_fit)) "default (too few visits or too little variation in distance to estimate)" else
    "estimated from the checklists"
c_fit <- exp(mean(log(fit_data$S) - z_used * log(fit_data$km)))
message(sprintf("Species-distance exponent z = %.2f (%s). Doubling the distance adds about %.0f%% more species.",
                z_used, z_source, (2^z_used - 1) * 100))
message(sprintf("Curve fit: z = %.2f (95%% CI %.2f to %.2f), R-squared = %.2f, %d benchmark visits of %.1f-%.1f km. See figures W2 and W2b.",
                coef(fit_lm)[2], z_ci[1], z_ci[2], fit_r2, nrow(fit_data), min(fit_data$km), max(fit_data$km)))

# --- Candidate curves, all fitted to the benchmark visits ---
curve_fns <- list(power = function(d) c_fit * d^z_used)
fit_sl <- lm(S ~ log(km), data = fit_data)
curve_fns$semilog <- function(d) pmax(coef(fit_sl)[1] + coef(fit_sl)[2] * log(d), 1)
fit_mm <- tryCatch(nls(S ~ Smax * km / (K + km), data = fit_data,
                       start = list(Smax = max(fit_data$S) * 1.2, K = median(fit_data$km) / 2)),
                   error = function(e) NULL)
if (!is.null(fit_mm)) curve_fns$michaelis <- function(d) coef(fit_mm)[["Smax"]] * d / (coef(fit_mm)[["K"]] + d)
rates_curve <- grid %>% filter(visit_id %in% fit_data$visit_id) %>%
  group_by(species_code) %>% summarise(rate = fit_rate(present, km), .groups = "drop")
curve_fns$encounter <- function(d) vapply(d, function(x) sum(1 - exp(-rates_curve$rate * x)), numeric(1))

curve_labels <- c(power = "Power: S = c x km^z", semilog = "Semi-log: S = a + b x log(km)",
                  michaelis = "Michaelis-Menten: levels off", encounter = "Encounter model: sum of species' chances")
curve_stats <- tibble(model = names(curve_fns)) %>%
  mutate(label    = curve_labels[model],
         rmse     = map_dbl(model, ~ sqrt(mean((fit_data$S - curve_fns[[.x]](fit_data$km))^2))),
         expected = map_dbl(model, ~ curve_fns[[.x]](STANDARD_KM)))
if (!RICHNESS_CURVE %in% names(curve_fns)) {
  warning("RICHNESS_CURVE '", RICHNESS_CURVE, "' is not available (", paste(names(curve_fns), collapse = ", "), "). Using 'power'.")
  RICHNESS_CURVE <- "power"
}
curve_f    <- curve_fns[[RICHNESS_CURVE]]
curve_desc <- switch(RICHNESS_CURVE,
                     power     = sprintf("power curve S = c x km^z, z = %.2f, %s", z_used, z_source),
                     semilog   = "semi-log curve S = a + b x log(km)",
                     michaelis = sprintf("Michaelis-Menten curve levelling off at %.0f species", coef(fit_mm)[["Smax"]]),
                     encounter = "encounter model: the sum over species of the chance each is recorded in that distance")
message("Species-distance curves (typical error in species per visit; species expected at ", km_lab(), "):")
print(as.data.frame(curve_stats %>% transmute(model, typical_error = round(rmse, 2), expected_species = round(expected, 1),
                                              used = if_else(model == RICHNESS_CURVE, "<- used", ""))), row.names = FALSE)

visits <- visits %>% mutate(S_std = S * curve_f(STANDARD_KM) / curve_f(km))
expected_S <- curve_f(STANDARD_KM)   # species expected on a regional visit of the standard distance

# =============================================================================
# 5. SITE-WINTERS: THE FARM AND THE BENCHMARK
# =============================================================================

site_winters <- visits %>%
  group_by(site, is_focal, winter) %>%
  summarise(n_visits = n(), mean_km = mean(km),
            S_std = mean(S_std), total_std = mean(total_std, na.rm = TRUE),
            n_count_visits = sum(!any_x), .groups = "drop") %>%
  mutate(total_std = if_else(n_count_visits == 0, NA_real_, total_std))

focal_sw <- site_winters %>% filter(is_focal, winter == FOCAL_WINTER)
if (nrow(focal_sw) == 0)
  stop("No farm visits in winter ", FOCAL_WINTER, ". Winters in the farm data: ",
       paste(sort(unique(site_winters$winter[site_winters$is_focal])), collapse = ", "))

bench_sw <- site_winters %>% filter(!is_focal, n_visits >= MIN_VISITS)
if (nrow(bench_sw) < 10)
  warning("Only ", nrow(bench_sw), " benchmark winter surveys - percentiles will be coarse. ",
          "Consider widening COUNTIES / WINTER_MONTHS or lowering MIN_VISITS.")

n_sites   <- n_distinct(bench_sw$site)
n_sw      <- nrow(bench_sw)
bench_yrs <- range(winter_start(bench_sw$winter))
bench_label <- sprintf("%d winter surveys from %d farms%s, winters %d/%02d–%d/%02d",
                       n_sw, n_sites,
                       if (all(is.na(COUNTIES))) "" else paste0(" in ", paste(COUNTIES, collapse = ", ")),
                       bench_yrs[1], (bench_yrs[1] + 1) %% 100, bench_yrs[2], (bench_yrs[2] + 1) %% 100)
std_note <- sprintf("All values are averages per visit; benchmark visits are adjusted to %s, %s.", km_lab(),
                    if (isTRUE(all.equal(STANDARD_KM, SITE_TRANSECT_KM))) "the length of the farm transect" else "the standard distance")
message("Benchmark: ", bench_label)
message(sprintf("Farm: %d visits in %s, %.1f km per visit. Regional expectation at %s: %.1f species per visit.",
                focal_sw$n_visits, FOCAL_WINTER, focal_sw$mean_km, km_lab(), expected_S))

focal_visits <- visits %>% filter(is_focal, winter == FOCAL_WINTER)
bench_visits <- visits %>% semi_join(bench_sw, by = c("site", "winter"))

# The standardisation is most reliable when the farm's visits and STANDARD_KM
# sit inside the range of distances walked on benchmark visits.
bench_km_range <- quantile(bench_visits$km, c(0.025, 0.975))
if (focal_sw$mean_km < bench_km_range[1] || focal_sw$mean_km > bench_km_range[2])
  warning(sprintf("The farm transect is %.1f km, outside the range of most benchmark visits (%.1f-%.1f km). The species-distance curve is being extrapolated to reach it, so treat the richness comparison with caution (see figure W2).",
                  focal_sw$mean_km, bench_km_range[1], bench_km_range[2]))
if (!isTRUE(all.equal(STANDARD_KM, SITE_TRANSECT_KM)) && (STANDARD_KM < bench_km_range[1] || STANDARD_KM > bench_km_range[2]))
  warning(sprintf("STANDARD_KM (%g) is outside the range of most benchmark visits (%.1f-%.1f km). Consider a value near the typical visit, e.g. %.1f km.",
                  STANDARD_KM, bench_km_range[1], bench_km_range[2], median(bench_visits$km)))

# Species means per winter survey (X counts ignored for that species)
sw_species <- grid %>%
  inner_join(visits %>% select(visit_id, site, is_focal, winter), by = "visit_id") %>%
  group_by(site, is_focal, winter, species_code) %>%
  summarise(mean_std = mean(count_std, na.rm = TRUE), peak = suppressWarnings(max(count, na.rm = TRUE)),
            visits_present = sum(present), n_visits = n(), .groups = "drop") %>%
  mutate(mean_std = if_else(is.nan(mean_std), NA_real_, mean_std),
         peak     = if_else(is.finite(peak), peak, NA_real_))

farm_sp  <- sw_species %>% filter(is_focal, winter == FOCAL_WINTER)
bench_sp <- sw_species %>% semi_join(bench_sw, by = c("site", "winter"))
site_present <- farm_sp %>% filter(visits_present > 0) %>% pull(species_code)

# =============================================================================
# 6. ENCOUNTER PROBABILITY AT THE STANDARD DISTANCE
# =============================================================================
# Encounter-rate model: P(recorded on a visit of d km) = 1 - exp(-rate x d).
# One rate per species, fitted to presence/absence on all benchmark visits
# (X counts included, since presence is all that is needed).

detection <- grid %>%
  semi_join(bench_visits, by = "visit_id") %>%
  group_by(species_code) %>%
  summarise(rate = fit_rate(present, km), bench_visits_present = sum(present), .groups = "drop") %>%
  mutate(p_std = 1 - exp(-rate * STANDARD_KM))

# =============================================================================
# 7. SPECIES COMPARISON TABLE
# =============================================================================

bench_summary <- bench_sp %>%
  group_by(species_code) %>%
  summarise(reg_mean = mean(mean_std, na.rm = TRUE), reg_median = median(mean_std, na.rm = TRUE),
            reg_q25 = quantile(mean_std, 0.25, na.rm = TRUE), reg_q75 = quantile(mean_std, 0.75, na.rm = TRUE),
            sites_recorded = n_distinct(site[visits_present > 0]), .groups = "drop")

species_comparison <- farm_sp %>%
  select(species_code, farm_mean = mean_std, farm_peak = peak, visits_present, n_visits) %>%
  full_join(bench_summary, by = "species_code") %>%
  left_join(bench_sp %>% select(species_code, b = mean_std) %>% nest(b = b), by = "species_code") %>%
  mutate(pct = map2_dbl(b, farm_mean, ~ if (is.na(.y) || is.null(.x)) NA_real_ else mid_rank_pct(.x$b, .y))) %>%
  select(-b) %>%
  left_join(detection %>% select(species_code, p_std), by = "species_code") %>%
  left_join(species_info, by = "species_code") %>%
  mutate(recorded = species_code %in% site_present,
         ratio    = farm_mean / reg_mean) %>%
  replace_na(list(p_std = 0, reg_mean = 0))

farm_wide <- species_comparison %>% filter(recorded, !is.na(farm_mean)) %>%
  mutate(farm_share = farm_mean / sum(farm_mean)) %>% select(species_code, farm_share)
species_comparison <- species_comparison %>%
  left_join(farm_wide, by = "species_code") %>%
  mutate(reg_share  = reg_mean / sum(reg_mean, na.rm = TRUE),
         farm_share = coalesce(farm_share, 0),
         diff_share = farm_share - reg_share,      # composition (used for the summary text and Bray-Curtis)
         diff_mean  = farm_mean - reg_mean)        # birds per standard visit: site minus regional mean

richness_pct <- mid_rank_pct(bench_sw$S_std, focal_sw$S_std)
total_pct    <- mid_rank_pct(bench_sw$total_std, focal_sw$total_std)
farm_total_species <- length(site_present)
message(sprintf("%s, %s: %d species over the winter; %.1f species and %.0f birds per %s visit (higher than %.0f%% and %.0f%% of benchmark winter surveys).",
                SITE_NAME, FOCAL_WINTER, farm_total_species, focal_sw$S_std, focal_sw$total_std, km_lab(),
                richness_pct * 100, total_pct * 100))

# =============================================================================
# 8. FIGURES
# =============================================================================

fig_title <- function(txt) txt

# --- W1: Farm species, mean count per standard distance ---
w1 <- species_comparison %>% filter(recorded, !is.na(farm_mean)) %>%
  mutate(english_name = fct_reorder(english_name, farm_mean)) %>%
  ggplot(aes(farm_mean, english_name)) + geom_col(fill = "steelblue") +
  labs(title = fig_title(paste0("Winter birds at ", SITE_NAME, ", ", FOCAL_WINTER)),
       subtitle = sprintf("Mean count per visit along the %g km transect (%d visits).",
                          SITE_TRANSECT_KM, focal_sw$n_visits),
       x = paste0("Mean birds per ", km_lab(), " visit"), y = NULL)
save_fig(w1, "figW1_farm_species.png", width = 18, height = max(10, farm_total_species * 0.42))

# --- W2: Species-distance curve, and the species expected at the farm's distance ---
curve <- tibble(km = seq(min(visits$km), max(c(visits$km, STANDARD_KM)), length.out = 100)) %>% mutate(S = curve_f(km))
w2 <- ggplot(bench_visits, aes(km, S)) +
  geom_point(colour = "grey60", alpha = 0.6, size = 1.6) +
  geom_line(data = curve, colour = "grey20", linewidth = 0.8) +
  geom_vline(xintercept = STANDARD_KM, linetype = "dashed", colour = "grey40") +
  annotate("point", x = STANDARD_KM, y = expected_S, shape = 21, size = 4, fill = "white", colour = "grey20", stroke = 1.1) +
  annotate_label_nb(x = STANDARD_KM, y = expected_S, hjust = 1.08, vjust = 1.4, size = 3.2,
                    fill = alpha("white", 0.8), colour = "grey20",
                    label = sprintf("Expected at %s:\n%.1f species", km_lab(), expected_S)) +
  geom_point(data = focal_visits, shape = 23, size = 4, fill = "steelblue", colour = "white") +
  labs(title = fig_title("Species recorded against distance walked"),
       subtitle = str_wrap(sprintf("Grey = benchmark visits; blue = %s visits in %s. Line = species-distance curve fitted to the benchmark (%s). Open circle = species expected on a regional visit of %s, the length of the farm transect.",
                                   SITE_NAME, FOCAL_WINTER, curve_desc, km_lab()), 95),
       x = "Distance walked (km)", y = "Species recorded on the visit")
save_fig(w2, "figW2_species_distance_curve.png", width = 18, height = 12)
print(w2)   # also shown in the R plot window

# --- W2b: How well the curve fits (log-log axes, where the curve is a straight line) ---
fit_band <- tibble(km = exp(seq(log(min(fit_data$km)), log(max(c(fit_data$km, STANDARD_KM))), length.out = 100)))
fit_band <- bind_cols(fit_band, as_tibble(exp(predict(fit_lm, fit_band, interval = "confidence"))))
fit_resid <- fit_data %>% mutate(resid = log(S / curve_f(km)),   # relative to the curve in use
                                 km_band = cut(km, breaks = unique(quantile(km, seq(0, 1, 0.2))), include.lowest = TRUE))
fit_text <- sprintf("z = %.2f (95%% CI %.2f to %.2f)\nR-squared = %.2f\n%d benchmark visits, %.1f to %.1f km%s",
                    coef(fit_lm)[2], z_ci[1], z_ci[2], fit_r2, nrow(fit_data), min(fit_data$km), max(fit_data$km),
                    if (STANDARD_KM > max(fit_data$km)) sprintf("\nFarm distance (%s) is beyond the data", km_lab()) else "")
w2b_top <- ggplot(fit_data, aes(km, S)) +
  geom_point(colour = "grey55", alpha = 0.6, size = 1.6) +
  geom_ribbon(data = fit_band, aes(x = km, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "steelblue", alpha = 0.25) +
  geom_line(data = fit_band, aes(km, fit), colour = "steelblue", linewidth = 0.9) +
  geom_vline(xintercept = STANDARD_KM, linetype = "dashed", colour = "grey40") +
  annotate_label_nb(x = min(fit_data$km), y = max(fit_data$S), hjust = 0, vjust = 1, size = 3.2,
                    fill = alpha("white", 0.8), colour = "grey20", label = fit_text) +
  scale_x_log10() + scale_y_log10() +
  labs(title = fig_title("Fit of the power curve"),
       subtitle = str_wrap(sprintf("Top: benchmark visits on log-log axes, where the curve S = c x km^z is a straight line. Blue = fitted line with 95%% confidence band; dashed = the farm's distance (%s). Bottom: how far visits fall from the line, grouped by distance. Boxes centred on zero at every distance mean the curve describes the data well; a drift up or down means it does not.",
                                   km_lab()), 95),
       x = "Distance walked (km, log scale)", y = "Species recorded (log scale)")
w2b_bottom <- ggplot(fit_resid, aes(km_band, exp(resid))) +
  geom_hline(yintercept = 1, colour = "steelblue", linewidth = 0.8) +
  geom_boxplot(fill = "grey85", colour = "grey30", outlier.size = 0.8, width = 0.6) +
  scale_y_continuous(labels = function(x) paste0(round((x - 1) * 100), "%")) +
  labs(x = "Distance walked (km), in five equal-sized groups of visits",
       y = "Species recorded relative\nto the fitted curve")
save_fig(w2b_top, "figW2b_species_distance_fit.png", width = 18, height = 12)
save_fig(w2b_bottom + labs(title = fig_title("Fit of the species-distance curve: departures by distance"),
                           subtitle = "Visits grouped by distance walked. 0% = exactly on the fitted curve."),
         "figW2c_species_distance_residuals.png", width = 18, height = 9)
print(w2b_top); print(w2b_bottom)

# --- W2d: The candidate curves compared ---
cmp_km <- seq(min(fit_data$km), max(c(fit_data$km, STANDARD_KM)) * 1.05, length.out = 150)
cmp <- map_dfr(names(curve_fns), ~ tibble(model = .x, km = cmp_km, S = curve_fns[[.x]](cmp_km))) %>%
  left_join(curve_stats, by = "model") %>%
  mutate(legend = sprintf("%s%s\n   typical error %.1f species; expects %.1f at %s",
                          label, if_else(model == RICHNESS_CURVE, "  (in use)", ""), rmse, expected, km_lab()))
w2d <- ggplot(fit_data, aes(km, S)) +
  annotate("rect", xmin = max(fit_data$km), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.12) +
  geom_point(colour = "grey60", alpha = 0.5, size = 1.4) +
  geom_line(data = cmp, aes(km, S, colour = legend), linewidth = 1) +
  geom_vline(xintercept = STANDARD_KM, linetype = "dashed", colour = "grey30") +
  geom_point(data = focal_visits, shape = 23, size = 3.5, fill = "steelblue", colour = "white") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(title = fig_title("Species-distance curves compared"),
       subtitle = str_wrap(sprintf("Grey points = benchmark visits; blue diamonds = %s visits; dashed line = the farm's distance (%s). Shaded = beyond the longest benchmark visit, where the curves are extrapolated. Typical error = root mean square difference between each curve and the visits. Where the curves disagree at the farm's distance, the expected species count is correspondingly uncertain.",
                                   SITE_NAME, km_lab()), 95),
       x = "Distance walked (km)", y = "Species recorded on the visit") +
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.key.height = unit(1, "cm"))
save_fig(w2d, "figW2d_species_distance_curves_compared.png", width = 18, height = 17)
print(w2d)

# --- W3: Richness and total birds, farm vs benchmark winter surveys ---
metric_levels <- c(paste0("Species per ", km_lab(), " visit"), paste0("Birds per ", km_lab(), " visit"))
w3_data <- bind_rows(
  bench_sw %>% transmute(metric = metric_levels[1], value = S_std),
  bench_sw %>% transmute(metric = metric_levels[2], value = total_std)) %>%
  filter(!is.na(value)) %>% mutate(metric = factor(metric, levels = metric_levels))
w3_farm <- tibble(metric = factor(metric_levels, levels = metric_levels),
                  value  = c(focal_sw$S_std, focal_sw$total_std),
                  pct    = c(richness_pct, total_pct)) %>%
  left_join(w3_data %>% count(metric, name = "n"), by = "metric") %>%
  mutate(lab = sprintf("Farm: %s\nhigher than %.0f%% of regional winter surveys\nn = %d winter surveys, %d farms",
                       ifelse(metric == metric_levels[1], sprintf("%.1f species", value), sprintf("%.0f birds", value)),
                       pct * 100, n, n_sites)) %>%
  filter(!is.na(value))

w3 <- ggplot(w3_data, aes(value)) +
  geom_histogram(bins = 20, fill = "grey75", colour = "white") +
  geom_vline(data = w3_farm, aes(xintercept = value), colour = "steelblue", linewidth = 1.1) +
  geom_label_nb(data = w3_farm, aes(x = if_else(pct >= 0.5, -Inf, Inf), y = Inf, label = lab,
                                    hjust = if_else(pct >= 0.5, 0, 1)), vjust = 1, size = 3.2,
                colour = "steelblue", fill = alpha("white", 0.75),
                label.padding = unit(0.4, "lines")) +
  # the strip under each panel acts as its x-axis label
  facet_wrap(~ metric, ncol = 1, scales = "free", strip.position = "bottom") +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(size = rel(1.05), colour = "grey20")) +
  labs(title = fig_title(paste0("Winter species richness and bird numbers at ", SITE_NAME, " vs the region")),
       subtitle = str_wrap(paste0("Each grey unit is one site in one winter. ", std_note, " Benchmark: ", bench_label, "."), 95),
       x = NULL, y = "Number of winter surveys")
save_fig(w3, "figW3_richness_and_abundance.png", width = 18, height = 16)

# --- W3b: Box plots with every winter survey shown ---
w3_box <- w3_data %>% group_by(metric) %>%
  summarise(ymin = quantile(value, 0.025), lower = quantile(value, 0.25), middle = median(value),
            upper = quantile(value, 0.75), ymax = quantile(value, 0.975), .groups = "drop")
# Drawn as two separate plots placed side by side, so each panel has its own
# y-axis title in the normal position (left of the numbers) in any ggplot2 version.
w3b_panel <- function(m, heading, ylab) {
  d <- w3_data %>% filter(metric == m); b <- w3_box %>% filter(metric == m); f <- w3_farm %>% filter(metric == m)
  ggplot(d, aes(x = "Region", y = value)) +
    geom_boxplot(data = b, inherit.aes = FALSE, stat = "identity", width = 0.5, alpha = 0.35,
                 fill = "tomato", colour = "grey30",
                 aes(x = "Region", ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax)) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 1), colour = "tomato3", size = 1.6, alpha = 0.7) +
    geom_hline(data = f, aes(yintercept = value), colour = "#2166AC", linetype = "dashed") +
    geom_point(data = f, aes(x = SITE_NAME, y = value), shape = 23, size = 6, fill = "#2166AC", colour = "white", stroke = 1.2) +
    scale_x_discrete(limits = c(SITE_NAME, "Region")) +
    labs(title = heading, x = NULL, y = ylab) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(face = "bold"))
}
w3b <- cbind(ggplotGrob(w3b_panel(metric_levels[1], "Species richness", "Number of species")),
             ggplotGrob(w3b_panel(metric_levels[2], "Bird abundance",  "Number of individual birds")), size = "max")
# record the caption (save_fig() only handles single ggplot objects)
figure_captions <- bind_rows(
  figure_captions[figure_captions$file != "figW3b_richness_and_abundance_boxplot.png", ],
  tibble(file  = "figW3b_richness_and_abundance_boxplot.png",
         title = paste0("Winter species richness and bird numbers: ", SITE_NAME, " vs the region"),
         notes = paste0("Each dot is one farm in one winter. Boxes show the middle 50% and the median; whiskers show the range containing 95% of winter surveys. ",
                        std_note, " Benchmark: ", bench_label, ".")))
ggsave(file.path(out_dir, "figW3b_richness_and_abundance_boxplot.png"), plot = w3b,
       width = 18, height = 12, units = "cm", dpi = 300, bg = "white")

# --- W4: Percentile heat strip by guild ---
pct_fill <- scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0.5,
                                 limits = c(0, 1), breaks = c(0.05, 0.5, 0.95), name = NULL,
                                 labels = c("Among the lowest\nin the region", "Typical", "Among the highest\nin the region"))
strip <- species_comparison %>%
  filter(!is.na(pct), recorded | p_std >= 0.25) %>%
  mutate(guild = factor(guild, levels = guild_levels),
         label = if_else(recorded, sprintf("%s  (region median %s)", comma(farm_mean, 0.1), comma(reg_median, 0.1)),
                         sprintf("not recorded  (region median %s)", comma(reg_median, 0.1))),
         english_name = fct_reorder(english_name, pct))
w4 <- ggplot(strip, aes(1, english_name, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label, colour = abs(pct - 0.5) > 0.3), size = 2.8) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  pct_fill + facet_grid(guild ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = fig_title("How each species on the farm compares with the region in winter"),
       subtitle = str_wrap(paste0("Colour = where the farm's mean count ranks among regional winter surveys. Text = farm mean per ", km_lab(),
                                  " visit (regional median). Includes regionally frequent species the farm did not record. Benchmark: ", bench_label, "."), 95),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title.position = "plot", axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
        legend.position = "bottom", legend.key.width = unit(2.2, "cm"))
save_fig(w4, "figW4_percentile_heat_strip.png", width = 18, height = max(12, nrow(strip) * 0.5 + 5))

# --- W4b: Percentile of each recorded species (equivalent of breeding fig 6b) ---
w4b_data <- species_comparison %>% filter(recorded, !is.na(pct)) %>%
  mutate(english_name = fct_reorder(english_name, pct))
w4b <- ggplot(w4b_data, aes(y = english_name)) +
  annotate("rect", xmin = 0.25, xmax = 0.75, ymin = -Inf, ymax = Inf, fill = "grey75", alpha = 0.35) +
  geom_vline(xintercept = 0.5, colour = "grey60", linetype = "dashed") +
  geom_segment(aes(x = 0.5, xend = pct, yend = english_name), colour = "grey70") +
  geom_point(aes(x = pct, colour = pct), size = 3) +
  scale_colour_gradient2(low = "tomato", mid = "grey50", high = "steelblue", midpoint = 0.5, guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = fig_title("How the site's winter counts rank against the region"),
       subtitle = str_wrap(paste0("Percentile of the site's mean count per ", km_lab(), " visit among regional winter surveys. ",
                                  "Shaded band = middle half of the region. Benchmark: ", bench_label, "."), 95),
       x = "Site percentile (0% = lowest in region, 100% = highest)", y = NULL)
save_fig(w4b, "figW4b_species_percentiles.png", width = 18, height = max(10, nrow(w4b_data) * 0.45))

# --- W5: Difference in mean count, site vs region (recorded species only) ---
# Bars: site mean count minus regional mean count (birds per standard visit), so
# each bar depends on that species alone. Labels: the same difference as a
# percentage of the regional mean (+100% = twice the regional mean).
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
w5 <- count_diff_plot(
  species_comparison %>% filter(recorded) %>% transmute(english_name, site = farm_mean, bench = reg_mean),
  fig_title("Differences in abundance compared to the regional winter bird community"),
  str_wrap(paste0("The 20 recorded species with the largest differences in number of birds. Bars: site mean count minus the regional mean count. ",
                  "Labels: the difference as a percentage of the regional mean (+100% = twice the regional mean; -50% = half). ",
                  std_note, " Benchmark: ", bench_label, "."), 95),
  "the region", "Difference in mean abundance per visit")
save_fig(w5, "figW5_divergence_regional.png", height = 16)

# --- W6: Abundance ratio ---
w6 <- species_comparison %>% filter(recorded, is.finite(ratio), ratio > 0) %>%
  mutate(log2_ratio = log2(ratio), english_name = fct_reorder(english_name, log2_ratio)) %>%
  ggplot(aes(log2_ratio, english_name)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", colour = "grey80") +
  geom_point(size = 2.5, colour = "steelblue") +
  labs(title = fig_title("Species abundance ratio: site vs region"),
       subtitle = str_wrap(paste0("0 = equal; +1 = twice as abundant at the site; −1 = half. ", std_note, " Benchmark: ", bench_label, "."), 95),
       x = "Log2(site mean / regional mean)", y = NULL)
save_fig(w6, "figW6_abundance_ratio.png", width = 18, height = max(10, farm_total_species * 0.42))

# --- W7: Regionally frequent species not recorded on the farm ---
det_axis <- "Regional encounter probability"
absent <- species_comparison %>% filter(!recorded, p_std > 0) %>% arrange(desc(p_std))
bocc_cols <- c("Red" = "#D95F02", "Amber" = "#E6AB02", "Green" = "#2E8B57", "Non-native" = "#7F6A93", "Not listed" = "grey55")
absent_native <- absent %>% filter(BOCC5 != "Non-native")   # figure and summary leave out non-native species
if (nrow(absent_native) > 0) {
  w7 <- absent_native %>% slice_head(n = 30) %>%
    mutate(english_name = fct_reorder(english_name, p_std), BOCC5 = factor(BOCC5, levels = names(bocc_cols))) %>%
    ggplot(aes(p_std, english_name, fill = BOCC5)) + geom_col(width = 0.6) +
    scale_fill_manual(values = bocc_cols, name = "BOCC5 status") +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = fig_title("Regional winter species not recorded at the site"),
         subtitle = str_wrap(paste0("The ", min(30, nrow(absent_native)), " most frequently recorded regional species absent from all ",
                                    focal_sw$n_visits, " farm visits in ", FOCAL_WINTER, " (non-native species not shown). Benchmark: ", bench_label, "."), 95),
         x = det_axis, y = NULL) + theme(legend.position = "bottom")
  save_fig(w7, "figW7_regional_absences.png", width = 18, height = 14)
}

# --- W8: How widespread the farm's species are regionally ---
w8 <- species_comparison %>% filter(recorded) %>% arrange(p_std) %>%
  mutate(english_name = fct_inorder(english_name)) %>%
  ggplot(aes(p_std, english_name)) + geom_col(fill = "steelblue", width = 0.6) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = fig_title("Regional encounter probability"),
       subtitle = str_wrap(paste0("Rarest regionally at the bottom of the scale. Benchmark: ", bench_label, "."), 95),
       x = det_axis, y = NULL)
save_fig(w8, "figW8_site_species_encounter_probability.png", width = 18, height = max(10, farm_total_species * 0.42))

# --- W9: The farm through time (if it has more than one winter) ---
farm_series <- site_winters %>% filter(is_focal) %>% arrange(winter_start(winter))
if (nrow(farm_series) >= 2) {
  bands <- w3_data %>% group_by(metric) %>%
    summarise(q25 = quantile(value, 0.25), q50 = median(value), q75 = quantile(value, 0.75), .groups = "drop")
  series_long <- bind_rows(
    farm_series %>% transmute(winter, metric = metric_levels[1], value = S_std),
    farm_series %>% transmute(winter, metric = metric_levels[2], value = total_std)) %>%
    filter(!is.na(value)) %>%
    mutate(metric = factor(metric, levels = metric_levels), winter = factor(winter, levels = farm_series$winter))
  w9 <- ggplot(series_long, aes(winter, value, group = 1)) +
    geom_rect(data = bands, inherit.aes = FALSE, aes(xmin = -Inf, xmax = Inf, ymin = q25, ymax = q75), fill = "grey80", alpha = 0.6) +
    geom_hline(data = bands, aes(yintercept = q50), linetype = "dashed", colour = "grey40") +
    geom_line(colour = "steelblue", linewidth = 0.9) + geom_point(colour = "steelblue", size = 3) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    labs(title = fig_title(paste0(SITE_NAME, " through time")),
         subtitle = str_wrap(paste0("Blue = the farm each winter. Grey band = middle 50% of regional winter surveys; dashed = regional median. ",
                                    std_note, " Benchmark: ", bench_label, "."), 95),
         x = "Winter", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(w9, "figW9_farm_through_time.png", width = 18, height = 15)
}

# --- W10: Species pool overlap, farm vs benchmark ---
bench_present <- detection %>% filter(bench_visits_present > 0) %>% pull(species_code)
shared_n     <- length(intersect(site_present, bench_present))
farm_only_n  <- length(setdiff(site_present, bench_present))
bench_only_n <- length(setdiff(bench_present, site_present))
overlap_results <- tibble(
  Comparison = "Site vs regional benchmark", Site_species = length(site_present),
  Benchmark_species = length(bench_present), Shared = shared_n, Site_only = farm_only_n,
  Benchmark_only = bench_only_n,
  Jaccard  = round(shared_n / (shared_n + farm_only_n + bench_only_n), 3),
  Sorensen = round(2 * shared_n / (length(site_present) + length(bench_present)), 3))
farm_only_species <- species_comparison %>% filter(recorded, !species_code %in% bench_present) %>%
  arrange(desc(farm_mean))

w10 <- tibble(category = factor(c("Recorded at both", "Site only", "Benchmark only"),
                                levels = c("Benchmark only", "Site only", "Recorded at both")),
              n_species = c(shared_n, farm_only_n, bench_only_n)) %>%
  ggplot(aes(x = "Site vs regional benchmark", y = n_species, fill = category)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = ifelse(n_species > 0, n_species, "")), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("Recorded at both" = "steelblue", "Site only" = "goldenrod", "Benchmark only" = "tomato"),
                    breaks = c("Recorded at both", "Site only", "Benchmark only")) +
  labs(title = fig_title("Species pool overlap: site vs regional benchmark"),
       subtitle = str_wrap(paste0("Jaccard: ", overlap_results$Jaccard, " | Sørensen: ", overlap_results$Sorensen,
                                  ". The benchmark pool is every species recorded across all benchmark visits, so it is much larger than any one site's list. Benchmark: ",
                                  bench_label, "."), 95),
       x = NULL, y = "Number of species", fill = NULL) +
  theme(legend.position = "bottom")
save_fig(w10, "figW10_species_overlap.png", width = 14, height = 12)

# --- Community similarity (Bray-Curtis on shares of all birds) ---
bc_similarity <- sum(pmin(species_comparison$farm_share, species_comparison$reg_share), na.rm = TRUE)
bench_bc <- bench_sp %>%
  filter(!is.na(mean_std)) %>%
  group_by(site, winter) %>% mutate(share = mean_std / sum(mean_std)) %>% ungroup() %>%
  left_join(species_comparison %>% select(species_code, reg_share), by = "species_code") %>%
  group_by(site, winter) %>%
  summarise(bc = sum(pmin(share, reg_share), na.rm = TRUE), .groups = "drop")
bc_pct <- mean(bench_bc$bc < bc_similarity, na.rm = TRUE)
message(sprintf("Bray-Curtis similarity to the regional winter community: %.3f (more typical than %.0f%% of benchmark winter surveys).",
                bc_similarity, bc_pct * 100))

# --- W11: The region's most widespread winter species, recorded or not ---
w11 <- species_comparison %>% filter(p_std > 0) %>% slice_max(p_std, n = 40, with_ties = FALSE) %>%
  mutate(english_name = fct_reorder(english_name, p_std),
         status = factor(if_else(recorded, "Recorded at site", "Absent from site"),
                         levels = c("Recorded at site", "Absent from site"))) %>%
  ggplot(aes(p_std, english_name, fill = status)) + geom_col() +
  scale_fill_manual(values = c("Recorded at site" = "steelblue", "Absent from site" = "tomato"), drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = fig_title("The region's most widespread winter species"),
       subtitle = str_wrap(paste0("Top 40 species by encounter probability on a ", km_lab(), " visit. Blue = recorded at the site in ",
                                  FOCAL_WINTER, "; red = not recorded. Benchmark: ", bench_label, "."), 95),
       x = det_axis, y = NULL, fill = NULL) +
  theme(legend.position = "bottom")
save_fig(w11, "figW11_regional_species_overview.png", width = 18, height = 22)

# =============================================================================
# 9. TABLES, CAPTIONS AND SUMMARY
# =============================================================================

species_table <- species_comparison %>%
  arrange(desc(recorded), desc(farm_mean)) %>%
  transmute(Species = english_name, Code = species_code, Guild = guild, BOCC5,
            `Recorded on farm` = if_else(recorded, "Yes", "No"),
            `Farm visits recorded` = visits_present,
            `Farm peak count` = farm_peak,
            `Farm mean per std visit` = round(farm_mean, 2),
            `Region mean` = round(reg_mean, 2), `Region median` = round(reg_median, 2),
            `Farm percentile (%)` = round(pct * 100),
            `Regional encounter probability (%)` = round(p_std * 100))
write_csv(species_table, file.path(out_dir, "winter_species_comparison.csv"))
write_csv(filter_log, file.path(out_dir, "winter_checklist_filter_log.csv"))
write_csv(site_winters %>% arrange(site, winter), file.path(out_dir, "winter_site_winters.csv"))

# --- Word tables ---
save_word_table <- function(ft, filename, title_text = NULL) {
  d <- read_docx()
  if (!is.null(title_text)) d <- body_add_par(d, title_text, style = "heading 2")
  print(body_add_flextable(d, ft), target = file.path(out_dir, filename))
}
ftab <- function(df, caption, digits2 = NULL, digits0 = NULL) {
  ft <- flextable(df)
  if (!is.null(digits2)) ft <- colformat_double(ft, j = digits2, digits = 2)
  if (!is.null(digits0)) ft <- colformat_double(ft, j = digits0, digits = 0)
  ft %>% theme_booktabs() %>% autofit() %>% set_caption(caption)
}

rec <- species_comparison %>% filter(recorded)
save_word_table(ftab(
  rec %>% arrange(english_name) %>%
    transmute(Species = english_name, Code = species_code, `Visits recorded` = visits_present,
              `Peak count` = farm_peak, `Mean count per visit` = farm_mean),
  sprintf("Winter birds at %s, %s: %d visits along the %g km transect.", SITE_NAME, FOCAL_WINTER, focal_sw$n_visits, SITE_TRANSECT_KM),
  digits2 = "Mean count per visit", digits0 = c("Visits recorded", "Peak count")),
  "tableW1_site_species.docx", "Table W1: Site species")

save_word_table(ftab(
  rec %>% arrange(desc(abs(log2(pmax(ratio, 1e-6))))) %>%
    transmute(Species = english_name, Code = species_code, `Site mean` = farm_mean, `Region mean` = reg_mean,
              `Region median` = reg_median, `Ratio site/region` = ratio, `Site percentile (%)` = round(pct * 100)),
  paste0("Species abundance comparison, mean birds per ", km_lab(), " visit. Benchmark: ", bench_label, "."),
  digits2 = c("Site mean", "Region mean", "Region median", "Ratio site/region"), digits0 = "Site percentile (%)"),
  "tableW2_abundance_comparison.docx", "Table W2: Abundance comparison")

save_word_table(ftab(
  overlap_results %>% rename(`Site species` = Site_species, `Benchmark species` = Benchmark_species,
                             `Site only` = Site_only, `Benchmark only` = Benchmark_only, `Sørensen` = Sorensen),
  paste0("Species pool overlap. Jaccard = shared / all species in either list; Sørensen = 2 x shared / (site + benchmark totals). ",
         "Bray-Curtis similarity of the site to the regional community: ", round(bc_similarity, 3),
         " (more typical than ", round(bc_pct * 100), "% of benchmark winter surveys). Benchmark: ", bench_label, ".")),
  "tableW3_species_overlap.docx", "Table W3: Species pool overlap and community similarity")

rng <- function(x) sprintf("%.1f–%.1f", quantile(x, 0.025, na.rm = TRUE), quantile(x, 0.975, na.rm = TRUE))
save_word_table(ftab(
  tibble(Measure = c(paste0("Species per ", km_lab(), " visit"), paste0("Birds per ", km_lab(), " visit")),
         Site = c(focal_sw$S_std, focal_sw$total_std),
         `Region median` = c(median(bench_sw$S_std), median(bench_sw$total_std, na.rm = TRUE)),
         `Region 95% range` = c(rng(bench_sw$S_std), rng(bench_sw$total_std)),
         `Site percentile (%)` = round(c(richness_pct, total_pct) * 100)),
  sprintf("Species richness and bird numbers per visit. Species expected on a regional visit of %s from the species-distance curve: %.1f (%s). The site recorded %d species over the whole winter. Benchmark: %s.",
          km_lab(), expected_S, curve_desc, farm_total_species, bench_label),
  digits2 = NULL, digits0 = "Site percentile (%)") %>% colformat_double(j = c("Site", "Region median"), digits = 1),
  "tableW4_richness_and_abundance.docx", "Table W4: Species richness and bird numbers")

if (nrow(absent) > 0)
  save_word_table(ftab(
    absent %>% slice_head(n = 30) %>%
      transmute(Species = english_name, Code = species_code, BOCC5,
                `Regional encounter probability (%)` = round(p_std * 100), `Region mean count` = reg_mean),
    paste0("The 30 most frequently recorded regional winter species not recorded at the site. Encounter probability = chance a ",
           km_lab(), " visit in the region records the species. Benchmark: ", bench_label, "."),
    digits2 = "Region mean count", digits0 = "Regional encounter probability (%)"),
    "tableW5_absent_regional_species.docx", "Table W5: Regional species absent from site")

if (nrow(farm_only_species) > 0) {
  save_word_table(ftab(
    farm_only_species %>% transmute(Species = english_name, Code = species_code, BOCC5,
                                    `Visits recorded` = visits_present, `Peak count` = farm_peak),
    paste0("Species recorded at the site but never in the regional benchmark. Benchmark: ", bench_label, "."),
    digits0 = c("Visits recorded", "Peak count")),
    "tableW6_site_only_species.docx", "Table W6: Site-only species")
} else message("No site-only species - table W6 skipped.")

save_word_table(ftab(
  species_comparison %>% filter(p_std > 0) %>% arrange(desc(p_std)) %>%
    transmute(Species = english_name, Code = species_code, `Recorded at site` = if_else(recorded, "Yes", "No"),
              `Regional encounter probability (%)` = round(p_std * 100), `Region mean count` = reg_mean),
  paste0("Regional winter species: encounter probability on a ", km_lab(), " visit. Benchmark: ", bench_label, "."),
  digits2 = "Region mean count", digits0 = "Regional encounter probability (%)"),
  "tableW7_regional_encounter_probability.docx", "Table W7: Regional species encounter probability")

message("Tables saved.")


write_csv(figure_captions, file.path(out_dir, "figure_captions.csv"))
cap_doc <- read_docx() %>%
  body_add_par(paste0("Figure captions — ", SITE_NAME, ", winter ", FOCAL_WINTER), style = "heading 1")
for (i in seq_len(nrow(figure_captions))) {
  cap_doc <- cap_doc %>%
    body_add_par(figure_captions$file[i], style = "heading 3") %>%
    body_add_fpar(fpar(ftext(figure_captions$title[i], fp_text(bold = TRUE, font.size = 11)))) %>%
    body_add_par(figure_captions$notes[i], style = "Normal")
}
print(cap_doc, target = file.path(out_dir, "figure_captions.docx"))

top_over  <- species_comparison %>% filter(recorded, diff_share > 0) %>% slice_max(diff_share, n = 5) %>% pull(english_name)
top_under <- species_comparison %>% filter(recorded, diff_share < 0) %>% slice_min(diff_share, n = 5) %>% pull(english_name)
top_absent <- head(absent_native$english_name, 5)
body <- function(...) fpar(ftext(paste0(...), fp_text(font.size = 11, font.family = "Arial")))
head1 <- function(x) fpar(ftext(x, fp_text(bold = TRUE, font.size = 13, font.family = "Arial", color = "#2E75B6")))

doc <- read_docx() %>%
  body_add_fpar(fpar(ftext(paste0("Winter Bird Community Analysis — ", SITE_NAME),
                           fp_text(bold = TRUE, font.size = 16, font.family = "Arial", color = "#2E75B6")))) %>%
  body_add_fpar(body("Winter ", FOCAL_WINTER, "  |  Benchmark: ", bench_label, ".")) %>%
  body_add_par("") %>%
  body_add_fpar(head1("1. Species richness and bird numbers")) %>%
  body_add_fpar(body(SITE_NAME, " recorded ", farm_total_species, " species over ", focal_sw$n_visits,
                     " visits in winter ", FOCAL_WINTER, ", averaging ", round(focal_sw$S_std, 1),
                     " species per visit along the ", SITE_TRANSECT_KM, " km transect. ",
                     "A regional survey of the same distance would be expected to record about ", round(expected_S, 1),
                     " species (from the species-distance curve). Against regional winter surveys adjusted to ", km_lab(),
                     " (median ", round(median(bench_sw$S_std), 1), " species), the farm was higher than ",
                     round(richness_pct * 100), "%",
                     if (!is.na(focal_sw$total_std)) paste0(". It averaged ", round(focal_sw$total_std), " birds per visit, higher than ",
                                                            round(total_pct * 100), "% of regional winter surveys (median ",
                                                            round(median(bench_sw$total_std, na.rm = TRUE)), ")") else "", ".")) %>%
  body_add_par("") %>%
  body_add_fpar(head1("2. Community composition")) %>%
  body_add_fpar(body("Relative to the region, the species making up a larger share of birds on the farm were: ",
                     paste(top_over, collapse = ", "), ". Those making up a smaller share were: ",
                     paste(top_under, collapse = ", "), ".")) %>%
  body_add_par("") %>%
  body_add_fpar(body("Community similarity to the regional winter bird community (Bray-Curtis, on each species' share of all birds) was ",
                     round(bc_similarity, 3), " (1 = identical make-up), more typical than ", round(bc_pct * 100),
                     "% of benchmark winter surveys. ", shared_n, " of the site's ", length(site_present),
                     " species were also in the regional benchmark pool of ", length(bench_present),
                     " species (Jaccard ", overlap_results$Jaccard, ", S\u00f8rensen ", overlap_results$Sorensen, ").",
                     if (nrow(farm_only_species) > 0) paste0(" Recorded at the site but never in the benchmark: ",
                                                             paste(farm_only_species$english_name, collapse = ", "), ".") else "")) %>%
  body_add_par("") %>%
  body_add_fpar(head1("3. Species not recorded")) %>%
  body_add_fpar(body(if (length(top_absent) > 0)
    paste0("The regionally most frequent winter species not recorded on the farm were: ", paste(top_absent, collapse = ", "), ".")
    else "All regionally recorded species were found on the farm.")) %>%
  body_add_par("") %>%
  body_add_fpar(head1("4. Methods")) %>%
  body_add_fpar(body("The farm was surveyed along a ", SITE_TRANSECT_KM, " km transect. Benchmark data are winter farmland ",
                     "bird surveys by Landmark Ecology held in eBird (",
                     paste(month.abb[WINTER_MONTHS], collapse = ", "), "; complete travelling checklists",
                     if (EXCLUDE_NOCTURNAL) paste0(" in daylight (within ", TWILIGHT_BUFFER_MIN, " minutes of sunrise and sunset)") else "", ", of ",
                     MIN_KM, "–", MAX_KM, " km",
                     if (speed_filter_on) paste0("; ", tolower(speed_label)) else "", "). Because the distance walked varies between benchmark surveys, each was adjusted to ",
                     km_lab(), ": counts were scaled in proportion to distance, and species richness was adjusted using a ",
                     "species-distance curve fitted to the benchmark (", curve_desc, "). ",
                     "The encounter probability of each species (the chance of recording it on a ", km_lab(), " visit) was estimated from presence and absence ",
                     "on all benchmark visits using an encounter-rate model. Counts entered as 'X' were treated as present ",
                     "and excluded from abundance calculations. Each site in each winter with at least ", MIN_VISITS,
                     " visits forms one benchmark unit (the average across its visits); the farm is ranked against these. ",
                     "The benchmark is drawn from surveyed farmland sites, not a random sample of the countryside."))
print(doc, target = file.path(out_dir, "winter_summary_report.docx"))

message("\nWinter analysis complete. Outputs saved to: ", out_dir)

# =============================================================================
# DATA CHECKS (RUN_CHECKS): tracking down odd points on figure W2
# Changes nothing in the analysis above. After it has run, in the console:
#   find_point(km = c(11, 13), species = c(5, 15))   look a point up by its position
#   inspect_visit(123)                               everything about one visit
# =============================================================================
if (RUN_CHECKS) {
  # --- 1. Every benchmark visit, with the details needed to judge it -------------
  visit_detail <- checklists %>%
    group_by(site, date) %>%
    summarise(checklist_ids = paste(checklist, collapse = ", "),
              start         = paste(unique(na.omit(time_str)), collapse = ", "),
              minutes       = sum(duration, na.rm = TRUE),
              county        = first(county),
              .groups = "drop")
  
  w2_points <- visits %>%
    filter(!is_focal) %>%
    left_join(visit_detail, by = c("site", "date")) %>%
    mutate(expected     = curve_f(km),                # species the curve expects at this distance
           vs_curve_pct = round((S / expected - 1) * 100),
           km_per_hour  = if_else(minutes > 0, round(km / (minutes / 60), 1), NA_real_),
           species_per_km = round(S / km, 1)) %>%
    select(visit_id, site, county, date, start, km, minutes, km_per_hour, n_checklists,
           species = S, expected, vs_curve_pct, species_per_km, has_X_counts = any_x, checklist_ids) %>%
    mutate(expected = round(expected, 1))
  
  # --- 2. The visits furthest from the curve --------------------------------------
  cat("\nVisits furthest from the curve (most extreme first):\n")
  w2_points %>% arrange(desc(abs(log(species / expected)))) %>% select(-checklist_ids) %>%
    head(15) %>% as.data.frame() %>% print(row.names = FALSE)
  
  # --- 3. Things that often explain an odd point ----------------------------------
  cat("\nFastest and slowest visits still in the benchmark (a typo in distance or duration is likely at the extremes):\n")
  bind_rows(slice_max(w2_points, km_per_hour, n = 8), slice_min(w2_points, km_per_hour, n = 8)) %>%
    distinct() %>% select(visit_id, site, date, km, minutes, km_per_hour, species) %>%
    as.data.frame() %>% print(row.names = FALSE)
  
  cat("\nVisits made of more than one checklist (distances are added together;",
      "check these are not duplicates or shared copies of the same walk):\n")
  w2_points %>% filter(n_checklists > 1) %>% arrange(desc(n_checklists)) %>%
    select(visit_id, site, date, n_checklists, km, species, checklist_ids) %>%
    as.data.frame() %>% print(row.names = FALSE)
  
  cat("\nLongest and shortest visits:\n")
  bind_rows(slice_max(w2_points, km, n = 5), slice_min(w2_points, km, n = 5)) %>%
    select(visit_id, site, date, km, minutes, species, expected) %>% as.data.frame() %>% print(row.names = FALSE)
  
  # --- 4. Find a point by where it sits on the graph ------------------------------
  # Give the rough range you can read off figure W2, e.g. find_point(km = c(11, 13), species = c(5, 15))
  find_point <- function(km = c(-Inf, Inf), species = c(-Inf, Inf)) {
    out <- w2_points %>% filter(between(.data$km, km[1], km[2]), between(.data$species, species[1], species[2]))
    print(as.data.frame(out), row.names = FALSE)
    invisible(out)
  }
  
  # --- 5. Everything about one visit ----------------------------------------------
  # Use the visit_id from the tables above, e.g. inspect_visit(123)
  inspect_visit <- function(id) {
    v <- w2_points %>% filter(visit_id == id)
    if (nrow(v) == 0) { cat("No benchmark visit with that id.\n"); return(invisible(NULL)) }
    cat(sprintf("\nVisit %d: %s, %s (%s)\n", id, v$site, format(v$date), v$county))
    cat(sprintf("  %.2f km in %s minutes (%s km/h), start %s, %d checklist(s)\n",
                v$km, v$minutes, v$km_per_hour, v$start, v$n_checklists))
    cat(sprintf("  %d species recorded; the curve expects %.1f at this distance (%+d%%)\n",
                v$species, v$expected, v$vs_curve_pct))
    ids <- str_split(v$checklist_ids, ", ")[[1]]
    cat("  Open in eBird:\n", paste0("    https://ebird.org/checklist/", ids, collapse = "\n"), "\n")
    ch <- checklists %>% filter(checklist %in% ids) %>%
      select(checklist, location, date, time_str, duration, km, protocol, complete)
    cat("\n  Checklists in this visit:\n"); print(as.data.frame(ch), row.names = FALSE)
    spp <- sp_obs %>% filter(checklist %in% ids) %>%
      group_by(Species = english_name) %>%
      summarise(Count = paste(count_raw, collapse = " + "), .groups = "drop") %>% arrange(Species)
    cat(sprintf("\n  Species list (%d):\n", nrow(spp))); print(as.data.frame(spp), row.names = FALSE)
    others <- w2_points %>% filter(site == v$site, visit_id != id) %>% arrange(date) %>%
      select(visit_id, date, km, minutes, species, expected)
    cat(sprintf("\n  Other visits to %s for comparison (%d):\n", v$site, nrow(others)))
    print(as.data.frame(others), row.names = FALSE)
    invisible(v)
  }
  
  # --- 6. Figure W2 again, with the most extreme visits labelled -------------------
  extreme <- w2_points %>% arrange(desc(abs(log(species / expected)))) %>% head(CHECK_N_LABEL) %>%
    mutate(lab = sprintf("%d: %s\n%s", visit_id, str_trunc(site, 22), format(date)))
  p_check <- ggplot(w2_points, aes(km, species)) +
    geom_point(colour = "grey60", alpha = 0.6) +
    geom_line(data = tibble(km = seq(min(w2_points$km), max(w2_points$km), length.out = 100)) %>%
                mutate(species = curve_f(km)), colour = "grey20") +
    geom_point(data = extreme, colour = "tomato", size = 2.5) +
    (if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(data = extreme, aes(label = lab), size = 2.8, colour = "tomato4",
                               lineheight = 0.9, min.segment.length = 0, max.overlaps = 30)
     else geom_text(data = extreme, aes(label = lab), size = 2.8, colour = "tomato4", vjust = -0.6, lineheight = 0.9)) +
    labs(title = paste("Check C1: figure W2 with the", CHECK_N_LABEL, "visits furthest from the curve labelled"),
         subtitle = "Label = visit_id: site, date. Use inspect_visit(visit_id) for the detail.",
         x = "Distance walked (km)", y = "Species recorded on the visit")
  print(p_check)
  ggsave(file.path(out_dir, "checkC1_w2_extreme_visits.png"), p_check, width = 22, height = 15, units = "cm", dpi = 200, bg = "white")
  
  # --- 7. Walking speeds of the checklists (before the speed filter), to help set the limits ---
  sp_chk <- speed_before %>% filter(!is.na(km_per_hour))
  if (nrow(sp_chk) > 0) {
    p_speed <- ggplot(sp_chk, aes(km_per_hour)) +
      geom_histogram(bins = 40, fill = "grey60", colour = "white") +
      geom_vline(xintercept = c(MIN_KM_PER_HOUR, MAX_KM_PER_HOUR)[!is.na(c(MIN_KM_PER_HOUR, MAX_KM_PER_HOUR))],
                 linetype = "dashed", colour = "tomato") +
      labs(title = "Check C2: walking speed of benchmark checklists",
           subtitle = "Dashed lines = MIN_KM_PER_HOUR / MAX_KM_PER_HOUR (if set); checklists outside them are removed",
           x = "km per hour", y = "Checklists")
    print(p_speed)
    ggsave(file.path(out_dir, "checkC2_walking_speed.png"), p_speed, width = 18, height = 11, units = "cm", dpi = 200, bg = "white")
    cat("\nWalking speed of checklists before the speed filter (km/h):\n")
    print(round(quantile(sp_chk$km_per_hour, c(0, .01, .05, .25, .5, .75, .95, .99, 1)), 2))
    cat(sum(is.na(speed_before$km_per_hour)), "checklists have no duration, so no speed (kept).\n")
  }
  
  write_csv(w2_points %>% arrange(desc(abs(log(species / expected)))),
            file.path(out_dir, "winter_w2_points_check.csv"))
  cat("\nAll benchmark visits saved to winter_w2_points_check.csv (most extreme first).\n",
      "To drop one checklist, add its ID to EXCLUDE_CHECKLISTS in the config.\n",
      "Next: find_point(km = c(lo, hi), species = c(lo, hi)) or inspect_visit(visit_id).\n")
}