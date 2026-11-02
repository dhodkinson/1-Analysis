rm(list=ls())

# =============================================================================
# Birds by crop, habitat and woodland type - winter and summer
#
# Uses the farm's own counts only (GeoPackage and/or CSV, BTO codes): no
# benchmark data. Companion to 7_Winter_Community-eBird_v1.R.
#
# For each season in RUN_SEASONS (W = winter, S = summer in the file names):
#   fig?13   the commonest species and the crops they were recorded in
#   fig?13b  the same, as birds per 10 hectares (needs crop areas)
#   fig?14   the same for the woodland habitats you list
#   fig?14b  woodland per 10 hectares (needs areas for the woodland habitats)
#   table?8  birds by crop
# =============================================================================

library(tidyverse)
library(scales)
library(officer)
library(flextable)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# --- Files ---
SPECIES_LOOKUP <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/1 Analysis/Data cleaning/BTO_species_codes_lookup.csv"
out_dir        <- "/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results/Crops"

# --- The farm ---
SITE_NAME  <- "Shimpling Park Farm"
# Farm survey data: GeoPackage and/or CSV (rows are combined). NA to skip either.
input_path <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'
if (!is.na(input_path)) { message("Layers in the GeoPackage:"); print(sf::st_layers(input_path)$name) }   # pick one of these
layer_name <- '2026_Counts_final'
csv_path   <- NA
CSV_COLS   <- c(species_code = "BTO_code",   # column names in the CSV
                transect     = "Site",
                count        = "Count",
                date_str     = "Date")

# --- Columns in the farm data ---
CROP_COLUMN     <- "Habitat"   # crop / habitat each record was in
# Flying-over records: the column that marks them and the value(s) that mean
# "flying over". Each season below says whether they are left out.
FLYOVER_COLUMN  <- "Status"
FLYOVER_VALUES  <- c("F")

# --- Crop areas (for birds per 10 ha). One row per field or per crop; areas are
# summed for each crop. NA = no density figures. Each season names its own crop
# column in this file (area_crop_col below). ---
CROP_AREA_FILE    <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Data/Shimpling_field_areas.csv'
CROP_AREA_HA_COL  <- "Area_ha"
# ... or take the areas from a polygon layer in the farm GeoPackage (input_path).
# When CROP_AREA_LAYER is set it is used instead of CROP_AREA_FILE. Each season's
# area_crop_col must then name a column in that layer. CROP_AREA_LAYER_HA_COL = NA
# works the area out from the polygon shapes; or name a column holding hectares.
CROP_AREA_LAYER        <- NA      # e.g. "2026_Crop_polygons"
CROP_AREA_LAYER_HA_COL <- NA

# --- Seasons ---
RUN_SEASONS <- c("winter", "summer")
# For each season:
#   months            months of the year (1 = January)
#   focal             which winter / year, e.g. "2025/26" or "2026"; NA = most recent
#   stat              "mean" = mean birds per visit; "max" = highest visit count
#                     (birds summed per species and crop on each visit, then the highest visit)
#   exclude_flyovers  leave out records marked as flying over
#   area_crop_col     the crop-name column for this season in CROP_AREA_FILE
#   relabel           c("name in the data" = "name to use") for this season only; NA = none.
#                     Needed where CROP_COLUMN holds one season's crop name (e.g. the
#                     spring oats fields, which were ploughed over winter).
#   crops             crops to show: set in CROP_INCLUDE_WINTER / CROP_INCLUDE_SUMMER below
#   exclude_species   species left out of the crop figures and table (English names
#                     or BTO codes); NA = none. Woodland figures are not affected.
#   woodland          habitats for the separate woodland figures; NA = none
# Crops to show in each season (names as in CROP_COLUMN, after any relabelling). NA = all.
CROP_INCLUDE_WINTER <- c("SAM3 Herbal Ley", "LG20 Field Beans (winter)", "AC66 Wheat (winter)",
                         "PG01: Permanent Grassland", "Vetch", "Ploughed")
CROP_INCLUDE_SUMMER <- c("SAM3 Herbal Ley", "LG20 Field Beans (winter)", "AC66 Wheat (winter)",
                         "PG01: Permanent Grassland", "Vetch", "AC19 Oats (spring)")

SEASONS <- list(
  winter = list(
    months = c(11, 12, 1, 2), focal = NA, stat = "mean", exclude_flyovers = FALSE,
    area_crop_col = "Winter_crop",
    # The crop layer gives these fields their spring crop all year; in winter they were bare plough
    relabel = c("AC19 Oats (spring)" = "Ploughed"),
    crops = CROP_INCLUDE_WINTER,
    exclude_species = NA,     # e.g. c("Blue Tit", "Great Tit", "Robin")
    woodland = c("Woodland: More mature", "Woodland: Less mature", "Woodland: Recent felling")),
  summer = list(
    months = c(4, 5, 6), focal = NA, stat = "max", exclude_flyovers = TRUE,
    area_crop_col = "Spring_crop",
    relabel = NA,   # the crop layer already labels these fields "AC19 Oats (spring)"
    crops = CROP_INCLUDE_SUMMER,
    exclude_species = NA,
    woodland = c("Woodland: More mature", "Woodland: Less mature", "Woodland: Recent felling")))
WOODLAND_LABEL <- "woodland type"   # how the woodland group is described in captions
N_SPECIES      <- 30                # species shown in each figure

# --- Figures ---
# Titles and subtitles are left off the images and written to
# crop_figure_captions.docx / .csv instead. Set TRUE to print them on the figures.
SHOW_FIGURE_TITLES <- FALSE
fig_title <- function(txt) txt

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

key <- function(x) tolower(str_squish(x))
# Winter spans two years ("2025/26", from July); other seasons are labelled by year
period_label <- function(d, is_winter) {
  y <- if_else(month(d) >= 7, year(d), year(d) - 1)
  if (is_winter) sprintf("%d/%02d", y, (y + 1) %% 100) else as.character(year(d))
}

# =============================================================================
# 1. LOAD THE FARM SURVEY DATA
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

farm_all <- farm_sources %>%
  map(~ transmute(.x, species_code = str_trim(as.character(species_code)),
                  date = parse_survey_date(date_str),
                  count = replace_na(suppressWarnings(as.numeric(count)), 0),
                  crop = if (CROP_COLUMN %in% names(.x)) as.character(.x[[CROP_COLUMN]]) else NA_character_,
                  flyover = if (!is.na(FLYOVER_COLUMN) && FLYOVER_COLUMN %in% names(.x))
                    str_trim(as.character(.x[[FLYOVER_COLUMN]])) %in% FLYOVER_VALUES else FALSE)) %>%
  bind_rows()
if (all(is.na(farm_all$crop))) stop("Column '", CROP_COLUMN, "' was not found (or is empty) in the farm data. Check CROP_COLUMN.")
if (any(is.na(farm_all$date))) warning(sum(is.na(farm_all$date)), " farm rows have unreadable dates and are excluded.")
farm_all <- farm_all %>% filter(!is.na(date), count > 0) %>%
  mutate(crop = coalesce(na_if(str_trim(crop), ""), "Not recorded"))

lookup_by_code <- read_csv(SPECIES_LOOKUP, col_types = cols(.default = "c")) %>%
  transmute(species_code = str_trim(bto_code), english_name = str_trim(species),
            BOCC5 = str_to_title(str_trim(bocc))) %>%
  mutate(BOCC5 = recode(BOCC5, "Non-Native" = "Non-native"),
         BOCC5 = if_else(BOCC5 %in% c("Red", "Amber", "Green", "Non-native"), BOCC5, "Not listed")) %>%
  distinct(species_code, .keep_all = TRUE)
farm_unknown <- setdiff(unique(farm_all$species_code), lookup_by_code$species_code)
if (length(farm_unknown) > 0)
  message("Farm species codes not in the BTO lookup (shown by code): ", paste(farm_unknown, collapse = ", "))
farm_all <- farm_all %>% left_join(lookup_by_code, by = "species_code") %>%
  mutate(english_name = coalesce(english_name, species_code))

area_raw <- NULL
if (!is.na(CROP_AREA_LAYER)) {
  gl <- sf::st_layers(input_path)$name
  if (!CROP_AREA_LAYER %in% gl)
    stop("CROP_AREA_LAYER '", CROP_AREA_LAYER, "' is not in the GeoPackage. Layers available: ", paste(gl, collapse = ", "))
  area_sf  <- sf::st_read(input_path, layer = CROP_AREA_LAYER, quiet = TRUE)
  area_raw <- sf::st_drop_geometry(area_sf) %>% as_tibble() %>% mutate(across(everything(), as.character))
  if (is.na(CROP_AREA_LAYER_HA_COL)) {
    area_raw$.area_from_shape <- as.character(as.numeric(sf::st_area(area_sf)) / 10000)
    CROP_AREA_HA_COL <- ".area_from_shape"
  } else CROP_AREA_HA_COL <- CROP_AREA_LAYER_HA_COL
  message(sprintf("Crop areas from layer '%s': %d polygons, %s ha in total. Columns: %s", CROP_AREA_LAYER, nrow(area_raw),
                  comma(sum(suppressWarnings(as.numeric(area_raw[[CROP_AREA_HA_COL]])), na.rm = TRUE), 0.1),
                  paste(setdiff(names(area_raw), ".area_from_shape"), collapse = ", ")))
} else if (!is.na(CROP_AREA_FILE)) {
  area_raw <- read_csv(CROP_AREA_FILE, col_types = cols(.default = "c"), show_col_types = FALSE)
  if (!CROP_AREA_HA_COL %in% names(area_raw))
    stop("Area file has no column '", CROP_AREA_HA_COL, "'. Columns found: ", paste(names(area_raw), collapse = ", "))
}

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



# =============================================================================
# 2. ONE SET OF FIGURES (crops or woodland) FOR ONE SEASON
# =============================================================================
# value for each species in each habitat: birds are summed on each visit, then
# averaged over all the season's visits ("mean") or the highest visit taken ("max")
habitat_values <- function(d, stat, n_visits) {
  d %>% group_by(crop, species_code, english_name, date) %>% summarise(t = sum(count), .groups = "drop") %>%
    group_by(crop, species_code, english_name) %>%
    summarise(v = if (stat == "max") max(t) else sum(t) / n_visits, .groups = "drop")
}
habitat_grid <- function(d, species, levels, fill_name, colour, lab_fun = identity) {
  d %>% filter(english_name %in% species) %>% select(crop, english_name, v) %>%
    complete(crop, english_name = species, fill = list(v = 0)) %>%
    mutate(crop = factor(crop, levels = intersect(levels, unique(crop))),
           english_name = factor(english_name, levels = rev(species))) %>%
    ggplot(aes(crop, english_name, fill = v)) + geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(v > 0, ifelse(v < 0.05, "<0.1", comma(v, 0.1)), "")), size = 2.8) +
    scale_x_discrete(labels = lab_fun) +
    scale_fill_gradient(low = "white", high = colour, trans = "sqrt", name = fill_name) +
    theme_minimal(base_size = 10) +
    theme(plot.title.position = "plot", panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

run_set <- function(obs, include, areas, n_visits, cfg, season, when, set = c("crop", "woodland")) {
  set <- match.arg(set)
  L <- toupper(substr(season, 1, 1))                       # W or S, for file names
  num <- if (set == "crop") "13" else "14"
  what <- if (set == "crop") tolower(CROP_COLUMN) else WOODLAND_LABEL
  colour <- if (set == "crop") "steelblue" else "#4C9F70"
  unit <- if (cfg$stat == "max") "Highest visit count" else "Mean birds per visit"
  unit_txt <- if (cfg$stat == "max") "the highest count on any one visit" else "mean birds per visit"
  
  if (!all(is.na(include))) {
    miss <- include[!key(include) %in% key(unique(obs$crop))]
    if (length(miss) > 0) warning(season, " ", set, " names not found in the data: ", paste(miss, collapse = ", "), call. = FALSE)
    obs <- obs %>% filter(key(crop) %in% key(include))
  }
  excl <- if (set == "crop") cfg$exclude_species else NA
  if (!all(is.na(excl))) {
    out <- key(obs$english_name) %in% key(excl) | obs$species_code %in% excl
    miss <- excl[!(key(excl) %in% key(farm_all$english_name) | excl %in% farm_all$species_code)]
    if (length(miss) > 0) warning("exclude_species not recognised (", season, "): ", paste(miss, collapse = ", "), call. = FALSE)
    if (any(out)) message("  Species left out of the ", season, " crop figures: ", paste(sort(unique(obs$english_name[out])), collapse = ", "),
                          sprintf(" (%d records, %s birds).", sum(out), comma(sum(obs$count[out]))))
    obs <- obs[!out, ]
  }
  if (nrow(obs) == 0) { warning("No ", season, " records for the ", set, " figures - skipped.", call. = FALSE); return(invisible(NULL)) }
  
  vals <- habitat_values(obs, cfg$stat, n_visits)
  summary <- vals %>% group_by(crop) %>%
    summarise(birds = sum(v), species = n_distinct(species_code),
              top_species = paste(head(english_name[order(-v)], 3), collapse = ", "), .groups = "drop") %>%
    mutate(share = birds / sum(birds))
  # column order: as listed in the config, or by number of birds
  levels <- if (!all(is.na(include))) summary$crop[order(match(key(summary$crop), key(include)))]
  else summary$crop[order(-summary$birds)]
  species <- vals %>% count(english_name, wt = v, sort = TRUE) %>% slice_head(n = N_SPECIES) %>% pull(english_name)
  notes <- paste0(SITE_NAME, ", ", season, " ", when, ", ", n_visits, " visits.",
                  if (cfg$exclude_flyovers) " Birds flying over are excluded." else "",
                  if (!all(is.na(excl))) paste0(" Species excluded: ", paste(excl, collapse = ", "), ".") else "")
  
  p <- habitat_grid(vals, species, levels, gsub(" ", "\n", unit, fixed = TRUE), colour) +
    labs(title = fig_title(paste0("The commonest ", season, " species by ", what)),
         subtitle = str_wrap(paste0("For the ", length(species), " most numerous species, ", unit_txt, " in each ", what,
                                    ". Not adjusted for the extent of each ", what, ". ", notes), 95),
         x = NULL, y = NULL)
  save_fig(p, sprintf("fig%s%s_species_by_%s.png", L, num, set), width = 18, height = 16)
  
  if (!is.null(areas)) {
    summary <- summary %>% mutate(k = key(crop)) %>% left_join(areas, by = "k") %>% select(-k) %>%
      mutate(birds_per_10ha = birds / area_ha * 10)
    dens <- summary %>% filter(!is.na(area_ha))
    no_area <- summary$crop[is.na(summary$area_ha)]
    if (length(no_area) > 0)
      message("  No area in ", cfg$area_crop_col, " for: ", paste(no_area, collapse = ", "), " (left out of the density figure).")
    if (nrow(dens) > 0) {
      pb <- vals %>% inner_join(dens %>% select(crop, area_ha), by = "crop") %>% mutate(v = v / area_ha * 10) %>%
        habitat_grid(species, levels, paste0(gsub(" ", "\n", unit, fixed = TRUE), "\nper 10 ha"), colour,
                     function(x) paste0(x, "\n(", comma(dens$area_ha[match(x, dens$crop)], 0.1), " ha)")) +
        labs(title = fig_title(paste0("The commonest ", season, " species by ", what, ": density")),
             subtitle = str_wrap(paste0("For the ", length(species), " most numerous species, ", unit_txt, " per 10 hectares of each ", what,
                                        ". Small areas give less reliable densities. ", notes), 95),
             x = NULL, y = NULL)
      save_fig(pb, sprintf("fig%s%sb_species_density_by_%s.png", L, num, set), width = 18, height = 16)
    }
  }
  summary <- summary %>% arrange(match(crop, levels))
  write_csv(summary, file.path(out_dir, sprintf("%s_birds_by_%s.csv", season, set)))
  
  if (set == "crop") {
    has_area <- "area_ha" %in% names(summary) && any(!is.na(summary$area_ha))
    tb <- summary %>% transmute(!!CROP_COLUMN := crop, !!unit := birds)
    if (has_area) tb <- tb %>% mutate(`Area (ha)` = summary$area_ha, `Birds per 10 ha` = summary$birds_per_10ha)
    tb <- tb %>% mutate(`Share of birds (%)` = round(summary$share * 100), Species = summary$species,
                        `Most numerous species` = summary$top_species)
    save_word_table(
      ftab(tb, paste0("Birds by ", what, ". ", unit, " = ", unit_txt, ", added up across species. ", notes),
           digits0 = c("Share of birds (%)", "Species")) %>%
        colformat_double(j = intersect(c(unit, "Area (ha)", "Birds per 10 ha"), names(tb)), digits = 1),
      sprintf("table%s8_birds_by_crop.docx", L), sprintf("Table %s8: Birds by %s, %s %s", L, what, season, when))
  }
  invisible(summary)
}

# =============================================================================
# 3. RUN EACH SEASON
# =============================================================================
results <- list()
for (season in RUN_SEASONS) {
  cfg <- SEASONS[[season]]
  if (is.null(cfg)) stop("RUN_SEASONS includes '", season, "' but SEASONS has no settings for it.")
  obs <- farm_all %>% filter(month(date) %in% cfg$months) %>%
    mutate(period = period_label(date, season == "winter"))
  if (nrow(obs) == 0) { warning("No farm records in the ", season, " months - skipped.", call. = FALSE); next }
  when <- if (is.na(cfg$focal)) max(obs$period) else cfg$focal
  obs <- obs %>% filter(period == when)
  n_visits <- n_distinct(obs$date)                 # counted before flyovers are removed
  message(sprintf("\n%s %s: %d records on %d visits (%s).", str_to_title(season), when, nrow(obs), n_visits,
                  paste(format(sort(unique(obs$date)), "%d %b"), collapse = ", ")))
  if (cfg$exclude_flyovers) {
    message(sprintf("  Flying-over records left out: %d (%s birds).", sum(obs$flyover), comma(sum(obs$count[obs$flyover]))))
    obs <- obs %>% filter(!flyover)
  }
  if (!all(is.na(cfg$relabel))) {
    hit <- key(obs$crop) %in% key(names(cfg$relabel))
    obs$crop[hit] <- unname(cfg$relabel[match(key(obs$crop[hit]), key(names(cfg$relabel)))])
    message("  Relabelled for ", season, ": ", paste(sprintf("%s -> %s", names(cfg$relabel), cfg$relabel), collapse = "; "),
            " (", sum(hit), " records).")
  }
  message("  ", CROP_COLUMN, " values: ", paste(sort(unique(obs$crop)), collapse = "; "))
  
  areas <- NULL
  if (!is.null(area_raw)) {
    if (!cfg$area_crop_col %in% names(area_raw)) {
      warning("The area file / layer has no column '", cfg$area_crop_col, "' - no ", season, " density figures. Columns available: ",
              paste(setdiff(names(area_raw), ".area_from_shape"), collapse = ", "), call. = FALSE)
    } else {
      areas <- area_raw %>% transmute(k = key(.data[[cfg$area_crop_col]]), area_ha = suppressWarnings(as.numeric(.data[[CROP_AREA_HA_COL]]))) %>%
        filter(!is.na(area_ha), area_ha > 0) %>% group_by(k) %>% summarise(area_ha = sum(area_ha), .groups = "drop")
    }
  }
  results[[season]] <- list(
    crops    = run_set(obs, cfg$crops, areas, n_visits, cfg, season, when, "crop"),
    woodland = if (!all(is.na(cfg$woodland))) run_set(obs, cfg$woodland, areas, n_visits, cfg, season, when, "woodland"))
}

# =============================================================================
# 4. CAPTIONS
# =============================================================================
write_csv(figure_captions, file.path(out_dir, "crop_figure_captions.csv"))
cap_doc <- read_docx() %>%
  body_add_par(paste0("Figure captions — ", SITE_NAME, ": birds by crop and woodland"), style = "heading 1")
for (i in seq_len(nrow(figure_captions))) {
  cap_doc <- cap_doc %>%
    body_add_par(figure_captions$file[i], style = "heading 3") %>%
    body_add_fpar(fpar(ftext(figure_captions$title[i], fp_text(bold = TRUE, font.size = 11)))) %>%
    body_add_par(figure_captions$notes[i], style = "Normal")
}
print(cap_doc, target = file.path(out_dir, "crop_figure_captions.docx"))

message("\nCrop analysis complete. Outputs saved to: ", out_dir)
