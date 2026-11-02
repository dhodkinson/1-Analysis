## ============================================================================
## Kernel-density heat maps of bird records
## CRS: EPSG:27700 (British National Grid)
##
## OUTPUT   PNG maps drawn over a basemap, with optional site-boundary and
##          transect overlays. Set `map_type` in the settings:
##            "abundance" - 1 density map (all species) + 1 per species
##            "richness"  - 1 smoothed species-richness map (all species)
##
## HOW TO USE
##   1. Edit section 2 (files & folders) and run sections 1-2 first. This
##      prints the visit dates available in your data, so you can choose
##      `filter_dates` in step 2.
##   2. Set SEASON ("summer" or "winter") in section 2 - one switch for the
##      visit dates, satellite image and output folder - then edit section 3
##      (SETTINGS). The things you change most often are grouped at the top:
##      `dataset`, `map_type` and the basemap.
##   3. Run the whole script (Source). Everything below section 3 should not
##      need editing.
##
## BASEMAP LICENCES (attribution is added to every map automatically)
##   - Sentinel-2 / OS OpenData (Road, Outdoor, Light) are fine for client
##     reports provided the attribution is shown.
##   - Do not use OS Premium-tier data, Google or Bing imagery in
##     documents given to third parties.
## ============================================================================

rm(list = ls())  # start from a clean session

# ----------------------------------------------------------------------------
# 1. PACKAGES
# ----------------------------------------------------------------------------
## Install once if any are missing:
## install.packages(c("sf", "spatstat.geom", "spatstat.explore", "terra",
##                    "maptiles", "rstac", "httr2", "tidyterra", "ggplot2",
##                    "dplyr", "purrr", "stringr"))
## Optional: install.packages("usethis") - used to open .Renviron for the OS key.
pkgs <- c("sf", "spatstat.geom", "spatstat.explore", "terra", "maptiles",
          "rstac", "httr2", "tidyterra", "ggplot2", "dplyr", "purrr", "stringr")
not_installed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed) > 0) {
  stop("Please install these packages first: ",
       paste(not_installed, collapse = ", "), call. = FALSE)
}
## Load order matters for function masking (terra/dplyr both mask some names),
## so keep this order.
invisible(lapply(pkgs, library, character.only = TRUE))

CRS_DATA <- 27700  # all data, basemaps and maps use British National Grid

# ----------------------------------------------------------------------------
# 2. FILES & FOLDERS  (edit, then run sections 1-2 to see the available dates)
# ----------------------------------------------------------------------------
## Input data: a GeoPackage (.gpkg) or a CSV (.csv). The type is detected
## from the file extension.
INPUT_PATH <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'

## ===== SEASON - THE ONE SWITCH FOR SUMMER / WINTER ===========================
## Set this to "summer" or "winter". It automatically controls:
##   - the visit dates used (each dataset's `filter_dates`, section 3B)
##   - the Sentinel-2 satellite image dates (`sentinel_seasons`, section 3D)
##   - the output folder (".../Analysis results/<Season>/Heatmaps")
##   - each dataset's output subfolder (`output_subdir`, section 3B)
## To add another season, give it an entry in each of those places.
SEASON <- "winter"
## =============================================================================

## Top folder for the project's analysis results.
ANALYSIS_DIR <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Analysis results'

## Maps are saved to <OUTPUT_BASE_DIR>/<output_subdir of the chosen dataset>,
## e.g. ".../Summer/Heatmaps/Counts 2026". The Summer/Winter part follows SEASON.
OUTPUT_BASE_DIR <- file.path(ANALYSIS_DIR, tools::toTitleCase(SEASON), "Heatmaps")

## Where downloaded basemaps are cached (shared by all datasets and seasons).
## Files are named by provider, dates/zoom and extent, so different basemaps
## never clash. NULL = a "_basemap_cache" folder inside OUTPUT_BASE_DIR.
BASEMAP_CACHE_DIR <- file.path(ANALYSIS_DIR, "_basemap_cache")

## Date column used for visit filtering, and how dates are written in the data
## (and in `filter_dates` below).
DATE_COL    <- "Date"
DATE_FORMAT <- "%Y-%m-%d"

## --- Print the dates available in the data ---------------------------------
## Looks at every layer (GeoPackage) or the CSV and lists the dates found, so
## you can copy the ones you want into `filter_dates` in section 3.
list_available_dates <- function(path, date_col, date_format) {
  report <- function(df, label) {
    if (!date_col %in% names(df)) {
      message("'", label, "': no '", date_col, "' column (nothing to filter by)")
      return(invisible(NULL))
    }
    dates <- sort(unique(as.Date(as.character(df[[date_col]]), format = date_format)))
    message("Available dates in '", label, "':")
    print(dates)
  }
  if (tolower(tools::file_ext(path)) == "gpkg") {
    for (layer in sf::st_layers(path)$name) {
      report(sf::st_drop_geometry(sf::st_read(path, layer = layer, quiet = TRUE)), layer)
    }
  } else {
    report(read.csv(path, stringsAsFactors = FALSE), basename(path))
  }
}
list_available_dates(INPUT_PATH, DATE_COL, DATE_FORMAT)

# ----------------------------------------------------------------------------
# 3. SETTINGS - EDIT THESE
# ----------------------------------------------------------------------------
CONFIG <- list(
  
  ## ===== A. CHANGE EVERY RUN ================================================
  
  ## Which dataset preset (from `datasets` below) to map.
  dataset = "counts",   # "counts" or "territories"
  ## (summer vs winter is set once, by SEASON in section 2)
  
  ## What to map:
  ##   "abundance" - density of records (weighted by the preset's weight_col);
  ##                 one map for all species + one per species
  ##   "richness"  - smoothed species richness (how many different species
  ##                 occur around each point); one map for all species
  map_type = "abundance",
  
  ## Basemap. One of:
  ##   "Sentinel2"                          free Copernicus satellite imagery
  ##   "OS.Road", "OS.Outdoor", "OS.Light"  OS Data Hub styles (need an API key)
  ##   any maptiles provider name, e.g. "OpenStreetMap", "CartoDB.Positron",
  ##   "Esri.WorldImagery", "OpenTopoMap"
  ## (OS "Leisure" is not supported: it is only published in British National
  ## Grid, which the Web-Mercator tile method used here cannot read.)
  tile_provider   = "Sentinel2",
  sentinel_season = SEASON,     # follows SEASON (section 2); entries are in
  # basemap$sentinel_seasons below
  
  ## ===== B. DATASET PRESETS ================================================
  ## Each preset bundles the settings that always change together.
  ##   layer_name    GeoPackage layer holding the records
  ##   weight_col    column to weight density by (NULL = unweighted)
  ##   filter_dates  visit dates to KEEP, as text. NULL = keep all dates.
  ##   output_subdir folder (inside OUTPUT_BASE_DIR) the maps are saved to;
  ##                 NULL = save directly in OUTPUT_BASE_DIR
  ## filter_dates and output_subdir can be either one value (used for every
  ## season) or a list with one entry per season, picked by SEASON, e.g.
  ##   filter_dates = list(summer = c("2026-04-14"), winter = c("2025-12-20"))
  ## Add more presets if needed.
  datasets = list(
    counts = list(
      layer_name    = '2026_Counts_final',
      weight_col    = "Count",
      filter_dates  = list(summer = c("2026-04-14", "2026-05-30"),
                           winter = c("2025-12-20", "2026-02-09")),
      output_subdir = list(summer = "Counts 2026",
                           winter = NULL)   # NULL = straight into OUTPUT_BASE_DIR
    ),
    territories = list(
      layer_name    = '2026_Territories_final',
      weight_col    = NULL,   # one row = one territory, so no weighting
      filter_dates  = NULL,   # this layer has no Date column
      output_subdir = "Territories 2026"
    )
  ),
  
  ## Column names (x_col / y_col are only used for CSV input).
  species_col = "Species",
  x_col       = "x_coord",
  y_col       = "y_coord",
  
  ## ===== C. KERNEL DENSITY ==================================================
  kde = list(
    grid_res_m    = 10,  # output cell size (m). Smaller = smoother but slower;
    # farm-scale maps look best at roughly 10-50 m
    sigma_m       = 50,  # bandwidth (m): how far each record's influence
    # spreads. NULL = automatic (bw.diggle)
    buffer_m      = 25,   # padding (m) around the data/overlays extent
    min_points    = 0,    # skip species with fewer records than this
    mask_quantile = 0     # cells below this quantile of density (0-1) become
    # transparent so the basemap shows through
  ),
  
  ## ===== C2. SPECIES RICHNESS (only used when map_type = "richness") ========
  ## Each species gets a "presence" surface: 1 at a record, fading to 0 over
  ## roughly kde$sigma_m metres (abundance is ignored - ten birds count the
  ## same as one). Summing these across species gives an approximate count of
  ## species present within about sigma_m of each location.
  richness = list(
    min_value   = 0.00,  # cells below this are transparent (basemap shows)
    show_legend = FALSE,  # TRUE = colour bar in the bottom-left corner
    legend_name = "Species richness",
    by_date     = FALSE   # TRUE = also one map per visit date (needs a date column)
  ),
  
  ## ===== D. BASEMAP DETAILS ================================================
  basemap = list(
    tile_zoom   = 14,     # tile zoom for non-Sentinel basemaps (higher = sharper)
    use_cache   = TRUE,   # reuse a cached basemap if one exists. Set FALSE once
    # (or delete files in BASEMAP_CACHE_DIR) to refetch
    os_api_key  = Sys.getenv("OS_DATA_KEY"),  # read from .Renviron, never
    # type the key into this script
    ## Sentinel-2: named date-range presets. Each is cached separately, so once
    ## fetched, switching between them is instant.
    sentinel_seasons = list(
      summer = list(start = "2026-04-01", end = "2026-06-30"),
      winter = list(start = "2025-12-01", end = "2026-02-28")
    ),
    sentinel_max_cloud_cover = 20   # %; the least-cloudy scene under this is used.
    # Raise it if no scenes are found.
  ),
  
  ## ===== E. HEAT-MAP APPEARANCE ============================================
  appearance = list(
    heat_alpha = 0.6,   # opacity of the density layer (lower = more basemap)
    ## Colour ramp. Either:
    ##  - a vector of 2+ colours, low -> high density, e.g.
    ##      c("white", "yellow", "red")      classic heat
    ##      c("white", "orange", "red")      softer
    ##      c("white", "gold", "firebrick")  warm, good in print
    ##      c("white", "pink", "red")        gentler
    ##      c("white", "red", "darkred")     restrained
    ##  - NULL to use the viridis scale named in `palette`
    ##    ("viridis", "magma", "inferno", "plasma", "turbo", ...)
    palette_colors = c("white", "orange", "red", "darkred"),
    palette        = "viridis",
    png_width_in   = 10,   # saved image width (height follows the map shape)
    png_dpi        = 300,
    save_pdf       = FALSE # TRUE = also save a PDF of every map
  ),
  
  ## ===== F. OVERLAYS (drawn from layers in the same GeoPackage) =============
  ## Colours can be names or hex codes; linetype is "solid", "dashed",
  ## "dotted", "dotdash"...; alpha is 1 (opaque) to 0 (invisible).
  ## To list layer names: sf::st_layers(INPUT_PATH)
  boundary = list(
    enabled    = TRUE,
    layer      = "Field_polygons",
    color      = "black",
    linewidth  = 0.3,
    linetype   = "solid",
    alpha      = 0.5,
    fill       = NA,     # e.g. "grey50" for a tint, NA for outline only
    fill_alpha = 0.15    # ignored when fill = NA
  ),
  transect = list(
    enabled    = TRUE,
    layer      = "survey_route",
    color      = "red",
    linewidth  = 0.8,
    linetype   = "solid",
    alpha      = 1
  )
)

## ONE-TIME SETUP for an OS Data Hub key (only needed for OS.* basemaps):
##   1. Run usethis::edit_r_environ()
##   2. Add the line   OS_DATA_KEY=your_actual_key_here
##   3. Save, then restart R (Session > Restart R).
## The key then loads automatically and never appears in this script.

# ============================================================================
# EVERYTHING BELOW HERE SHOULD NOT NEED EDITING
# ============================================================================

# ----------------------------------------------------------------------------
# 4. HELPER FUNCTIONS
# ----------------------------------------------------------------------------

## ---- 4a. Reading data ------------------------------------------------------

## Read one layer from a GeoPackage and make sure it is in EPSG:27700.
## (A layer with no CRS is assumed to already be in 27700.)
read_sf_layer <- function(path, layer) {
  x <- sf::st_read(path, layer = layer, quiet = TRUE)
  if (is.na(sf::st_crs(x))) {
    sf::st_crs(x) <- CRS_DATA
  } else if (!isTRUE(sf::st_crs(x)$epsg == CRS_DATA)) {
    x <- sf::st_transform(x, CRS_DATA)
  }
  x
}

## Load the bird records (GeoPackage layer or CSV with x/y columns) as sf points.
load_points <- function(path, layer, x_col, y_col, species_col) {
  if (tolower(tools::file_ext(path)) == "gpkg") {
    pts <- read_sf_layer(path, layer)
  } else {
    dat <- read.csv(path, stringsAsFactors = FALSE)
    stopifnot(all(c(x_col, y_col, species_col) %in% names(dat)))
    dat <- dplyr::filter(dat, !is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
    pts <- sf::st_as_sf(dat, coords = c(x_col, y_col), crs = CRS_DATA)
  }
  if (!species_col %in% names(pts)) {
    stop("Species column '", species_col, "' not found. Columns are: ",
         paste(names(pts), collapse = ", "), call. = FALSE)
  }
  pts
}

## Load an optional overlay (boundary / transect) from the GeoPackage.
## Returns NULL (with a warning) if it is disabled, impossible or missing,
## so a typo in a layer name never stops the whole run.
load_overlay <- function(spec, label, path) {
  if (!isTRUE(spec$enabled)) return(NULL)
  if (tolower(tools::file_ext(path)) != "gpkg") {
    warning("The ", label, " overlay needs a GeoPackage input - skipping it.",
            call. = FALSE)
    return(NULL)
  }
  available <- sf::st_layers(path)$name
  if (!spec$layer %in% available) {
    warning(sprintf(
      "%s layer '%s' not found in %s. Available layers: %s. Skipping it - fix the layer name or set enabled = FALSE.",
      label, spec$layer, basename(path), paste(available, collapse = ", ")
    ), call. = FALSE)
    return(NULL)
  }
  message("Loaded ", label, " layer '", spec$layer, "'")
  read_sf_layer(path, spec$layer)
}

## Pick the value for the current season from a dataset setting. If `x` is a
## named list (e.g. list(summer = ..., winter = ...)) return the entry for
## `season`; otherwise `x` applies to every season and is returned unchanged.
by_season <- function(x, season, what) {
  if (is.list(x) && !is.null(names(x))) {
    if (!season %in% names(x)) {
      stop("'", what, "' has no entry for SEASON = '", season, "' (it has: ",
           paste(names(x), collapse = ", "), ").", call. = FALSE)
    }
    return(x[[season]])
  }
  x
}

## Keep only records whose date is in `keep` (text dates). NULL = keep all.
filter_by_dates <- function(pts, keep, date_col, date_format) {
  if (is.null(keep)) return(pts)
  if (!date_col %in% names(pts)) {
    warning("filter_dates is set but this layer has no '", date_col,
            "' column - skipping the date filter.", call. = FALSE)
    return(pts)
  }
  record_dates <- as.Date(as.character(pts[[date_col]]), format = date_format)
  kept <- pts[record_dates %in% as.Date(keep, format = date_format), ]
  message(sprintf("Date filter: kept %d of %d records (dates: %s)",
                  nrow(kept), nrow(pts), paste(keep, collapse = ", ")))
  kept
}

## ---- 4b. Map extent ------------------------------------------------------

## One bounding box covering every supplied layer (NULLs ignored), padded by
## `buffer_m`. Using the points alone could clip parts of the boundary or
## route that extend beyond where records happen to fall.
study_bbox <- function(layers, buffer_m) {
  layers <- Filter(Negate(is.null), layers)
  bb <- vapply(layers, function(l) as.numeric(sf::st_bbox(l)), numeric(4))
  sf::st_bbox(c(xmin = min(bb[1, ]) - buffer_m, ymin = min(bb[2, ]) - buffer_m,
                xmax = max(bb[3, ]) + buffer_m, ymax = max(bb[4, ]) + buffer_m),
              crs = sf::st_crs(CRS_DATA))
}

## ---- 4c. Basemap ---------------------------------------------------------

## Providers that go through the OS Maps API (Web Mercator "_3857" layers).
OS_LAYERS <- c(OS.Road = "Road_3857", OS.Outdoor = "Outdoor_3857", OS.Light = "Light_3857")

## Stop early, with a fix-it prompt, if an OS basemap is chosen but no API key is set.
check_os_key <- function(provider, key) {
  if (!provider %in% names(OS_LAYERS) || nzchar(key)) return(invisible(TRUE))
  msg <- paste0(
    "tile_provider = '", provider, "' needs an OS Data Hub API key, but OS_DATA_KEY is empty.\n",
    "Fix: add the line  OS_DATA_KEY=your_actual_key_here  to your .Renviron, ",
    "then restart R and run again."
  )
  if (interactive() && requireNamespace("usethis", quietly = TRUE)) {
    message(msg)
    if (isTRUE(utils::askYesNo("Open your .Renviron now?"))) usethis::edit_r_environ()
    stop("Stopped: add OS_DATA_KEY to .Renviron, restart R, then re-run.", call. = FALSE)
  }
  stop(msg, call. = FALSE)
}

## Sentinel-2: fetch ESA's pre-rendered true-colour image of the least-cloudy
## scene in the date range, via Microsoft's free Planetary Computer catalogue
## (no key needed). Returns an RGB raster cropped to the area of interest.
fetch_sentinel2 <- function(aoi, start_date, end_date, max_cloud_cover) {
  bbox_wgs84 <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))
  ## Full RFC 3339 timestamps avoid rstac's date-parsing errors.
  datetime <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")
  
  items <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
    rstac::stac_search(collections = "sentinel-2-l2a", bbox = bbox_wgs84,
                       datetime = datetime, limit = 100) |>
    rstac::post_request()
  if (length(items$features) == 0) {
    stop("No Sentinel-2 scenes found for this area and date range.", call. = FALSE)
  }
  
  cloud <- vapply(items$features, function(f) f$properties[["eo:cloud_cover"]], numeric(1))
  ok <- cloud <= max_cloud_cover
  if (!any(ok)) {
    stop(sprintf(paste("No Sentinel-2 scenes under %.0f%% cloud (lowest available: %.0f%%).",
                       "Raise basemap$sentinel_max_cloud_cover or widen the date range."),
                 max_cloud_cover, min(cloud)), call. = FALSE)
  }
  best <- which(ok)[which.min(cloud[ok])]
  item <- items$features[[best]]
  message("Using Sentinel-2 scene from ", item$properties$datetime,
          " (", round(cloud[best], 1), "% cloud cover)")
  
  ## Image URLs must be "signed" via a public endpoint; the signature expires
  ## after about an hour, which is why the result is cached to disk below.
  signed <- httr2::request("https://planetarycomputer.microsoft.com/api/sas/v1/sign") |>
    httr2::req_url_query(href = item$assets$visual$href) |>
    httr2::req_perform() |>
    httr2::resp_body_json()
  
  r <- terra::rast(signed$href)
  r <- terra::project(r, paste0("EPSG:", CRS_DATA))
  r <- terra::crop(r, terra::vect(aoi))
  ## `* 1` forces every pixel to be read into memory NOW. Otherwise terra keeps
  ## a lazy link to the expiring signed URL and re-fetches (then fails) later.
  r * 1
}

## Standard slippy-map tiles (OS Data Hub styles or any maptiles provider),
## reprojected to British National Grid.
fetch_tiles <- function(provider, aoi, zoom, os_key) {
  if (provider %in% names(OS_LAYERS)) {
    os_layer <- OS_LAYERS[[provider]]
    tile_source <- maptiles::create_provider(
      name     = paste0("OS_Maps_API_", os_layer),
      url      = paste0("https://api.os.uk/maps/raster/v1/zxy/", os_layer,
                        "/{z}/{x}/{y}.png?key={apikey}"),
      citation = "Contains OS data (c) Crown copyright and database rights"
    )
    tiles <- maptiles::get_tiles(aoi, provider = tile_source, zoom = zoom,
                                 crop = TRUE, apikey = os_key)
  } else {
    tiles <- maptiles::get_tiles(aoi, provider = provider, zoom = zoom, crop = TRUE)
  }
  terra::project(tiles, paste0("EPSG:", CRS_DATA))
}

## Get the basemap for the run, using the disk cache when possible.
## The cache file name encodes provider + (dates or zoom) + rounded extent, so
## changing any of those automatically fetches a fresh image.
get_basemap <- function(cfg, bbox, cache_dir) {
  bm       <- cfg$basemap
  provider <- cfg$tile_provider
  is_s2    <- identical(provider, "Sentinel2")
  aoi      <- sf::st_as_sfc(bbox)
  
  if (is_s2) {
    if (!cfg$sentinel_season %in% names(bm$sentinel_seasons)) {
      stop("sentinel_season '", cfg$sentinel_season, "' is not one of: ",
           paste(names(bm$sentinel_seasons), collapse = ", "), call. = FALSE)
    }
    season <- bm$sentinel_seasons[[cfg$sentinel_season]]
    message(sprintf("Sentinel-2 season: '%s' (%s to %s)",
                    cfg$sentinel_season, season$start, season$end))
    detail <- c(season$start, season$end)
  } else {
    if (identical(provider, "OS.Leisure")) {
      stop("'OS.Leisure' is not supported (it is only published in British National ",
           "Grid). Use OS.Road, OS.Outdoor or OS.Light.", call. = FALSE)
    }
    check_os_key(provider, bm$os_api_key)
    detail <- bm$tile_zoom
  }
  
  key <- paste(c(provider, detail, round(as.numeric(bbox))),
               collapse = "_")
  cache_path <- file.path(cache_dir, paste0("basemap_", gsub("[^A-Za-z0-9_]", "", key), ".tif"))
  
  if (isTRUE(bm$use_cache) && file.exists(cache_path)) {
    message("Using cached basemap: ", cache_path)
    return(terra::rast(cache_path))
  }
  
  message("Fetching basemap...")
  basemap <- if (is_s2) {
    fetch_sentinel2(aoi, season$start, season$end, bm$sentinel_max_cloud_cover)
  } else {
    fetch_tiles(provider, aoi, bm$tile_zoom, bm$os_api_key)
  }
  if (isTRUE(bm$use_cache)) {
    terra::writeRaster(basemap, cache_path, overwrite = TRUE)
    message("Cached basemap saved to: ", cache_path)
  }
  basemap
}

## Credit line printed on every map, as the data licences require.
basemap_attribution <- function(provider) {
  year <- format(Sys.Date(), "%Y")
  if (identical(provider, "Sentinel2")) {
    sprintf("Contains modified Copernicus Sentinel data %s", year)
  } else if (provider %in% names(OS_LAYERS)) {
    sprintf("Contains OS data © Crown copyright and database rights %s", year)
  } else if (grepl("^OpenStreetMap", provider)) {
    "© OpenStreetMap contributors"
  } else if (grepl("^Esri", provider)) {
    "Basemap © Esri"
  } else if (grepl("^CartoDB", provider)) {
    "© OpenStreetMap contributors © CARTO"
  } else if (grepl("^Stadia", provider)) {
    "© Stadia Maps © OpenStreetMap contributors"
  } else if (grepl("^OpenTopoMap", provider)) {
    "© OpenStreetMap contributors, SRTM | © OpenTopoMap (CC-BY-SA)"
  } else {
    ""   # unknown provider: add the right credit yourself
  }
}

## ---- 4d. Kernel density ----------------------------------------------------

## Kernel-density raster of a set of points over `window`.
## If `weight_col` is given, each point counts as its value in that column
## (e.g. number of birds); otherwise every point counts as 1.
compute_kde <- function(pts, window, kde, weight_col) {
  coords <- sf::st_coordinates(pts)
  ppp_obj <- spatstat.geom::ppp(coords[, 1], coords[, 2], window = window, checkdup = FALSE)
  
  weights <- NULL
  if (!is.null(weight_col) && weight_col %in% names(pts)) {
    weights <- pts[[weight_col]]
    weights[is.na(weights)] <- 0
  }
  
  sigma <- if (is.null(kde$sigma_m)) spatstat.explore::bw.diggle(ppp_obj) else kde$sigma_m
  dens  <- spatstat.explore::density.ppp(ppp_obj, sigma = sigma, eps = kde$grid_res_m,
                                         weights = weights, positive = TRUE)
  
  r <- terra::rast(dens)
  terra::crs(r) <- paste0("EPSG:", CRS_DATA)
  
  ## Hide the faintest cells so the basemap shows cleanly through them.
  if (!is.null(kde$mask_quantile) && kde$mask_quantile > 0) {
    cutoff <- stats::quantile(terra::values(r), probs = kde$mask_quantile, na.rm = TRUE)
    r[r < cutoff] <- NA
  }
  r
}

## Approximate species richness surface over `window`.
## For each species: sum of unit-peak Gaussian bumps (one per record, width
## sigma), capped at 1 -> a 0-1 "presence" surface that does not depend on how
## many individuals were recorded. Adding the species' surfaces together gives
## roughly "number of species present within ~sigma of this spot".
compute_richness <- function(pts, window, kde, species_col, min_value) {
  sigma <- if (is.null(kde$sigma_m)) {
    all_xy <- sf::st_coordinates(pts)
    spatstat.explore::bw.diggle(spatstat.geom::ppp(all_xy[, 1], all_xy[, 2],
                                                   window = window, checkdup = FALSE))
  } else kde$sigma_m
  
  presence <- lapply(split(pts, pts[[species_col]]), function(sp_pts) {
    xy <- sf::st_coordinates(sp_pts)
    p  <- spatstat.geom::ppp(xy[, 1], xy[, 2], window = window, checkdup = FALSE)
    ## edge = FALSE keeps the peak calibrated; x 2*pi*sigma^2 rescales the
    ## density so a single record peaks at exactly 1.
    d <- spatstat.explore::density.ppp(p, sigma = sigma, eps = kde$grid_res_m,
                                       edge = FALSE, positive = TRUE)
    terra::rast(spatstat.geom::eval.im(pmin(1, d * 2 * pi * sigma^2)))
  })
  
  r <- Reduce(`+`, presence)
  terra::crs(r) <- paste0("EPSG:", CRS_DATA)
  r[r < min_value] <- NA
  r
}

## ---- 4e. Drawing & saving ------------------------------------------------

## Colour scale for the density layer. NA cells stay fully transparent.
heat_scale <- function(app, name = "Density") {
  if (!is.null(app$palette_colors)) {
    scale_fill_gradientn(colors = app$palette_colors, name = name, na.value = NA)
  } else {
    scale_fill_viridis_c(option = app$palette, name = name, na.value = NA)
  }
}

## Add a boundary/transect layer to a plot. An optional fill is drawn first
## (with its own opacity), then the outline on top, so the two can be styled
## independently. `spec$fill` may be missing (lines have no fill).
add_overlay <- function(p, data, spec) {
  if (is.null(data)) return(p)
  if (!is.null(spec$fill) && !is.na(spec$fill)) {
    p <- p + geom_sf(data = data, fill = spec$fill, color = NA,
                     alpha = spec$fill_alpha, inherit.aes = FALSE)
  }
  p + geom_sf(data = data, fill = NA, color = spec$color, linewidth = spec$linewidth,
              linetype = spec$linetype, alpha = spec$alpha, inherit.aes = FALSE)
}

## Draw a surface raster (density or richness) over the basemap and save it.
## `ctx` carries everything shared by all maps (basemap, extent, overlays,
## attribution) so it is prepared only once.
save_map <- function(surface, file_stub, ctx, cfg, out_dir,
                     legend_name = "Density", show_legend = FALSE) {
  app <- cfg$appearance
  bb  <- ctx$bbox
  
  p <- ggplot() +
    tidyterra::geom_spatraster_rgb(data = ctx$basemap, interpolate = TRUE) +
    tidyterra::geom_spatraster(data = surface, alpha = app$heat_alpha, interpolate = TRUE) +
    heat_scale(app, legend_name)
  
  p <- add_overlay(p, ctx$boundary, cfg$boundary)
  p <- add_overlay(p, ctx$transect, cfg$transect)   # on top of the boundary
  
  ## Attribution: plain text just inside the bottom-right corner of the map.
  inset_x <- (bb[["xmax"]] - bb[["xmin"]]) * 0.01
  inset_y <- (bb[["ymax"]] - bb[["ymin"]]) * 0.01
  
  p <- p +
    coord_sf(crs = CRS_DATA, expand = FALSE) +
    annotate("text", x = bb[["xmax"]] - inset_x, y = bb[["ymin"]] + inset_y,
             label = ctx$attribution, hjust = 1, vjust = 0, size = 2.3, color = "grey20") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          ## Legend (if shown) sits inside the map, bottom-left, on a pale panel.
          legend.position = if (show_legend) "inside" else "none",
          legend.position.inside = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.background = element_rect(fill = "#FFFFFFB3", colour = NA),
          plot.margin = margin(0, 0, 0, 0))
  
  ## Size the image to the map's own aspect ratio so there is no white space.
  width  <- app$png_width_in
  height <- width * (bb[["ymax"]] - bb[["ymin"]]) / (bb[["xmax"]] - bb[["xmin"]])
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), plot = p,
         width = width, height = height, dpi = app$png_dpi, bg = "white")
  if (isTRUE(app$save_pdf)) {
    ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), plot = p,
           width = width, height = height, bg = "white")
  }
  message("Saved: ", file_stub)
  invisible(p)
}

## Abundance map: kernel density of the records (weighted if weight_col is set).
make_heatmap <- function(pts, name, file_stub, ctx, cfg, weight_col, out_dir) {
  if (nrow(pts) < cfg$kde$min_points) {
    message(sprintf("Skipping '%s': %d record(s), fewer than min_points (%d)",
                    name, nrow(pts), cfg$kde$min_points))
    return(invisible(NULL))
  }
  surface <- compute_kde(pts, ctx$window, cfg$kde, weight_col)
  save_map(surface, file_stub, ctx, cfg, out_dir)
}

## Species-richness map for a set of records (all species together).
make_richness_map <- function(pts, file_stub, ctx, cfg, out_dir) {
  rc <- cfg$richness
  surface <- compute_richness(pts, ctx$window, cfg$kde, cfg$species_col, rc$min_value)
  save_map(surface, file_stub, ctx, cfg, out_dir,
           legend_name = rc$legend_name, show_legend = rc$show_legend)
}

# ----------------------------------------------------------------------------
# 5. RUN
# ----------------------------------------------------------------------------

## ---- 5a. Resolve the chosen dataset preset ---------------------------------
if (!CONFIG$dataset %in% names(CONFIG$datasets)) {
  stop("CONFIG$dataset ('", CONFIG$dataset, "') must be one of: ",
       paste(names(CONFIG$datasets), collapse = ", "), call. = FALSE)
}
ds <- CONFIG$datasets[[CONFIG$dataset]]

## Settings that can differ by season (filter_dates, output_subdir) -> this season's value.
filter_dates <- by_season(ds$filter_dates,  SEASON, paste0(CONFIG$dataset, "$filter_dates"))
output_subdir <- by_season(ds$output_subdir, SEASON, paste0(CONFIG$dataset, "$output_subdir"))

## output_subdir = NULL means "no subfolder": save directly in OUTPUT_BASE_DIR.
out_dir   <- if (is.null(output_subdir)) OUTPUT_BASE_DIR else file.path(OUTPUT_BASE_DIR, output_subdir)
cache_dir <- if (is.null(BASEMAP_CACHE_DIR)) file.path(OUTPUT_BASE_DIR, "_basemap_cache") else BASEMAP_CACHE_DIR
dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("Season: '%s' | Dataset: '%s' (layer '%s', weighted by %s)", SEASON,
                CONFIG$dataset, ds$layer_name,
                if (is.null(ds$weight_col)) "nothing" else ds$weight_col))
message("Output folder: ", out_dir)

## ---- 5b. Load data, overlays and filter by date ---------------------------
pts_sf <- load_points(INPUT_PATH, ds$layer_name, CONFIG$x_col, CONFIG$y_col, CONFIG$species_col)
pts_sf <- filter_by_dates(pts_sf, filter_dates, DATE_COL, DATE_FORMAT)
if (nrow(pts_sf) == 0) stop("No records left after filtering - check filter_dates.", call. = FALSE)

boundary_sf <- load_overlay(CONFIG$boundary, "boundary", INPUT_PATH)
transect_sf <- load_overlay(CONFIG$transect, "transect", INPUT_PATH)

## ---- 5c. Extent and basemap (shared by every map) -------------------------
bbox <- study_bbox(list(pts_sf, boundary_sf, transect_sf), CONFIG$kde$buffer_m)

ctx <- list(
  bbox        = bbox,
  window      = spatstat.geom::as.owin(bbox),   # KDE study window = map extent
  basemap     = get_basemap(CONFIG, bbox, cache_dir),
  boundary    = boundary_sf,
  transect    = transect_sf,
  attribution = basemap_attribution(CONFIG$tile_provider)
)

## ---- 5d. Make the maps ------------------------------------------------------
if (!CONFIG$map_type %in% c("abundance", "richness")) {
  stop("CONFIG$map_type must be 'abundance' or 'richness'", call. = FALSE)
}

if (CONFIG$map_type == "abundance") {
  ## Overview map, then one map per species.
  make_heatmap(pts_sf, "All species", "heatmap_ALL_SPECIES", ctx, CONFIG, ds$weight_col, out_dir)
  
  by_species <- split(pts_sf, pts_sf[[CONFIG$species_col]])
  purrr::iwalk(by_species, function(sp_pts, sp_name) {
    stub <- paste0("heatmap_", stringr::str_replace_all(sp_name, "[^A-Za-z0-9]+", "_"))
    make_heatmap(sp_pts, sp_name, stub, ctx, CONFIG, ds$weight_col, out_dir)
  })
  
} else {
  ## Richness: one map for all species (a single species has no "richness").
  make_richness_map(pts_sf, "richness_ALL_SPECIES", ctx, CONFIG, out_dir)
  
  ## Optionally also one map per visit date.
  if (isTRUE(CONFIG$richness$by_date)) {
    if (!DATE_COL %in% names(pts_sf)) {
      warning("richness$by_date = TRUE but there is no '", DATE_COL,
              "' column in this dataset - skipping per-date maps.", call. = FALSE)
    } else {
      visit_date <- as.character(as.Date(as.character(pts_sf[[DATE_COL]]), format = DATE_FORMAT))
      purrr::iwalk(split(pts_sf, visit_date), function(d_pts, d_name) {
        make_richness_map(d_pts, paste0("richness_", d_name), ctx, CONFIG, out_dir)
      })
    }
  }
}

message("Done. Outputs written to: ", normalizePath(out_dir))