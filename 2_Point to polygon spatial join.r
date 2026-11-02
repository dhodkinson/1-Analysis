## =================================================================
## Attach crop/habitat and field code to bird records (spatial join)
## =================================================================
##
## WHAT THIS SCRIPT DOES
## ----------------------
## For each input point layer (counts and territories), finds the
## polygon in the crop layer that each record falls inside and copies
## its Habitat and field code onto the record. The results are written
## back to the GeoPackage and to an attribute-only CSV.
##
## Field codes now come from the crop layer itself, so there is a single
## join. Field_polygons is no longer used here (it is still used by the
## heatmap script for drawing field boundaries).
##
## CHECKS (the script stops or warns rather than writing bad data)
##   - invalid polygon geometries are reported and repaired
##   - overlapping crop polygons are listed (a point in an overlap would
##     otherwise be duplicated by the join)
##   - STOPS if any record is duplicated by the join
##   - reports records that fall outside every polygon (no field code)
##   - warns where a record's field code does not match its transect
##     (A-L should be Shimpling Park Farm, M-X Lavenham Airfield)
## =================================================================

rm(list = ls()) # NEVER REMOVE

library(sf)
library(dplyr)

## -----------------------------------------------------------------
## 0. CONFIG - edit these
## -----------------------------------------------------------------
gpkg_path <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/GIS/Shimpling.gpkg'
csv_dir   <- '/Users/david/Library/CloudStorage/Dropbox/Landmark Ecology/Projects/2506 Shimpling Park Farm repeat/Data/'

year_prefix <- "2026"

# Input point layers to process, named by the label used in the output
# names: "<year>_<label>_final" (GeoPackage layer and CSV).
input_layers <- c(
  Territories = "2026_Territories_clean",
  Counts      = "2026_Counts_clean"
)

crop_layer <- "2026_Crop_polygons"

# Columns to take from the crop layer. The field code column may be
# called "Field Code", "Field.Code", "Field_Code" or "FieldCode"; it is
# found automatically and always written out as "Field.Code" so the
# later scripts (3-8) need no changes.
habitat_col <- "Habitat"

# Which field codes belong to which transect (used for the site check)
site_codes <- list(
  "Shimpling Park Farm" = LETTERS[1:12],   # A-L
  "Lavenham Airfield"   = LETTERS[13:24]   # M-X
)

# Overlaps smaller than this (m2) are treated as digitising slivers and
# only listed together in the summary, not one by one.
overlap_report_min_m2 <- 1

## -----------------------------------------------------------------
## 1. Read and check the crop polygons
## -----------------------------------------------------------------
polys <- st_read(gpkg_path, layer = crop_layer, quiet = TRUE, fid_column_name = "fid")

fc_candidates <- c("Field Code", "Field.Code", "Field_Code", "FieldCode", "field_code", "Field code")
fc_col <- intersect(fc_candidates, names(polys))[1]
if (is.na(fc_col)) {
  stop("No field code column found in '", crop_layer, "'. Columns are: ",
       paste(setdiff(names(polys), attr(polys, "sf_column")), collapse = ", "))
}
if (!habitat_col %in% names(polys)) {
  stop("No '", habitat_col, "' column found in '", crop_layer, "'.")
}
message("Using '", fc_col, "' as the field code column in ", crop_layer, ".")

# Keep only the two columns we need, tidy the text, and give each polygon an id
polys <- polys %>%
  mutate(
    poly_id    = fid,   # the polygon's fid, as shown in QGIS
    Habitat    = trimws(as.character(.data[[habitat_col]])),
    Field.Code = toupper(trimws(as.character(.data[[fc_col]]))),
    Habitat    = na_if(Habitat, ""),
    Field.Code = na_if(Field.Code, "")
  ) %>%
  select(poly_id, Habitat, Field.Code)

# Invalid geometries
invalid <- which(!st_is_valid(polys))
if (length(invalid) > 0) {
  message("Repairing ", length(invalid), " invalid polygon(s): poly_id ",
          paste(polys$poly_id[invalid], collapse = ", "),
          " (", paste(polys$Habitat[invalid], collapse = "; "), ")")
  polys <- st_make_valid(polys)
}

# Overlapping polygons
ov <- suppressWarnings(st_intersection(polys %>% select(poly_id), polys %>% select(poly_id_b = poly_id))) %>%
  filter(poly_id < poly_id_b) %>%
  mutate(area_m2 = as.numeric(st_area(.))) %>%
  filter(area_m2 > 0)
if (nrow(ov) > 0) {
  ov_tab <- ov %>% st_drop_geometry() %>%
    left_join(st_drop_geometry(polys), by = "poly_id") %>%
    left_join(st_drop_geometry(polys) %>% rename(poly_id_b = poly_id, Habitat_b = Habitat, Field.Code_b = Field.Code),
              by = "poly_id_b") %>%
    arrange(desc(area_m2))
  big <- ov_tab %>% filter(area_m2 >= overlap_report_min_m2)
  message(nrow(ov_tab), " overlapping polygon pair(s) in ", crop_layer,
          " (", nrow(big), " of ", overlap_report_min_m2, " m2 or more):")
  if (nrow(big) > 0) print(as.data.frame(big %>% mutate(area_m2 = round(area_m2, 1))), row.names = FALSE)
  message("Points inside an overlap would be duplicated; the check after the join will stop the script if that happens.")
} else {
  message("No overlapping polygons in ", crop_layer, ".")
}

## -----------------------------------------------------------------
## 2. Join and check every point layer (nothing is written until all pass)
## -----------------------------------------------------------------
results <- list()
for (label in names(input_layers)) {
  input_layer  <- input_layers[[label]]
  output_layer <- paste0(year_prefix, "_", label, "_final")
  csv_path     <- file.path(csv_dir, paste0(year_prefix, "_", label, "_final.csv"))
  message("\n==== ", input_layer, " -> ", output_layer, " ====")
  
  points <- st_read(gpkg_path, layer = input_layer, quiet = TRUE)
  if (!"record_id" %in% names(points)) stop("'", input_layer, "' has no record_id column.")
  
  # Drop any Habitat / field code columns left over from an earlier join
  points <- points %>% select(-any_of(c("Habitat", "Field.Code", "Field Code", "Field_Code")))
  
  # Make sure both layers use the same coordinate system
  polys_use <- if (st_crs(polys) != st_crs(points)) st_transform(polys, st_crs(points)) else polys
  
  joined <- st_join(points, polys_use, join = st_within)
  
  # --- Check: no record may be duplicated by the join ---
  dups <- joined %>% st_drop_geometry() %>% count(record_id) %>% filter(n > 1)
  if (nrow(dups) > 0) {
    dup_detail <- joined %>% st_drop_geometry() %>%
      filter(record_id %in% dups$record_id) %>%
      select(any_of(c("record_id", "Species", "Site", "Date")), poly_id, Habitat, Field.Code)
    print(as.data.frame(dup_detail), row.names = FALSE)
    stop(nrow(dups), " record(s) in '", input_layer, "' fall inside more than one crop polygon ",
         "(listed above; poly_id is the polygon's fid in QGIS). Fix the overlap(s) in ", crop_layer,
         " and re-run. Nothing has been written for any layer.")
  }
  
  # --- Report: records outside every polygon ---
  outside <- joined %>% st_drop_geometry() %>% filter(is.na(poly_id))
  if (nrow(outside) > 0) {
    message(nrow(outside), " record(s) fall outside every crop polygon (no Habitat or field code):")
    print(as.data.frame(outside %>% count(Site, name = "records")), row.names = FALSE)
  }
  
  # --- Report: inside a polygon that has no field code (e.g. neighbouring land) ---
  no_code <- joined %>% st_drop_geometry() %>% filter(!is.na(poly_id), is.na(Field.Code))
  if (nrow(no_code) > 0) {
    message(nrow(no_code), " record(s) are inside a crop polygon with no field code (excluded from field tables):")
    print(as.data.frame(no_code %>% count(Site, Habitat, name = "records")), row.names = FALSE)
  }
  
  # --- Warn: field code on the wrong transect ---
  if ("Site" %in% names(joined)) {
    site_tab <- joined %>% st_drop_geometry() %>%
      mutate(Site_clean = trimws(Site)) %>%
      filter(!is.na(Field.Code), Field.Code %in% unlist(site_codes))
    expected_site <- setNames(rep(names(site_codes), lengths(site_codes)), unlist(site_codes))
    mismatch <- site_tab %>% filter(Site_clean != expected_site[Field.Code])
    if (nrow(mismatch) > 0) {
      warning(nrow(mismatch), " record(s) in '", input_layer,
              "' have a field code that belongs to the other transect:\n",
              paste(sprintf("  %s (%s, field %s)", mismatch$record_id, mismatch$Site_clean, mismatch$Field.Code),
                    collapse = "\n"), call. = FALSE)
    }
  }
  
  # --- Summary ---
  message("Records per field code:")
  print(table(joined$Field.Code, useNA = "ifany"))
  
  results[[label]] <- list(data = joined %>% select(-poly_id), n_in = nrow(points),
                           output_layer = output_layer, csv_path = csv_path)
}

## -----------------------------------------------------------------
## 3. All layers passed the checks: write the outputs
## -----------------------------------------------------------------
for (res in results) {
  st_write(res$data, gpkg_path, layer = res$output_layer, append = FALSE, quiet = TRUE)
  write.csv(st_drop_geometry(res$data), res$csv_path, row.names = FALSE)
  message("Written: layer '", res$output_layer, "' (", nrow(res$data), " records, input had ",
          res$n_in, ") and ", basename(res$csv_path))
}

message("\nDone.")
