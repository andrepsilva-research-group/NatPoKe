## checkup_IUCNrangeAreas.R ##
## Inês Silva ##
## 21 Nov 2026 ##
## Goal? Check whether species IUCN ranges intersect with both tropical and boreal areas

source(file.path("src", "libraries.R"))
source(file.path("src", "customFunctions2.R"))


##########
# STEP 1 # LOAD TARGET SPECIES
##########

target_species_df <- read.csv("./data/species_by_region.csv", stringsAsFactors = FALSE)
target_species_df <- target_species_df %>%
  mutate(sci_name = gsub(" ", ".", sci_name))

# Extract scientific names.
# Adjust sci_name if your CSV uses a different column name.
target_species <- target_species_df %>%
  mutate(sci_name = gsub(" ", ".", sci_name)) %>%
  distinct(sci_name) %>%
  pull(sci_name)


##########
# STEP 2 # LOAD IUCN RANGES
##########

iucn <- vect("./data/externaldata/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
iucn$sci_name <- gsub(" ", ".", iucn$sci_name)
# Keep only target species
iucn <- iucn[iucn$sci_name %in% target_species, ]
# Remove empty geometries
iucn <- iucn[!is.empty(iucn), ]
# Keep only extant distribution areas
iucn <- iucn[iucn$presence %in% c(1,3),]
invisible(gc())


##########
# STEP 3 # LOAD BIOME REGIONS   #
##########

biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Boreal Forests/Taiga")

# Load biome polygons using your existing function
biome_list <- lapply(biomes, function(target_biome) {
  
  extent_sf <- load_biome(target_biome)
  
  # Convert to terra vector if necessary
  if (inherits(extent_sf, "sf")) {
    extent_sf <- vect(extent_sf)
  }
  # Ensure same CRS as IUCN ranges
  extent_sf <- project(extent_sf, crs(iucn))
  
  extent_sf
})

names(biome_list) <- c("Tropical", "Boreal")


# Robinson projection
robinson_crs <- "ESRI:54030"

iucn_area <- project(iucn, robinson_crs)
# Convert IUCN ranges to sf
iucn_sf <- st_as_sf(iucn_area)

# Repair invalid geometries
iucn_sf <- st_make_valid(iucn_sf)

# Repair biome geometries too
biome_sf <- lapply(biome_list, function(x) {
  
  x <- st_as_sf(x)
  x <- st_make_valid(x)
  x <- st_collection_extract(x, "POLYGON", warn = FALSE)
  x <- x[!st_is_empty(x), ]
  x
})

names(biome_sf) <- c("Tropical", "Boreal")

# Project to equal-area CRS for calculations
#equal_area_crs <- "ESRI:54009"

#iucn_area <- st_transform(iucn_sf, equal_area_crs)

biome_area <- lapply(biome_sf, st_transform, crs = robinson_crs)

################################
# STEP 5: CALCULATE INTERSECTIONS
################################

intersection_results <- list()

for (species in target_species) {
  
  message("Processing: ", species)
  
  tryCatch({
    
    # Select species range
    species_range <- iucn_sf %>%
      filter(sci_name == species)
    
    if (nrow(species_range) == 0) next
    
    # Combine species polygons
    species_range <- st_union(st_geometry(species_range))
    
    species_range <- st_make_valid(species_range)
    
    # Calculate total range area
    total_area <- as.numeric(sum(st_area(species_range)) / 1e6)
    
    for (biome_name in names(biome_area)) {
      
      biome_polygon <- st_union(st_geometry(biome_area[[biome_name]]))
      
      # Default intersection area
      intersection_area <- 0
      
      # Calculate intersection
      intersection <- suppressWarnings(st_intersection(species_range, biome_polygon))
      
      # Calculate area if intersection exists
      if (length(intersection) > 0 &&
          !all(st_is_empty(intersection))) {
        
        intersection_area <- as.numeric(sum(st_area(intersection)) / 1e6)
      }
      
      # Store results
      intersection_results[[paste(species, biome_name, sep = "_")]] <- data.frame(sci_name = species,
                                                                                  biome = biome_name,
                                                                                  total_range_km2 = total_area,
                                                                                  intersection_km2 = intersection_area,
                                                                                  percentage_range = ifelse(total_area > 0, 100 * intersection_area / total_area,
                                                                                                            NA_real_),
                                                                                  intersects = intersection_area > 0)
    }
    
  }, error = function(e) {
    
    message("ERROR processing ", species, ": ", conditionMessage(e))
    
  })
}



##########
# STEP 6 # BUILD SUMMARY TABLE  
##########

intersection_df <- bind_rows(intersection_results)

# Convert to wide format: one row per species
intersection_summary <- intersection_df %>%
  dplyr::select(sci_name, biome, intersection_km2, percentage_range, intersects) %>%
  pivot_wider(names_from = biome,
              values_from = c(intersection_km2, percentage_range, intersects),
              values_fill = list(intersection_km2 = 0, percentage_range = 0, intersects = FALSE))

# Add total range area
range_areas <- intersection_df %>%
  dplyr::select(sci_name, total_range_km2) %>%
  distinct()

intersection_summary <- intersection_summary %>%
  left_join(range_areas, by = "sci_name") %>%
  #left_join(target_species_df, by = "sci_name") %>% 
  distinct() %>% 
  mutate(
    intersects_both = intersects_Tropical & intersects_Boreal,
    intersection_category = case_when(
      intersects_both ~ "Both",
      intersects_Tropical ~ "Tropical only",
      intersects_Boreal ~ "Boreal only",
      TRUE ~ "Neither"
    )
  ) %>%
  select(
    sci_name,
    total_range_km2,
    intersection_km2_Tropical,
    percentage_range_Tropical,
    intersects_Tropical,
    intersection_km2_Boreal,
    percentage_range_Boreal,
    intersects_Boreal,
    intersects_both,
    intersection_category,
    #BIOME_NAME, CONTINENT
  ) %>%
  arrange(desc(intersects_both), sci_name)
#%>% 
 # rename(biomeModelled = BIOME_NAME,
  #       continentModelled = CONTINENT)


################################
# STEP 7: SAVE SUMMARY TABLE   #
################################

dir.create("E:/metaRange_May26/outputs/checkup_reachableAreas/IUCN_biome_intersections",
           recursive = TRUE,showWarnings = FALSE)

write.csv(intersection_summary,
  "E:/metaRange_May26/outputs/checkup_reachableAreas/IUCN_biome_intersections/IUCN_biome_intersection_summary.csv",
  row.names = FALSE)

# Print species intersecting both biomes
both_species <- intersection_summary %>%
  filter(intersects_both)

print(both_species)


################################
# STEP 8: SAVE INTERSECTION MAPS
################################

plot_dir <- "E:/metaRange_May26/outputs/checkup_reachableAreas/IUCN_biome_intersections"

dir.create(plot_dir, recursive = TRUE, showWarnings = TRUE)


library(rnaturalearth)
# Download whole-world land polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Project to Robinson
world <- st_transform(world, crs = robinson_crs)
# Full biome polygons
tropical_biome <- biome_area$Tropical
boreal_biome   <- biome_area$Boreal

# Bounding boxes (as sf objects)
tropical_bbox <- st_as_sf(st_as_sfc(st_bbox(tropical_biome))) %>%
  mutate(biome = "Tropical biome extent")

boreal_bbox <- st_as_sf(st_as_sfc(st_bbox(boreal_biome))) %>%
  mutate(biome = "Boreal biome extent")

# Convert spatial data to sf for plotting
biome_sf <- lapply(biome_area, st_as_sf)

for (species in both_species$sci_name) {
  
  message("Creating map: ", species)
  
  species_range <- iucn_sf[iucn_area$sci_name == species,]
  
  #species_range <- aggregate(species_range)
  
  # Convert range to sf
  #range_sf <- st_as_sf(species_range)
  
  # Intersect with each biome
  tropical_intersection <- st_intersection(species_range, biome_area$Tropical)
  
  boreal_intersection <- st_intersection(species_range, biome_area$Boreal)
  
  p <- ggplot() +
    
    # Whole world background (excluded from legend)
    geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.15) +
    
    # Full biome polygons
    geom_sf(data = tropical_biome, aes(fill = "Full tropical biome"), color = NA, alpha = 0.35) +
    
    geom_sf(data = boreal_biome, aes(fill = "Full boreal biome"), color = NA, alpha = 0.35) +
    
    # Full IUCN range
    geom_sf(data = species_range, aes(fill = "Extant IUCN range"), color = "grey35", linewidth = 0.25) +
    
    # Tropical intersection
    geom_sf(data = tropical_intersection, aes(fill = "Tropical intersection"), color = "#557f30", alpha = 0.85, linewidth = 0.3) +
    
    # Boreal intersection
    geom_sf(data = boreal_intersection, aes(fill = "Boreal intersection"), color = "#016e84", alpha = 0.85, linewidth = 0.3) +
    
    # Biome bounding boxes
    geom_sf(data = tropical_bbox, aes(color = "Tropical biome extent"),
            fill = NA, linewidth = 0.6, linetype = "dashed") +
    
    geom_sf(data = boreal_bbox, aes(color = "Boreal biome extent"),
            fill = NA, linewidth = 0.6, linetype = "dashed") +
    
    # Legend colours
    scale_fill_manual(name = "Range category",
      values = c("Extant IUCN range" = "grey75",
                 "Full tropical biome" = "#92c956",
                 "Full boreal biome" = "#00aaad",
                 "Tropical intersection" = "#557f30",
                 "Boreal intersection" = "#016e84"),
      breaks = c("Extant IUCN range",  "Full tropical biome",
                 "Full boreal biome", "Tropical intersection", "Boreal intersection")) +
    # Bounding box legend
    scale_color_manual(name = "Biome bounding boxes",
      values = c("Tropical biome extent" = "#557f30", "Boreal biome extent" = "#016e84")) +
    labs(title = gsub("\\.", " ", species),
         subtitle = "IUCN range intersections with tropical and boreal biomes") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "italic", size = 14),
      plot.subtitle = element_text(size = 10),
      panel.grid = element_line(color = "grey90", linewidth = 0.2),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9)
    ) +
    coord_sf(crs = robinson_crs, datum = NA)
  
  # Save individual plot
  ggsave(filename = file.path(plot_dir, paste0(species, "_biome_intersections.png")),
    plot = p, bg ="white", width = 10, height = 7, dpi = 300)
}
