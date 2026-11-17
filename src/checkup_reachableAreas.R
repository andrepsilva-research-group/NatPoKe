## checkup_reachableAreas.R ##
## Inês Silva ##
## 17 Nov 2026 ##
## GOAL? Create a piece of code to maks current suitability maps to include only
## reachable areas for each species, based on the patches() function to find
## different clusters of land and intersecting them with the current species
## IUCN range.
## **Reachable areas** are defined here as all areas currently connected with
## the species range. This effectively removes islands and other territories
## separated by ocean/waterways.

##############################
# STEP-BY-STEP DEMONSTRATION #
##############################

# test species = Felis Chaus

# load pre-made suitability map to test
land_rast <- rast("E:/metaRange_May26/outputs/checkup_fullModellingArea/Asia_ssp126_20261110/Inputs/Felis.chaus_tropical_ssp126_cropped_reprojectedm.tif")
land_rast <- r_utm

# load iucn's species ranges 
iucn <- vect("./data/externaldata/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
iucn$sci_name <- gsub(" ", ".", iucn$sci_name)
invisible(gc())
# felix chaus
this_iucn_range <- iucn[iucn$sci_name == "Felis.chaus",]

# check that both the land_rast and the iucr ange use the same CRS
if (!terra::same.crs(land_rast, this_iucn_range)) {
  this_iucn_range <- terra::project(this_iucn_range, terra::crs(land_rast))
}


##########
# STEP 1 # Creat a binary land raster
##########

land_binary <- terra::ifel(!is.na(land_rast) & land_rast > 0, 1, NA)
plot(land_binary$`2015`)

##########
# STEP 2 # Identify connected land patches
##########

land_patches <- terra::patches(land_binary, directions = 8) 
plot(land_patches$`2015`, main = "Reachable Land for Felis Chaus, 8 directions")

##########
# STEP 3 # Identify reachable areas (land patches that intersect with the IUCN range)
##########

range_rast <- terra::rasterize(this_iucn_range, land_patches,
                               field = 1, background = NA, touches = TRUE)

range_patches <- terra::mask(land_patches, range_rast)

# these are the patches IDs
accessible_patch_ids <- unique(terra::values(range_patches, mat = FALSE))
accessible_patch_ids <- accessible_patch_ids[!is.na(accessible_patch_ids)]
accessible_patch_ids

# Stop if the IUCN range does not overlap any land cells
if (length(accessible_patch_ids) == 0) {
  stop("IUCN range does not intersect the land raster.")
}

# keep only land patches connected to the current range
accessibility_mask <- terra::subst(land_patches,
                                   from = accessible_patch_ids,
                                   to = rep(1, length(accessible_patch_ids)),
                                   others = 0)

# convert the ocean and other NA cells to zero
accessibility_mask <- terra::ifel(is.na(land_patches), 0, accessibility_mask)
plot(accessibility_mask$`2015`)


# Step-by-step process converted to function

# Modelling area with reachable land (land area connected to IUCN range)

make_reachableAreas_mask <- function(land_rast, iucn_range, directions = 4) {
  
  # Check that both inputs use the same CRS
  if (!terra::same.crs(land_rast, iucn_range)) {
    iucn_range <- terra::project(iucn_range, terra::crs(land_rast))
  }
  
  # Step 1 # Create a binary land raster
  ## Land = 1; water = NA
  land_binary <- terra::ifel(!is.na(land_rast[[1]]) & land_rast[[1]] > 0, 1, NA)
  
  # Step 2 # Identify connected land patches
  land_patches <- terra::patches(land_binary,
                                 directions = 4)
  
  # Step 3 # Identify reachable areas
  # rasterise the current IUCN range
  range_rast <- terra::rasterize(iucn_range,
                                 land_patches,
                                 field = 1, background = NA, touches = TRUE)
  
  # identify land patches intersecting the IUCN range
  range_patches <- terra::mask(land_patches, range_rast)
  
  accessible_patch_ids <- unique(terra::values(range_patches, mat = FALSE))
  accessible_patch_ids <- accessible_patch_ids[!is.na(accessible_patch_ids)]
  
  # Stop if the IUCN range does not overlap any land cells
  if (length(accessible_patch_ids) == 0) {
    stop("IUCN range does not intersect the land raster.")
  }
  
  # keep only land patches connected to the current range
  accessibility_mask <- terra::subst(land_patches,
                                     from = accessible_patch_ids,
                                     to = rep(1, length(accessible_patch_ids)),
                                     others = 0)
  
  # convert ocean and other NA cells to zero
  accessibility_mask <- terra::ifel(is.na(land_patches), 0, accessibility_mask)
  
  names(accessibility_mask) <- "accessibility"
  
  return(accessibility_mask)
}



#####################################################
# WHAT'S THE DIFFERENCE BETWEEN 4 AND 8 DIRECTIONS? #
#####################################################

land_binary <- terra::ifel(!is.na(land_rast) & land_rast > 0, 1, NA)

# Select 2015 layer
land_2015 <- land_binary[["2015"]]


run_connectivity <- function(land_binary, iucn_range, directions) {
  
  # Identify connected land patches
  land_patches <- terra::patches(
    land_binary,
    directions = directions
  )
  
  # Rasterize IUCN range
  range_rast <- terra::rasterize(
    iucn_range,
    land_patches,
    field = 1,
    background = NA,
    touches = TRUE
  )
  
  # Identify patches intersecting the IUCN range
  range_patches <- terra::mask(
    land_patches,
    range_rast
  )
  
  accessible_patch_ids <- unique(
    terra::values(range_patches, mat = FALSE)
  )
  
  accessible_patch_ids <- accessible_patch_ids[
    !is.na(accessible_patch_ids)
  ]
  
  if (length(accessible_patch_ids) == 0) {
    stop("IUCN range does not intersect any land patches.")
  }
  
  # Keep patches connected to the IUCN range
  accessibility_mask <- terra::subst(
    land_patches,
    from = accessible_patch_ids,
    to = rep(1, length(accessible_patch_ids)),
    others = 0
  )
  
  # Convert ocean cells to zero
  accessibility_mask <- terra::ifel(
    is.na(land_patches),
    0,
    accessibility_mask
  )
  
  return(list(
    patches = land_patches,
    accessibility = accessibility_mask,
    accessible_ids = accessible_patch_ids
  ))
}


# --------------------------------------------------
# STEP 3: Run both connectivity options
# --------------------------------------------------

results_4 <- run_connectivity(
  land_2015,
  this_iucn_range,
  directions = 4
)

results_8 <- run_connectivity(
  land_2015,
  this_iucn_range,
  directions = 8
)


# --------------------------------------------------
# STEP 4: Create figures
# --------------------------------------------------

plot_connectivity <- function(results, directions, land_binary) {
  
  # Number of connected patches
  n_patches <- terra::global(
    results$patches,
    "max",
    na.rm = TRUE
  )[1, 1]
  
  n_accessible <- length(results$accessible_ids)
  
  # Create patch colour palette
  patch_colours <- hcl.colors(
    max(1, as.integer(n_patches)),
    palette = "Dynamic"
  )
  
  # Save figure
  filename <- paste0(
    "ReachableAreas_1km_", directions, "_directions.png"
  )
  
  png(
    filename,
    width = 3600,
    height = 1400,
    res = 300
  )
  
  par(
    mfrow = c(1, 3),
    mar = c(3, 3, 5, 5),
    oma = c(0, 0, 3, 0)
  )
  
  # Panel A: Binary land
  terra::plot(
    land_binary,
    main = "A. Binary land raster",
    col = "grey40",
    axes = FALSE,
    legend = FALSE
  )
  
  # Panel B: Connected patches
  terra::plot(
    results$patches,
    main = paste0(
      "B. Connected land patches\n",
      "Total patches: ",
      format(n_patches, big.mark = ",")
    ),
    col = patch_colours,
    axes = FALSE,
    legend = FALSE
  )
  
  # Panel C: Accessibility mask
  terra::plot(
    results$accessibility,
    main = paste0(
      "C. Accessible land patches\n",
      "Connected patches: ",
      n_accessible
    ),
    col = c("grey90", "#238B8E"),
    axes = FALSE,
    legend = FALSE
  )
  
  mtext(
    paste0(
      "Reachable Areas for Felis Chaus using 1km and",
      directions,
      "-direction neighbourhood"
    ),
    outer = TRUE,
    cex = 1.6,
    font = 2
  )
  
  dev.off()
}


# --------------------------------------------------
# STEP 5: Export both figures
# --------------------------------------------------

plot_connectivity(results_4, 4, land_2015)
plot_connectivity(results_8, 8, land_2015)



comparison <- data.frame(
  Directions = c(4, 8),
  Total_patches = c(
    global(results_4$patches, "max", na.rm = TRUE)[1, 1],
    global(results_8$patches, "max", na.rm = TRUE)[1, 1]
  ),
  Accessible_patches = c(
    length(results_4$accessible_ids),
    length(results_8$accessible_ids)
  ),
  Accessible_land_cells = c(
    global(results_4$accessibility, "sum", na.rm = TRUE)[1, 1],
    global(results_8$accessibility, "sum", na.rm = TRUE)[1, 1]
  )
)

print(comparison)
