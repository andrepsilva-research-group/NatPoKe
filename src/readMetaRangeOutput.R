## Name: readMetaRangeOutput.R ##
## Author: Inês Silva ##
## Date: November 25th 2025 ##
## Goal? Read metaRange ouputs and perform some basic diagnostics including:
    ##
## .csv file with complete metaRange runs for all regions and scenarios alongsode metadata file.
## diagnostic population trends using the TNIND over time, a trait dataframe for all
## modeled species for all regions and average suitability over time plots

source("./src/libraries.R")
source("./src/customFunctions2.R")  


# -------------------------------------------------------------------------
# CONFIGURATION: edit paths, run date and analysis selections here only
# -------------------------------------------------------------------------

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

src_dir       <- file.path(project_dir, "src")
data_dir      <- file.path(project_dir, "data")
external_dir  <- file.path(data_dir, "externaldata")

output_root   <- "E:/metaRange_May26/outputs/checkup_reachableAreas"
run_date      <- "20261118"

run_table_file <- file.path(data_dir, "run_table_longBurnIn.csv")
traits_file    <- file.path(external_dir, "mammalTraits_2025-12-11.csv")
iucn_file      <- file.path(external_dir, "MAMMALS_TERRESTRIAL_ONLY",
                            "MAMMALS_TERRESTRIAL_ONLY.shp")

diagnostics_dir <- file.path(output_root, "diagnostics")
SDMsFigures_dir <- file.path(output_root, "SDMsFigures")

# Create output directories once, before any section writes to them.
dir.create(output_root, recursive = TRUE, showWarnings = TRUE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = TRUE)
dir.create(SDMsFigures_dir, recursive = TRUE, showWarnings = TRUE)


# Read the run table once and construct all run-specific directories once.
runs <- read.csv(run_table_file, stringsAsFactors = FALSE)

run_paths <- runs %>%
  mutate(
    runname = paste(region, scenario, run_date, sep = "_"),
    input_dir = file.path(output_root, runname, "Inputs"),
    output_dir = file.path(output_root, runname, "Outputs")
  ) %>%
  distinct(runname, .keep_all = TRUE)

selected_sps <- c(
  # tropical Africa
  "Loxodonta africana",
  "Caracal caracal",
  "Mandrillus sphinx",
  
  # tropical South America
  "Odocoileus virginianus",
  "Ozotoceros bezoarticus",
  
  # boreal North America
  "Bison bison",
  "Canis latrans",
  "Gulo gulo",
  
  # boreal Europe + Asia
  "Panthera tigris", 
  "Bison bonasus", 
  "Capreolus pygargus",
  "Alces alces",
  
  # tropical asia
  "Axis axis",
  "Felis chaus",
  "Canis lupus"
)



# Load Natural Earth countries
world_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

# Function to retrieve the geographical region from the run name
get_region_map <- function(runname, raster_crs) {
  
  # Extract region from run name
  region_name <- strsplit(runname, "_")[[1]][1]
  
  # Match run names to Natural Earth continent names
  continents <- switch(
    region_name,
    "Africa" = "Africa",
    "Asia" = "Asia",
    "SouthAmerica" = "South America",
    "South America" = "South America",
    "NorthAmerica" = "North America",
    "North America" = "North America",
    "Europe+Asia" = c("Europe", "Asia"),
    NULL
  )
  
  if (is.null(continents)) {
    warning("Unknown region in runname: ", runname)
    return(NULL)
  }
  
  # Select countries belonging to the region
  region_sf <- world_sf %>%
    filter(continent %in% continents)
  
  # Dissolve countries into one regional polygon
  region_sf <- sf::st_sf(
    region = region_name,
    geometry = sf::st_union(sf::st_geometry(region_sf))
  )
  
  # Transform to raster CRS
  region_sf <- sf::st_transform(
    region_sf,
    crs = raster_crs
  )
  
  return(region_sf)
}


################
# FULL DATASET # 
################

TNIND_paths <- list()

##########
# STEP 1 # collect TNIND paths if not already provided
##########

if (length(TNIND_paths) == 0) {
  
  runs <- read.csv(run_table_file, stringsAsFactors = FALSE)
  
  TNIND_paths <- character(0)
  
  for (i in seq_len(nrow(runs))) {
    
    target_region <- runs$region[i]
    target_biome <- runs$biome[i]
    future_scenario <- runs$scenario[i]
    
    runname <- paste(target_region, future_scenario, run_date,
      #format(Sys.time(), "%Y%m%d"),
      sep = "_")
    
    tnind_file <- file.path(output_root, runname, "Outputs", paste0("TNIND_yr_", runname, ".csv"))
    
    if (!file.exists(tnind_file)) {
      warning("TNIND file not found (skipping): ", tnind_file)
      next
    }
    
    TNIND_paths[runname] <- tnind_file
  }
}

# read all dfs
TNIND_all_runs <- lapply(TNIND_paths, data.table::fread)

# combine into unique df
TNIND_yr <- data.table::rbindlist(TNIND_all_runs, use.names = TRUE, fill = TRUE)

# correction for biome column
TNIND_yr[biome == "TropicalSubtropicalMoistBroadleafForests",
         biome := "Tropical & Subtropical Moist Broadleaf Forests"]

TNIND_yr[biome == "BorealForestsTaiga",
         biome := "Boreal Forests/Taiga"]

# corresction for region name
TNIND_yr[region == "SouthAmerica", region := "South America"]
TNIND_yr[region == "NorthAmerica", region := "North America"]

# clean up
rm(TNIND_all_runs)
invisible(gc())

# get correspondence between species names and functional group
combined_traits_data <- read_csv(traits_file) %>% 
  dplyr::filter(BIOME_NAME %in% c("Tropical & Subtropical Moist Broadleaf Forests", "Boreal Forests/Taiga")) %>% 
  mutate(
    CONTINENT = case_when(
      BIOME_NAME == "Boreal Forests/Taiga" & CONTINENT == "Europe" ~ "Europe+Asia",
      TRUE ~ CONTINENT),
    trophic_level = case_when(
      trophic_level == 1 ~ "Herbivore",
      trophic_level == 2 ~ "Omnivore",
      trophic_level == 3 ~ "Carnivore",
      TRUE ~ as.character(trophic_level))
  )


# clean up the dataset
TNIND_yr <- TNIND_yr %>% 
  mutate(species = pretty_species_names(species)) %>% # if running twice it throws a warning - It's ok!
  left_join(
    dplyr::select(combined_traits_data, sci_name, BIOME_NAME, CONTINENT, trophic_level),
    by = c("species" = "sci_name",
           "biome" = "BIOME_NAME", # keep biome & continent here or a many-to-many warning will appear
           "region" = "CONTINENT"
           )) %>%
  # simplify replicates numbering
  mutate(rep_num = str_extract(rep, "^[0-9]+")) %>% 
  # correction for tigers that are from asia but asian boreal forest are modelled together with europe
  mutate(trophic_level = case_when(
      species == "Panthera tigris" ~ "Carnivore",
      species == "Ursus thibetanus" ~ "Omnivore",
      TRUE ~ trophic_level)
    ) %>%
  # deal with integer 64 columns (=big big numbers)
  mutate(across(where(bit64::is.integer64), as.numeric))

##########
# STEP 3 #  Write complete dataset into .csv (RAW DATA)
##########

write_csv(TNIND_yr, 
          file = file.path(output_root, paste0("completeMetaRangeRun_", #format(Sys.time(), "%Y%m%d"),
                                               run_date,
                                               ".csv")))

# metadata file
metadata <- tibble::tribble(
  ~variable, ~description, ~units, ~data_type,
  
  "TNIND",
  "Total number of individuals across the whole landscape.",
  "Number of individuals",
  "Numeric",
  
  "MNIND",
  "Mean number of individuals per occupied cell.",
  "Individuals per cell",
  "Numeric",
  
  "mean_repRate",
  "Mean reproductive rate of the species across the landscape.",
  "Specify according to model parameterisation",
  "Numeric",
  
  "mean_carrCap",
  "Mean of the maximum population size supported by local habitat conditions",
  "Number of individuals per cell",
  "Numeric",
  
  "occupancy",
  "Number of cells occupied by the species.",
  "Number of cells",
  "Numeric",
  
  "future_scenario",
  "Future climate scenario used in the simulation. ssp126 represents SSP1-2.6 and ssp585 represents SSP5-8.5.",
  "Categorical",
  "Character",
  
  "biome",
  "Biome in which the species population was modelled. Tropical & Subtropical Moist Broadleaf Forests or Boreal Forests/Taiga.",
  "Categorical",
  "Character",
  
  "region",
  "Geographic region in which the species population was modelled. Regions include Asia, Africa, South America, Europe+Asia, and North America",
  "Categorical",
  "Character",
  
  "species",
  "Scientific name of the modelled species.",
  "Text",
  "Character",
  
  "timestep",
  "Sequential time step of the simulation. One timestep corresponds to one simulation year, with 1 to 25 representing the burn-in period, and remaining ones modelled years.",
  "Years / model timesteps",
  "Numeric",
  
  "rep",
  "Identifier of the stochastic simulation replicate.",
  "Text",
  "Character",
  
  "trophic_level",
  "Functional trophic group of the species, Herbivore, Omnivore, or Carnivore.",
  "Categorical",
  "Character",
  
  "rep_num",
  "Numeric identifier of the simulation replicate.",
  "Integer",
  "Integer")

write_csv(metadata, file = file.path(output_root, paste0("completeMetaRangeRun_", run_date, "_metadata.csv")))

#################################
# AVERAGE SUITABILITY OVER TIME #
#################################

# Average-suitability outputs are saved in SDMsFigures_dir (defined above).

inputFolder_paths <- stats::setNames(run_paths$input_dir, run_paths$runname)

for (i in seq_len(nrow(runs))) {
  
  target_region <- runs$region[i]
  target_biome <- runs$biome[i]
  future_scenario <- runs$scenario[i]
  
  runname <- paste(
    target_region,
    future_scenario,
    run_date,
    #format(Sys.time(), "%Y%m%d"),
    sep = "_"
  )
  
  input_folder <- file.path(
    output_root, 
    runname, "Inputs")
  
  
  if (!file.exists(input_folder)) {
    warning("Input folder not found (skipping): ", input_folder)
    next
  }
  
  inputFolder_paths[runname] <- input_folder
}

sps_suitability_ALL <- list()

# loop over folders and produce one plot per folder
for (dir in inputFolder_paths) {
  
  # list tif files in this folder
  files <- list.files(dir, pattern = "_reprojectedm\\.tif$", full.names = TRUE)
  if (length(files) == 0) next
  
  message("Average Suitability Over Time plot for: ", dir)
  
  # process all species inside this folder
  folder_df <- lapply(files, function(f) {
    
    fname <- basename(f)
    
    # get species name 
    species_name <- str_extract(fname, "^(.*?)_(?=(boreal|tropical))")
    species_name <- gsub("_$", "", species_name)  # remove trailing _
    
    # load raster
    r <- rast(f)
    
    # get global mean suitability value per layer
    df <- global(r, "mean", na.rm = TRUE) %>%
      as.data.frame()
    
    # add years & sps names
    df$year <- as.numeric(names(r))
    df$species <- suppressMessages(pretty_species_names(species_name)) #pretty_species_names() is custom function
    df
  }) %>%
    bind_rows()
  
  
  # folder name for saving plot
  folder_name <- basename(dirname(dir))
  
  # create object with everyone's suitbaility values
  sps_suitability_ALL[[folder_name]] <- folder_df
  
  # plot average suitability
  p <- ggplot(folder_df, aes(x = year, y = mean, group = species)) +
    geom_line() +
    facet_wrap(~ species, scales = "free_y", ncol = 3) +
    scale_x_continuous(
      breaks = seq(min(folder_df$year), max(folder_df$year), by = 20)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "italic")) +
    labs(x = "Year", y = "Mean suitability", title = paste("Average suitability –", folder_name))
  
  #print(p)
  n_species <- length(unique(folder_df$species))
  ncol <- 3
  nrow <- ceiling(n_species / ncol)
  
  ggsave(filename = paste0(folder_name, "_suitabilityOverTime.png"),
         path = diagnostics_dir,
         plot = p, width = 12,
         height = nrow * 3,
         dpi = 300)
  invisible(gc())
}


# average suitability per biome, continent and functional group ----------------
# average suitability per biome, continent and functional group ----------------
sps_suitability <- bind_rows(sps_suitability_ALL, .id = "folder") %>%
  mutate(region = sub("_ssp.*", "", folder),
         scenario = sub(".*_(ssp\\d+)_.*", "\\1", folder),
         biome = case_when(region %in% c("Asia", "Africa", "South America") ~ "Tropical & Subtropical Moist Broadleaf Forests",
                           region %in% c("Europe+Asia", "North America") ~ "Boreal Forests/Taiga",
                           TRUE ~ NA_character_)) %>% 
  # add trait df to get trophic levels
  left_join(
    dplyr::select(combined_traits_data, sci_name, BIOME_NAME, CONTINENT, trophic_level),
    by = c("species" = "sci_name",
           "biome" = "BIOME_NAME", # keep biome & continent here or a many-to-many warning will appear
           "region" = "CONTINENT")) %>% 
  # correction for tigers that are from asia but asian boreal forest are modelled together with europe
  mutate(trophic_level = case_when(
    species == "Panthera tigris" ~ "Carnivore",
    species == "Ursus thibetanus" ~ "Omnivore",
    TRUE ~ trophic_level))


trophic_suit <- sps_suitability %>% 
  group_by(scenario, biome, region, trophic_level, year) %>% 
  # mean suitability per trophic group
  dplyr::summarise(avg_suit = round(mean(mean, na.rm = TRUE), 3)) %>% 
  dplyr::filter(year %in% c(2015, 2100)) %>% 
  # Put years into separate columns
  pivot_wider(names_from = year, values_from = avg_suit, names_prefix = "suit_") %>% 
  # Calculate absolute and percentage change
  mutate(perc_change = (suit_2100 - suit_2015) * 100) %>% 
  arrange(scenario, biome, region, trophic_level) %>% 
  rename("Scenario" = "scenario", "Biome" = "biome", "Region" = "region", 
         "Trophic Level" = "trophic_level", "Mean Suitability in 2015" = "suit_2015",
         "Mean Suitability in 2100" = "suit_2100", "Change between 2100 and 2015 (%)" = "perc_change")

# save average suitability values per trophic group as csv
write.csv(trophic_suit,
          file = file.path(SDMsFigures_dir,
                           paste0("AverageSuitabilityPerTrophicGroup", Sys.Date(), ".csv")), row.names = FALSE)


# average suitability per biome ------------------------------------------------

biome_suit <- sps_suitability %>% 
  group_by(scenario, biome, #region,
           year) %>% 
  # mean suitability per trophic group
  dplyr::summarise(avg_suit = round(mean(mean, na.rm = TRUE), 3)) %>% 
  dplyr::filter(year %in% c(2015, 2100)) %>% 
  # Put years into separate columns
  pivot_wider(names_from = year, values_from = avg_suit, names_prefix = "suit_") %>% 
  # Calculate absolute and percentage change
  mutate(perc_change = (suit_2100 - suit_2015) * 100) %>% 
  arrange(scenario, biome#, region
  ) %>% 
  rename("Scenario" = "scenario", "Biome" = "biome", #"Region" = "region", 
         "Mean Suitability in 2015" = "suit_2015", "Mean Suitability in 2100" = "suit_2100",
         "Change between 2100 and 2015 (%)" = "perc_change")

# save average suitability values per biome csv
write.csv(biome_suit,
          file = file.path(SDMsFigures_dir,
                           paste0("AverageSuitabilityPerBiome", Sys.Date(), ".csv")), row.names = FALSE)


##########################################
# DIAGNOSTIC POPULATION TRENDS OVER TIME # 
##########################################

##########
# STEP 1 # Build and excel file
##########

# Diagnostic outputs are saved in diagnostics_dir (defined above).

# get number of sps per combin
n_sps <- TNIND_yr %>% 
  group_by(future_scenario, biome, region) %>% 
  summarise(n_species = n_distinct(species), .groups = "drop")

# population trends ------------------------------------------------------------
## differences within EACH replicate
TNIND_diff <- TNIND_yr %>%
  arrange(species, biome, rep_num, timestep) %>%
  group_by(future_scenario, biome, region, species, rep_num) %>%
  mutate(diff_TNIND = TNIND - lag(TNIND)) %>%
  ungroup() %>% 
  dplyr::select(future_scenario, biome, region, species, timestep, rep_num, TNIND, diff_TNIND) %>% 
  arrange(future_scenario, biome, region, species)

## mean TNIND and diff ACROSS replicates
TNIND_mean <- TNIND_diff %>%
  group_by(future_scenario, biome, region, species, timestep) %>%
  summarise(mean_TNIND = round(mean(TNIND, na.rm = TRUE), 0)) %>% 
  mutate(mean_diff_TNIND = mean_TNIND - lag(mean_TNIND, n = 1),
         rep_num = "Mean") %>%
  ungroup() %>% 
  dplyr::select(future_scenario, biome, region, species, timestep, mean_TNIND, mean_diff_TNIND, rep_num)

# save diagnostics in .xlsx
# create workbook
wb <- createWorkbook()
addWorksheet(wb, "spsModelled")
writeData(wb, "spsModelled", n_sps)
addWorksheet(wb, "perReplicate")
writeData(wb, "perReplicate", TNIND_diff)
# add across-replicate sheet 
addWorksheet(wb, "acrossReplicates")
writeData(wb, "acrossReplicates", TNIND_mean)
saveWorkbook(wb, file.path(diagnostics_dir, paste0("metaRangeRun_", 
                                                   run_date,
                                                   "_diagnostics.xlsx")), overwrite = TRUE)
invisible(gc())

##########
# STEP 2 # Plot population trends over time 
##########

### per biome × region × scenario (all species together)
combo_list <- TNIND_diff %>%
  distinct(biome, region, future_scenario)

for (i in seq_len(nrow(combo_list))) {
  
  b <- combo_list$biome[i]
  r <- combo_list$region[i]
  s <- combo_list$future_scenario[i]
  
  # filter data for this combination
  df <- TNIND_mean %>%
    filter(biome == b,
           region == r,
           future_scenario == s)
  
  # skip if nothing there
  if (nrow(df) == 0) next
  
  # nice title
  biome_title <- gsub("([A-Z])", " \\1", b) |> trimws()
  
  # plot (all species together)
  p <- ggplot(df, aes(x = timestep, y = mean_TNIND)) +
    geom_line() +
    # add line for burn-in
    geom_vline(xintercept = 25, linetype="dotted", 
               color = "blue") +
    # add line for new burnin
    geom_vline(xintercept = 100, linetype="dotted", 
               color = "red") +
    # add yet another line for new burnin
    geom_vline(xintercept = 300, linetype="dotted", 
               color = "green") +
    facet_wrap(~species, scales = "free_y") +
    labs(title = paste0(biome_title, " – ", r, " – ", s),
         subtitle = "Species population trends") +
    theme_minimal() +
    theme(strip.text = element_text(face = "italic"))
  
  # safe biome name for files
  b <- gsub("[/& ]", "", b)
  
  # filename
  fname <- file.path(diagnostics_dir,
                     paste0(gsub(" ", "", b), "_", gsub(" ", "", r), "_", gsub(" ", "", s),
                            "_speciesPopulationTrends.png"))
  # save plot
  ggsave(filename = fname, plot = p,
         bg = "white", width = 350, height = 210, units = "mm", dpi = 300)
  rm(df, p)
  gc()
}


### per species ----------------------------------------------------------------
# get unique combination to plot (biome + region + species + scenario)
combo_list <- TNIND_diff %>%
  distinct(biome, region, species, future_scenario)

# go through each combination
for (i in seq_len(nrow(combo_list))) {
  biome_to_plot   <- combo_list$biome[i]
  region_to_plot  <- combo_list$region[i]
  sp              <- combo_list$species[i]
  scen_to_plot    <- combo_list$future_scenario[i]
  
  # filter for that combo
  sp_mean <- TNIND_mean %>%
    filter(biome == biome_to_plot,
           region == region_to_plot,
           species == sp,
           future_scenario == scen_to_plot,
           timestep > 25)       # remove burn in
  
  sp_data <- TNIND_diff %>%
    filter(biome == biome_to_plot,
           region == region_to_plot,
           species == sp,
           future_scenario == scen_to_plot)
  
  # skip if data missing
  if (nrow(sp_data) == 0 | nrow(sp_mean) == 0) next
  
  # right plot - mean TNIND across replicates
  p_right <- ggplot(sp_mean, aes(x = timestep, y = mean_TNIND)) +
    geom_line(color = "black", size = 0.8) +
    labs(
      title = paste(sp, "- Mean across replicates"),
      x = "Timestep", y = "Mean TNIND") +
    theme_minimal()
  
  # left plot - TNIND per replicate separately
  p_left <- ggplot(sp_data, aes(
    x = timestep, y = TNIND,
    group = factor(rep_num),
    color = factor(rep_num))) +
    geom_line(size = 0.8, alpha = 0.7) +
    geom_vline(xintercept = 25, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 100, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 300, linetype = "dashed", color = "green") +
    labs(
      title = paste(sp, "- Replicates"),
      x = "Timestep", y = "TNIND", color = "Replicate") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # combine both
  combined_plot <- p_left + p_right + 
    plot_annotation(title = paste("Temporal dynamics of TNIND", "\nBiome:", biome_to_plot,
                                  "| Region:", region_to_plot, "| Scenario:", scen_to_plot))
  
  # safe biome name for files
  biome_to_plot <- gsub("[/& ]", "", biome_to_plot)
  
  # filename 
  filename <- paste0("TNIND_", gsub(" ", "_", sp), "_", gsub(" ", "_", biome_to_plot), "_",
                     gsub(" ", "_", region_to_plot), "_", gsub(" ", "_", scen_to_plot), ".png")
  
  # save plots
  ggsave(filename, combined_plot, path = diagnostics_dir,
         width = 16, height = 6, dpi = 300)
  
  message("Saved: ", filename)
  
  gc(rm(sp_data, sp_mean, p_left, p_right, combined_plot,
        biome_to_plot, region_to_plot, scen_to_plot, sp))
}

############################
# COMPLETE TRAIT DATAFRAME #
############################

traitdf_paths <- file.path(run_paths$input_dir, "metaRangeSpeciesDataframe.csv")
names(traitdf_paths) <- run_paths$runname
traitdf_paths <- traitdf_paths[file.exists(traitdf_paths)]

# read all dfs
traits_all_runs <- lapply(traitdf_paths, data.table::fread)

# combine into unique df
all_traits <- data.table::rbindlist(traits_all_runs, use.names = TRUE, fill = TRUE) %>% 
  select(!Index) %>%   # remove Index column if present 
  distinct(Species, .keep_all = TRUE) %>% # ensure one row per species
  mutate(Species = pretty_species_names(as.character(Species)))

# clean up
rm(traits_all_runs)
invisible(gc())

# Save merged .xlsx
write_xlsx(all_traits, file.path(output_root, "completeTraitDataframe_allSps.xlsx"))

###################################################
# DIAGNOSTIC POPULATION TRENDS SPATIALLY EXPLICIT # 
###################################################

# input/output folders
inputFolder_paths <- stats::setNames(run_paths$input_dir, run_paths$runname)
outputFolder_paths <- stats::setNames(run_paths$output_dir, run_paths$runname)

# Keep only runs whose input directory exists.
valid_runs <- dir.exists(inputFolder_paths)
inputFolder_paths <- inputFolder_paths[valid_runs]
outputFolder_paths <- outputFolder_paths[valid_runs]

# load IUCN ranges
iucn <- vect(iucn_file)
iucn$sci_name <- gsub(" ", ".", iucn$sci_name)
invisible(gc())

# chose timesteps to check
timestep_1 <- 1
timestep_2 <- 2
#timestep_100 <- 100
timestep_300 <- 300

# go through folders to find correct rasters
for (dir in inputFolder_paths) {
  
  # Identify run
  runname <- basename(dirname(dir))
  
  message("")
  message("==============================================")
  message("RUN: ", runname)
  message("==============================================")
  
  
  ##########
  # STEP 2 # Find suitability rasters
  ########## 
  
  files <- list.files(dir, pattern = "_reprojectedm\\.tif$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No suitability rasters found in: ", dir)
    next}
  
  # Keep ONLY selected species
  files_selected <- files[
    sapply(files, function(f) {
      
      fname <- basename(f)
      species_name <- str_extract(fname, "^(.*?)_(?=(boreal|tropical))")
      species_name <- gsub("_$", "", species_name)
      species_name %in% gsub(" ", ".", selected_sps)
    })
  ]
  
  if (length(files_selected) == 0) {
    message("No selected species found in: ", dir)
    next
  }
  
  for (f in files_selected) {
    
    # Extract species name
    fname <- basename(f)
    
    species_name <- str_extract(fname, "^(.*?)_(?=(boreal|tropical))")
    
    species_name <- gsub("_$", "", species_name)
    species_label <- suppressMessages(pretty_species_names(species_name))
    
    message("")
    message("Processing: ", species_label, " | ", runname)
    # Load suitability rasters
    r <- rast(f)
    
    message("Suitability layers: ", nlyr(r))
    
    if (max(timestep_1, timestep_2) > nlyr(r)) {
      
      warning(species_label,
              ": requested timestep exceeds number of suitability raster layers.")
      next
    }
    
    # Extract correct suitability timesteps
    r1 <- r[[timestep_1]]
    r2 <- r[[timestep_1]]
    r300 <- r[[timestep_1]] # these need to be layer numbers and not the timesteps!!! timestep 100 is still layer 1 
    
    
    names(r1) <- paste0("Timestep ", timestep_1)
    names(r2) <- paste0("Timestep ", timestep_2)
    names(r300) <- paste0("Timestep ", timestep_300)
    
    ############
    # STEP 2.1 # Crop first layer by the IUCN range (just like the burnin in metaRange)
    ############
    
    species_range <- iucn[iucn$sci_name == species_name, ]
    region <- str_split(runname, "_")[[1]]
    continent_sf <- load_select_continents(region[1])
    
    if (nrow(species_range) == 0) {
      warning("No IUCN range found for: ", species_name)
      next
    }
    
    # Project geometries to suitability raster CRS
    species_range <- terra::project(species_range, terra::crs(r1))
    continent_sf <- sf::st_transform(continent_sf, terra::crs(r1))
    
    # 1. Crop and mask by continent
    r1 <- terra::crop(r1, terra::vect(continent_sf))
    r1 <- terra::mask(r1, terra::vect(continent_sf))
    
    # 2. Crop and mask by species range
    r1 <- terra::crop(r1, species_range)
    r1 <- terra::mask(r1, species_range)
    
    
    # STEP 3 # Find outputs folders to retrieve abundance and repRate
    output_dir <- unname(outputFolder_paths[runname])
    
    if (
      length(output_dir) == 0 ||
      is.na(output_dir) ||
      !file.exists(output_dir)
    ) {
      
      warning("Output folder not found for: ", runname,"\nExpected: ", output_dir)
      next
    }
    
    
    message("Output folder: ", output_dir)
    
    ##########
    # STEP 4 # FIND ABUNDANCE RASTERS
    ##########
    
    abundance_files <- list.files(output_dir, pattern = "_abundance\\.tif$",
                                  full.names = TRUE)
    
    if (length(abundance_files) == 0) {
      warning("No abundance rasters found in: ", output_dir)
      next
    }
    # keep only selected timesteps
    abundance_025 <- abundance_files[grepl(paste0("^1_.*_Mammals001_", species_name, "_abundance\\.tif$"), basename(abundance_files))]
    abundance_026 <- abundance_files[grepl(paste0("^1_.*_Mammals002_", species_name, "_abundance\\.tif$"), basename(abundance_files))]
    abundance_300 <- abundance_files[grepl(paste0("^1_.*_Mammals300_", species_name, "_abundance\\.tif$"), basename(abundance_files))]
    
    if (length(abundance_025) == 0 || length(abundance_026) == 0 || length(abundance_300) == 0) {
      warning("Could not find both abundance rasters for ", species_name)
      next
    }
    
    abundance_025 <- abundance_025[1]
    abundance_026 <- abundance_026[1]
    abundance_300 <- abundance_300[1]
    
    message("Abundance timestep 1: ", basename(abundance_025))
    message("Abundance timestep 2: ", basename(abundance_026))
    message("Abundance timestep 300: ", basename(abundance_300))
    
    # load ratsers
    a1 <- rast(abundance_025)
    a2 <- rast(abundance_026)
    a300 <- rast(abundance_300)
    
    if (nlyr(a1) > 1) {a1 <- a1[[1]]}
    if (nlyr(a2) > 1) {a2 <- a2[[1]]}
    if (nlyr(a300) > 1) {a300 <- a300[[1]]}
    
    ########## 
    # STEP 5 # FIND REPRODUCTION RATE RASTERS
    ##########
    
    reproduction_files <- list.files(output_dir, pattern = "_reproductionRate\\.tif$",
                                     full.names = TRUE)
    
    message("Number of reproduction rate files found: ", length(reproduction_files))
    
    if (length(reproduction_files) == 0) {
      warning("No reproduction rate rasters found in: ",output_dir)
      next
    }
    
    
    reproduction_025 <- reproduction_files[grepl(paste0("^1_.*_Mammals001_", species_name, "_reproductionRate\\.tif$"), basename(reproduction_files))]
    reproduction_026 <- reproduction_files[grepl(paste0("^1_.*_Mammals002_", species_name, "_reproductionRate\\.tif$"), basename(reproduction_files))]
    reproduction_300 <- reproduction_files[grepl(paste0("^1_.*_Mammals300_", species_name, "_reproductionRate\\.tif$"), basename(reproduction_files))]
    
    if (length(reproduction_025) == 0) {
      warning("Could not find reproduction rate Mammals025 for ", species_name)
      next
    }
    
    if (length(reproduction_026) == 0) {
      warning("Could not find reproduction rate Mammals026 for ", species_name)
      next
    }
    
    if (length(reproduction_300) == 0) {
      warning("Could not find reproduction rate Mammals300 for ", species_name)
      next
    }
    
    reproduction_025 <- reproduction_025[1]
    reproduction_026 <- reproduction_026[1]
    reproduction_300 <- reproduction_300[1]
    
    message("Reproduction rate timestep 1: ", basename(reproduction_025))
    message("Reproduction rate timestep 2: ", basename(reproduction_026))
    message("Reproduction rate timestep 300: ", basename(reproduction_300))
    
    # load rasters
    rr1 <- rast(reproduction_025)
    rr2 <- rast(reproduction_026)
    rr300 <- rast(reproduction_300)
    
    if (nlyr(rr1) > 1) {rr1 <- rr1[[1]]}
    if (nlyr(rr2) > 1) {rr2 <- rr2[[1]]}
    if (nlyr(rr300) > 1) {rr300 <- rr300[[1]]}
    
    ##########
    # STEP 6 # Convert all rasters to dataframe
    ##########
    
    # suitability --------------------------------------------------------------
    df1 <- as.data.frame(r1, xy = TRUE, na.rm = TRUE)
    
    names(df1)[3] <- "value"
    df1$type <- "Suitability"
    df1$timestep <- paste0("Timestep ",timestep_1)
    
    df2 <- as.data.frame(r2, xy = TRUE,na.rm = TRUE)
    
    names(df2)[3] <- "value"
    df2$type <- "Suitability"
    df2$timestep <- paste0("Timestep ", timestep_2)
    
    df3 <- as.data.frame(r300, xy = TRUE,na.rm = TRUE)
    
    names(df3)[3] <- "value"
    df3$type <- "Suitability"
    df3$timestep <- paste0("Timestep ", timestep_300)
    
    # abundance ----------------------------------------------------------------
    df4 <- as.data.frame(a1, xy = TRUE, na.rm = TRUE)
    
    names(df4)[3] <- "value"
    df4$type <- "Abundance"
    df4$timestep <- paste0("Timestep ", timestep_1)
    
    df5 <- as.data.frame(a2, xy = TRUE, na.rm = TRUE)
    
    names(df5)[3] <- "value"
    df5$type <- "Abundance"
    df5$timestep <- paste0("Timestep ", timestep_2)
    
    df6 <- as.data.frame(a300, xy = TRUE, na.rm = TRUE)
    
    names(df6)[3] <- "value"
    df6$type <- "Abundance"
    df6$timestep <- paste0("Timestep ", timestep_300)
    
    # reproduction rate --------------------------------------------------------
    df7 <- as.data.frame(rr1, xy = TRUE, na.rm = TRUE)
    
    names(df7)[3] <- "value"
    df7$type <- "Reproduction rate"
    df7$timestep <- paste0("Timestep ", timestep_1)
    
    df8 <- as.data.frame(rr2, xy = TRUE, na.rm = TRUE)
    
    names(df8)[3] <- "value"
    df8$type <- "Reproduction rate"
    df8$timestep <- paste0("Timestep ", timestep_2)
    
    df9 <- as.data.frame(rr300, xy = TRUE, na.rm = TRUE)
    
    names(df9)[3] <- "value"
    df9$type <- "Reproduction rate"
    df9$timestep <- paste0("Timestep ", timestep_300)
    
    # combine everything into one
    df <- bind_rows(df1, df2, df3,
                    df4, df5, df6,
                    df7, df8, df9)
    
    message("Total plotting cells: ", nrow(df))
    
    ##########
    # STEP 7 # Plot
    ##########
    
    df$type <- factor(df$type,
                      levels = c("Suitability", "Abundance", "Reproduction rate"))
    
    df$timestep <- factor(df$timestep,
                          levels = c("Timestep 1", "Timestep 2", "Timestep 300"))
    
    # Get geographical background from run name
    region_map <- get_region_map(runname, terra::crs(r1))
    
    
    p <- ggplot() +
      # GEOGRAPHICAL BACKGROUND
      geom_sf(data = region_map, fill = "grey95", color = "grey70", linewidth = 0.25) +
      # SUITABILITY
      geom_raster(data = df[df$type == "Suitability", ], 
                  aes(x = x, y = y, fill = value))+
      scale_fill_viridis_c(option = "viridis", limits = c(0, 1), na.value = "transparent",
                           name = "Suitability",  guide = guide_colourbar(order = 1)) +
      ggnewscale::new_scale_fill() +
      # ABUNDANCE
      geom_raster(data = df[df$type == "Abundance", ], aes(x = x,y = y,fill = value)) +
      scale_fill_viridis_c(option = "viridis", na.value = "transparent", name = "Abundance",
                           guide = guide_colourbar(order = 2)) +
      ggnewscale::new_scale_fill() +
      # REPRODUCTION RATE
      geom_raster(data = df[df$type == "Reproduction rate", ], aes(x = x, y = y, fill = value)) +
      scale_fill_viridis_c(option = "viridis", na.value = "transparent", name = "Reproduction rate",
                           guide = guide_colourbar(order = 3)) +
      facet_grid(type ~ timestep) +
      #coord_equal(expand = FALSE) +
      theme_bw() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), "bottom", "right")) +
      labs(title = species_label)
    
    # print or save
    print(p)
    
    output_plot_dir <- diagnostics_dir
    
    ggsave(filename = file.path(output_plot_dir, paste0("SuitAbundRepRate_",species_label, runname, ".png")),
           plot = p,
           bg = "white",
           width = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 300, 210),
           height = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 150, 180),
           units = "mm",
           dpi = 1200
    )
    message("Plot saved successfully.")
    
    
    invisible(gc())
  }
}

########################################
# CHECKING WHAT HAPPENS IN THE BURN IN #
########################################

output_plot_dir <- diagnostics_dir

for (i in seq_len(nrow(run_paths))) {
  
  runname <- run_paths$runname[i]
  dir <- run_paths$output_dir[i]
  
  if (!dir.exists(dir)) next
  
  # STEP 1: Find all abundance rasters
  abundance_files <- list.files(dir, pattern = "_abundance\\.tif$", full.names = TRUE)
  
  if (length(abundance_files) == 0) next
  
  # Go through each species
  for (species_name in selected_sps) {
    
    species_file <- gsub(" ", ".", species_name)
    
    # Select timesteps 001 to 300
    species_abundance <- abundance_files[grepl(paste0("^1_.*_Mammals(00[1-9]|0[1-9][0-9]|[12][0-9]{2}|300)_", species_file, "_abundance\\.tif$"), basename(abundance_files))]
    
    if (length(species_abundance) == 0) {
      message("No abundance files found for ", species_name, " in ", runname)
      next
    }
    
    # Get timestep number from filename
    timestep <- str_extract(basename(species_abundance), "Mammals[0-9]{3}")
    timestep <- as.numeric(sub("Mammals", "", timestep))
    
    # Keep timestep 001, then every 12 timesteps to 300
    selected_timesteps <- c(1, seq(12, 300, by = 12))
    
    keep <- timestep %in% selected_timesteps
    
    species_abundance <- species_abundance[keep]
    timestep <- timestep[keep]
    
    if (length(species_abundance) == 0) {
      message("No selected timesteps found for ", species_name, " in ", runname)
      next
    }
    
    # Order timesteps correctly
    ord <- order(timestep)
    species_abundance <- species_abundance[ord]
    timestep <- timestep[ord]
    
    # Load and convert rasters to dataframes
    abundance_df <- list()
    
    for (j in seq_along(species_abundance)) {
      
      r <- rast(species_abundance[j])
      
      if (nlyr(r) > 1) {
        r <- r[[1]]
      }
      
      df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
      
      names(df)[3] <- "value"
      df$timestep <- sprintf("%03d", timestep[j])
      
      abundance_df[[j]] <- df
    }
    # Geographical background
    #region_map <- get_region_map(runname, terra::crs(r))
    
    abundance_df <- bind_rows(abundance_df)
    
    abundance_df$timestep <- factor(abundance_df$timestep, levels = sprintf("%03d", selected_timesteps))
    
    # Plot the 25 panels in a 5 x 5 grid
    p <- ggplot(abundance_df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_viridis_c(name = "Abundance", option = "viridis", na.value = "transparent") +
      facet_wrap(~timestep, ncol = 5, drop = FALSE) +
      coord_equal() +
      labs(title = species_name, x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    # Save figure per species
    safe_species <- gsub(" ", "_", species_name)
    
    output_file <- file.path(output_plot_dir,paste0("AbundEvolution_", runname, "_", safe_species, "_abundance_001_300.png"))
    
    ggsave(output_file,
      p, width = 12, height = ifelse( grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 6, 12),
      dpi = 300, bg = "white")
    
    message("Saved: ", output_file)
  }
}

###############################
# DISPERSAL CHANGES OVER TIME #
###############################

for (i in seq_len(nrow(run_paths))) {
  
  runname <- run_paths$runname[i]
  dir <- run_paths$output_dir[i]
  
  if (!dir.exists(dir)) next
  
  # Find all dispersal change rasters
  dispersal_files <- list.files(dir, pattern = "_dispersalChange\\.tif$", full.names = TRUE)
  
  if (length(dispersal_files) == 0) next
  
  # Go through all species
  for (species_name in selected_sps) {
    
    species_file <- gsub(" ", ".", species_name)
    
    # Select timesteps 001, 002, 200 and 300
    species_dispersal <- dispersal_files[grepl(paste0("^1_.*_Mammals(001|002|200|300)_", species_file, "_dispersalChange\\.tif$"), basename(dispersal_files))]
    
    if (length(species_dispersal) == 0) {
      message("No dispersal files found for ", species_name, " in ", runname)
      next
    }
    
    # Get timestep number from filename
    timestep <- str_extract(basename(species_dispersal), "Mammals[0-9]{3}")
    timestep <- as.numeric(sub("Mammals", "", timestep))
    
    # Order timesteps correctly
    ord <- order(timestep)
    species_dispersal <- species_dispersal[ord]
    timestep <- timestep[ord]
    
    # Load and convert rasters to dataframes
    dispersal_df <- list()
    
    for (j in seq_along(species_dispersal)) {
      
      r <- rast(species_dispersal[j])
      
      if (nlyr(r) > 1) {
        r <- r[[1]]
      }
      
      df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
      
      names(df)[3] <- "value"
      df$timestep <- sprintf("%03d", timestep[j])
      
      dispersal_df[[j]] <- df
    }
    
    dispersal_df <- bind_rows(dispersal_df)
    
    # Plot dispersal change in one row
    p <- ggplot(dispersal_df, aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_viridis_c(name = "Dispersal change", option = "viridis", na.value = "transparent") +
      facet_wrap(~timestep, ncol = 4, drop = FALSE) +
      coord_equal() +
      labs(title = species_name, x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    # Save figure per species
    safe_species <- gsub(" ", "_", species_name)
    
    output_file <- file.path(
      output_plot_dir,
      paste0("DispersalChangeEvolution_", runname, "_", safe_species, "_dispersalChange_001_300.png"))
    
    ggsave(output_file, p,
      width = 12, height = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 4, 4),
      dpi = 300, bg = "white")
    
    gc()
    message("Saved: ", output_file)
  }
}

message("✅ All diagnostic plots & files were created successfully!")