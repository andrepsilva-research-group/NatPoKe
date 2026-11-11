## Name: Checkup_suitAbundRepRate.R ##
## Author: Inês Silva ##
## Date: 02 October 2026 ##
## Description: Check what is happening to the suitability abundance and reproduction rates
## before and after the burn-in period. This is to check if they have too much
## suitable habitat outside the IUCN range, understand where sps are going and
## why we have explosions of population numbers. ##

##########
# STEP 1 # Choose a few selected species
##########

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

# chose timesteps to check
timestep_1 <- 1
timestep_2 <- 2
#timestep_100 <- 100
timestep_300 <- 300


# read the run table
runs <- read.csv("./data/run_table_longBurnIn.csv", stringsAsFactors = FALSE)

# input/output folders
inputFolder_paths <- character(0)
outputFolder_paths <- character(0)

for (i in seq_len(nrow(runs))) {
  
  target_region <- runs$region[i]
  target_biome <- runs$biome[i]
  future_scenario <- runs$scenario[i]
  
  runname <- paste(
    target_region,
    future_scenario,
    "20261103",
    sep = "_"
  )
  
  input_folder <- file.path("E:/metaRange_May26/outputs/checkup_newDispersal", runname, "Inputs")
  
  output_folder <- file.path("E:/metaRange_May26/outputs/checkup_newDispersal", runname, "Outputs")
  
  if (!file.exists(input_folder)) {
    
    warning("Input folder not found (skipping): ", input_folder)
    next
  }
  
  inputFolder_paths[runname] <- input_folder
  outputFolder_paths[runname] <- output_folder
}


# load IUCN ranges
iucn <- vect("./data/externaldata/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
iucn$sci_name <- gsub(" ", ".", iucn$sci_name)
invisible(gc())


# treshold from which we consider cells as suitable
threshold <- 0.6

selected_sps <- "Felis chaus"

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
    
    next
  }
  
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
    
    species_range <- iucn[
      iucn$sci_name == species_name,
    ]
    
    if (nrow(species_range) == 0) {
      
      warning("No IUCN range found for: ", species_name)
      next
    }
    
    # Project IUCN range to suitability CRS
    species_range <- project(species_range, crs(r1))
    # Crop + mask suitability
    r1 <- crop(r1, species_range)
    r1 <- mask(r1, species_range)
    
    
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
    
    abundance_files <- list.files(
      output_dir,
      pattern = "_abundance\\.tif$",
      full.names = TRUE
    )
    
    if (length(abundance_files) == 0) {
      
      warning(
        "No abundance rasters found in: ",
        output_dir
      )
      
      next
    }
    # keep only selected timesteps
    abundance_025 <- abundance_files[
      grepl(paste0("^2_.*_Mammals001_", species_name, "_abundance\\.tif$"),
        basename(abundance_files))
    ]
    
    abundance_026 <- abundance_files[
      grepl(paste0("^2_.*_Mammals002_", species_name, "_abundance\\.tif$"),
        basename(abundance_files))
    ]
    
    
    abundance_300 <- abundance_files[
      grepl(paste0("^2_.*_Mammals300_", species_name, "_abundance\\.tif$"),
            basename(abundance_files))
    ]
    
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
    
    reproduction_files <- list.files(
      output_dir,
      pattern = "_reproductionRate\\.tif$",
      full.names = TRUE
    )
    
    message(
      "Number of reproduction rate files found: ",
      length(reproduction_files)
    )
    
    if (length(reproduction_files) == 0) {
      
      warning("No reproduction rate rasters found in: ",output_dir)
      next
    }
    

    reproduction_025 <- reproduction_files[
      grepl(paste0("^2_.*_Mammals001_", species_name, "_reproductionRate\\.tif$"),
        basename(reproduction_files))
    ]
    
    reproduction_026 <- reproduction_files[
      grepl(paste0("^2_.*_Mammals002_", species_name, "_reproductionRate\\.tif$"),
        basename(reproduction_files))
    ]
    
    reproduction_300 <- reproduction_files[
      grepl(paste0("^2_.*_Mammals300_", species_name, "_reproductionRate\\.tif$"),
            basename(reproduction_files))
    ]
    
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
    df <- bind_rows(
      df1,
      df2,
      df3,
      df4,
      df5,
      df6,
      df7,
      df8,
      df9
    )
    
    message("Total plotting cells: ", nrow(df))
    
    ##########
    # STEP 7 # Plot
    ##########
    
    df$type <- factor(df$type,
      levels = c("Suitability", "Abundance", "Reproduction rate"))
    
    df$timestep <- factor(df$timestep,
                      levels = c("Timestep 1", "Timestep 2", "Timestep 300"))
    
    
    p <- ggplot() +
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
      coord_equal(expand = FALSE) +
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
    
    output_plot_dir <- file.path("E:/metaRange_May26/outputs/checkup_newDispersal/diagnostics")
    
    ggsave(filename = file.path(output_plot_dir, paste0("SuitAbundRepRate2_",species_label, runname, ".png")),
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


################################################################################
####################### EXTRA CHECK UP FOR BURN IN STEPS #######################
################################################################################

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

inputFolder_paths <- list.dirs(
  "E:/metaRange_May26/outputs/checkup_newDispersal",
  recursive = FALSE
)

output_plot_dir <- "E:/metaRange_May26/outputs/checkup_newDispersal/diagnostics"

for (dir in inputFolder_paths) {
  
  runname <- basename(dir)
  
  output_dir <- file.path(dir, "Outputs")
  
  if (!dir.exists(output_dir)) next
  
  # -----------------------------------------------------
  # FIND ALL ABUNDANCE RASTERS
  # -----------------------------------------------------
  
  abundance_files <- list.files(
    output_dir,
    pattern = "_abundance\\.tif$",
    full.names = TRUE
  )
  
  if (length(abundance_files) == 0) next
  
  # -----------------------------------------------------
  # LOOP THROUGH SPECIES
  # -----------------------------------------------------
  
  for (species_name in selected_sps) {
    
    species_file <- gsub(" ", ".", species_name)
    
    # ---------------------------------------------------
    # FIND ABUNDANCE FILES FOR TIMESTEPS 001–300
    # ---------------------------------------------------
    
    species_abundance <- abundance_files[
      grepl(
        paste0(
          "^2_.*_Mammals(00[1-9]|0[1-9][0-9]|[12][0-9]{2}|300)_",
          species_file,
          "_abundance\\.tif$"
        ),
        basename(abundance_files)
      )
    ]
    
    if (length(species_abundance) == 0) {
      message(
        "No abundance files found for ",
        species_name,
        " in ",
        runname
      )
      next
    }
    
    # ---------------------------------------------------
    # EXTRACT TIMESTEP FROM FILENAME
    # ---------------------------------------------------
    
    timestep <- str_extract(
      basename(species_abundance),
      "Mammals[0-9]{3}"
    )
    
    timestep <- as.numeric(sub("Mammals", "", timestep))
    
    # ---------------------------------------------------
    # KEEP TIMESTEP 001, THEN EVERY 12 TIMESTEPS TO 300
    # Total: 25 panels
    # ---------------------------------------------------
    
    selected_timesteps <- c(1, seq(12, 300, by = 12))
    
    keep <- timestep %in% selected_timesteps
    
    species_abundance <- species_abundance[keep]
    timestep <- timestep[keep]
    
    if (length(species_abundance) == 0) {
      message(
        "No selected timesteps found for ",
        species_name,
        " in ",
        runname
      )
      next
    }
    
    # ---------------------------------------------------
    # ORDER CORRECTLY
    # ---------------------------------------------------
    
    ord <- order(timestep)
    
    species_abundance <- species_abundance[ord]
    timestep <- timestep[ord]
    
    # ---------------------------------------------------
    # LOAD AND CONVERT RASTERS
    # ---------------------------------------------------
    
    abundance_df <- list()
    
    for (i in seq_along(species_abundance)) {
      
      r <- rast(species_abundance[i])
      
      if (nlyr(r) > 1) {
        r <- r[[1]]
      }
      
      df <- as.data.frame(
        r,
        xy = TRUE,
        na.rm = TRUE
      )
      
      names(df)[3] <- "value"
      
      df$timestep <- sprintf("%03d", timestep[i])
      
      abundance_df[[i]] <- df
    }
    
    abundance_df <- bind_rows(abundance_df)
    
    # ---------------------------------------------------
    # MAKE TIMESTEP A FACTOR IN CORRECT ORDER
    # ---------------------------------------------------
    
    abundance_df$timestep <- factor(
      abundance_df$timestep,
      levels = sprintf("%03d", selected_timesteps)
    )
    
    # ---------------------------------------------------
    # PLOT: 25 PANELS, 5 × 5 GRID
    # ---------------------------------------------------
    
    p <- ggplot(
      abundance_df,
      aes(x = x, y = y, fill = value)
    ) +
      geom_raster() +
      scale_fill_viridis_c(
        name = "Abundance",
        option = "viridis",
        na.value = "transparent"
      ) +
      facet_wrap(
        ~ timestep,
        ncol = 5,
        drop = FALSE
      ) +
      coord_equal() +
      labs(
        title = species_name,
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    # ---------------------------------------------------
    # SAVE FIGURE
    # ---------------------------------------------------
    
    safe_species <- gsub(" ", "_", species_name)
    
    output_file <- file.path(
      output_plot_dir,
      paste0(
        "AbundEvolution_",
        runname,
        "_",
        safe_species,
        "_abundance_001_300.png"
      )
    )
    
    ggsave(
      output_file,
      p,
      width = 12,
      height = ifelse(
        grepl(
          "Europe\\+Asia|North America",
          runname,
          ignore.case = TRUE
        ),
        8,
        12
      ),
      dpi = 300,
      bg = "white"
    )
    
    message("Saved: ", output_file)
  }
}

#########################
### DISPERSAL CHANGES ###
#########################

for (dir in inputFolder_paths) {
  
  runname <- basename(dir)
  
  output_dir <- file.path(dir, "Outputs")
  
  if (!dir.exists(output_dir)) next
  
  # -----------------------------------------------------
  # FIND ALL ABUNDANCE RASTERS
  # -----------------------------------------------------
  
  dispersal_files <- list.files(output_dir, pattern = "_dispersalChange\\.tif$", full.names = TRUE)
  
  if (length(dispersal_files) == 0) next
  
  # -----------------------------------------------------
  # LOOP THROUGH SPECIES
  # -----------------------------------------------------
  
  for (species_name in selected_sps) {
    
    species_file <- gsub(" ", ".", species_name)
    
    # ---------------------------------------------------
    # FIND DISPERSAL CHANGES FILES FOR TIMESTEPS 001–300
    # ---------------------------------------------------
    
    species_dispersal <- dispersal_files[
      grepl(
        paste0(
          "^1_.*_Mammals(001|002|200|300)_",
          species_file,
          "_dispersalChange\\.tif$"
        ),
        basename(dispersal_files)
      )
    ]
    
    if (length(species_dispersal) == 0) {
      message(
        "No abundance files found for ",
        species_name,
        " in ",
        runname
      )
      next
    }
    
    # ---------------------------------------------------
    # EXTRACT TIMESTEP FROM FILENAME
    # ---------------------------------------------------
    
    timestep <- str_extract(
      basename(species_dispersal),
      "Mammals[0-9]{3}"
    )
    
    timestep <- as.numeric(sub("Mammals", "", timestep))
    
    # ---------------------------------------------------
    # KEEP TIMESTEP 001, THEN EVERY 12 TIMESTEPS TO 300
    # Total: 25 panels
    # ---------------------------------------------------
    
    #selected_timesteps <- c(1, seq(12, 300, by = 12))
    
    #keep <- timestep %in% selected_timesteps
    
    #species_dispersal <- species_dispersal[keep]
    #timestep <- timestep[keep]
    
    if (length(species_dispersal) == 0) {
      message(
        "No selected timesteps found for ",
        species_name,
        " in ",
        runname
      )
      next
    }
    
    # ---------------------------------------------------
    # ORDER CORRECTLY
    # ---------------------------------------------------
    
    ord <- order(timestep)
    
    species_dispersal <- species_dispersal[ord]
    timestep <- timestep[ord]
    
    # ---------------------------------------------------
    # LOAD AND CONVERT RASTERS
    # ---------------------------------------------------
    
    dispersal_df <- list()
    
    for (i in seq_along(species_dispersal)) {
      
      r <- rast(species_dispersal[i])
      
      if (nlyr(r) > 1) {
        r <- r[[1]]
      }
      
      df <- as.data.frame(
        r,
        xy = TRUE,
        na.rm = TRUE
      )
      
      names(df)[3] <- "value"
      
      df$timestep <- sprintf("%03d", timestep[i])
      
      dispersal_df[[i]] <- df
    }
    
    dispersal_df <- bind_rows(dispersal_df)
    
    # ---------------------------------------------------
    # MAKE TIMESTEP A FACTOR IN CORRECT ORDER
    # ---------------------------------------------------
    
    # dispersal_df$timestep <- factor(
    #   dispersal_df$timestep,
    #   levels = sprintf("%03d", selected_timesteps)
    # )
    
    # ---------------------------------------------------
    # PLOT: 25 PANELS, 5 × 5 GRID
    # ---------------------------------------------------
    
    p <- ggplot(
      dispersal_df,
      aes(x = x, y = y, fill = value)
    ) +
      geom_raster() +
      scale_fill_viridis_c(
        name = "Abundance",
        option = "viridis",
        na.value = "transparent"
      ) +
      facet_wrap(
        ~ timestep,
        ncol = 4,
        drop = FALSE
      ) +
      coord_equal() +
      labs(
        title = species_name,
        x = NULL,
        y = NULL
      ) +
      theme_minimal() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    p
    
    # ---------------------------------------------------
    # SAVE FIGURE
    # ---------------------------------------------------
    
    safe_species <- gsub(" ", "_", species_name)
    
    output_file <- file.path(
      output_plot_dir,
      paste0(
        "DispersalChangeEvolution_",
        runname,
        "_",
        safe_species,
        "_001_300.png"
      )
    )
    
    ggsave(output_file, p, width = 12, height = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 2, 4), dpi = 300, bg = "white")
    gc()
    message("Saved: ", output_file)
  }
}


# for (dir in inputFolder_paths) {
#   
#   runname <- basename(dir)
#   
#   output_dir <- file.path(
#     dir,
#     "Outputs"
#   )
#   
#   if (!dir.exists(output_dir)) next
#   
#   # -----------------------------------------------------
#   # FIND ALL ABUNDANCE RASTERS
#   # -----------------------------------------------------
#   
#   abundance_files <- list.files(
#     output_dir,
#     pattern = "_abundance\\.tif$",
#     full.names = TRUE
#   )
#   
#   if (length(abundance_files) == 0) next
#   
#   # -----------------------------------------------------
#   # LOOP THROUGH SPECIES
#   # -----------------------------------------------------
#   
#   for (species_name in selected_sps) {
#     
#     # Convert species to filename format
#     species_file <- gsub(
#       " ",
#       ".",
#       species_name
#     )
#     
#     # ---------------------------------------------------
#     # FIND ABUNDANCE FILES 001–025
#     # ---------------------------------------------------
#     
#     species_abundance <- abundance_files[
#       grepl(
#         paste0(
#           "^1_.*_Mammals(00[1-9]|0[1-9][0-9]|100)_",
#           species_file,
#           "_abundance\\.tif$"
#         ),
#         basename(abundance_files)
#       )
#     ]
#     
#     if (length(species_abundance) == 0) {
#       message(
#         "No abundance files found for ",
#         species_name,
#         " in ",
#         runname
#       )
#       next
#     }
#     
#     # ---------------------------------------------------
#     # EXTRACT TIMESTEP FROM FILENAME
#     # ---------------------------------------------------
#     
#     timestep <- str_extract(basename(species_abundance), "Mammals[0-9]{3}")
#     
#     timestep <- as.numeric(sub("Mammals", "", timestep))
#     
#     # ---------------------------------------------------
#     # KEEP ONLY:
#     # 001, 005, 010, 015, ..., 100
#     # ---------------------------------------------------
#     
#     keep <- timestep == 1 |
#       (timestep >= 5 &
#          timestep <= 100 &
#          timestep %% 5 == 0)
#     
#     species_abundance <- species_abundance[keep]
#     timestep <- timestep[keep]
#     
#     # ---------------------------------------------------
#     # ORDER CORRECTLY
#     # ---------------------------------------------------
#     
#     ord <- order(timestep)
#     
#     species_abundance <- species_abundance[ord]
#     timestep <- timestep[ord]
#     
#     # ---------------------------------------------------
#     # LOAD AND CONVERT RASTERS
#     # ---------------------------------------------------
#     
#     abundance_df <- list()
#     
#     for (i in seq_along(species_abundance)) {
#       
#       r <- rast(species_abundance[i])
#       
#       if (nlyr(r) > 1) {
#         r <- r[[1]]
#       }
#       
#       df <- as.data.frame(
#         r,
#         xy = TRUE,
#         na.rm = TRUE
#       )
#       
#       names(df)[3] <- "value"
#       
#       df$timestep <- sprintf(
#         "%03d",
#         timestep[i]
#       )
#       
#       abundance_df[[i]] <- df
#     }
#     
#     abundance_df <- bind_rows(
#       abundance_df
#     )
#     
#     # ---------------------------------------------------
#     # MAKE TIMESTEP A FACTOR IN CORRECT ORDER
#     # ---------------------------------------------------
#     
#     selected_timesteps <- c(1, seq(5, 100, by = 5))
#     
#     abundance_df$timestep <- factor(abundance_df$timestep,
#                                     levels = sprintf("%03d", selected_timesteps))
#     
#     # ---------------------------------------------------
#     # PLOT
#     # ---------------------------------------------------
#     
#     p <- ggplot(abundance_df, aes(x = x, y = y, fill = value)) +
#       geom_raster() +
#       scale_fill_viridis_c(name = "Abundance", option = "viridis", na.value = "transparent") +
#       facet_wrap(~ timestep, ncol = 5) +
#       coord_equal() +
#       labs(title = species_name, #subtitle = runname,
#         x = NULL, y = NULL) +
#       theme_minimal() +
#       theme(
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid = element_blank(),
#         strip.text = element_text(face = "bold"),
#         plot.title = element_text(face = "bold"),
#         legend.position = "right")
#     
#     safe_species <- gsub(" ", "_", species_name)
#     
#     output_file <- file.path(output_plot_dir,
#       paste0("AbundEvolution_", runname,"_", safe_species,"_abundance_001_025.png")
#     )
#     
#     ggsave(output_file,
#            p,
#            width = 12,
#            height = ifelse(grepl("Europe\\+Asia|North America", runname, ignore.case = TRUE), 5, 10),
#            dpi = 300,  bg = "white")
#     
#     message("Saved: ", output_file)
#   }
# }
