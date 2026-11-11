## checkup_flyingDumbo.R ##
## Inês Silva ##
## 29 Oct 2026 ##
## Goal? Produce figures from model inputs and outputs to showcase how some
## species seem to be reaching certain areas that they shouldn't. The elefant
## is a very obvious example of this, as they are initiated in mainland africa
## and then cross the ocean (suitability = NA) and settle in Madagáscar.

source(file.path("src", "libraries.R"))
gc()
# Plot landscape going in as the metaRange environment

elef <- rast("E:/metaRange_May26/outputs/checkup_longBurnIn2/Africa_ssp126_20261003/Inputs/Loxodonta.africana_tropical_ssp126_cropped_reprojectedm.tif")

plot(elef$`2015`, main = "Loxodonta africana suitability map (Current)")

# plot species range
# load iucn's species ranges (to initiate species only within their range)
iucn <- vect("./data/externaldata/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
iucn$sci_name <- gsub(" ", ".", iucn$sci_name)
invisible(gc())

plot(iucn[iucn$sci_name == "Loxodonta.africana"])
gc()


# plot first timestep when species is initiated
library(terra)
library(sf)
library(rnaturalearth)

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(rnaturalearth)

# Load raster: first timestep when species is initiated
elef1 <- rast(
  "E:/metaRange_May26/outputs/checkup_longBurnIn2/Africa_ssp126_20261003/Outputs/1_TropicalSubtropicalMoistBroadleafForests_Africa_Mammals001_Loxodonta.africana_abundance.tif"
)

# Load Africa boundaries
africa <- ne_countries(
  continent = "Africa",
  scale = "small",
  returnclass = "sf"
)

# Select African elephant range
elef_range <- iucn[iucn$sci_name == "Loxodonta.africana", ]

# Match all layers to the raster CRS
africa <- st_transform(africa, crs(elef1))
elef_range <- project(elef_range, crs(elef1))


# Plot
ggplot()  +
  geom_spatraster(
    data = elef1
  ) +
  geom_sf(
    data = africa,
    fill = NA,
    colour = "grey60",
    linewidth = 0.2
  ) +
  scale_fill_viridis_c(
    name = "Abundance",
    na.value = "transparent"
  ) +
  coord_sf(
    crs = st_crs(elef1),
    expand = FALSE
  ) +
  labs(
    title = "Loxodonta africana",
    subtitle = "Initial abundance"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  )
