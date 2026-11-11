## checkup_dispersalDistanceMonthly.R ##
## Inês Silva ##
## 03 Nov 2026 ##
## GOAL? Before implementing a new formulation of the dispersal process that
## runs iteratively ove 12 monthls in a year, we need to check what happens to
## dispersal disatnces when species have small distances, e.g. 1 cell.


library(metaRange)
library(ggplot2)
library(dplyr)
library(purrr)
library(viridis)

# --------------------------------------------------
# 1. Annual dispersal parameters
# --------------------------------------------------

annual_mean <- 1
annual_max  <- 2

n_iterations <- 12

# --------------------------------------------------
# 2. Define scenarios
# --------------------------------------------------

scenarios <- data.frame(
  scenario = c(
    "Original annual",
    "Divide by 12",
    "Divide by sqrt(12)",
    "Mean unchanged, maximum = 1"
  ),
  
  mean_distance = c(
    annual_mean,
    annual_mean / 12,
    annual_mean / sqrt(n_iterations),
    annual_mean
  ),
  
  max_distance = c(
    annual_max,
    annual_max / 12,
    annual_max / sqrt(n_iterations),
    1
  )
)

# Maximum radius must be at least one cell
scenarios$max_radius <- pmax(
  1,
  ceiling(scenarios$max_distance)
)

print(scenarios)

# --------------------------------------------------
# 3. Calculate kernels
# --------------------------------------------------

kernels <- map(seq_len(nrow(scenarios)), function(i) {
  
  calculate_dispersal_kernel(
    max_dispersal_dist = as.integer(
      scenarios$max_radius[i]
    ),
    kfun = negative_exponential_function,
    mean_dispersal_dist = scenarios$mean_distance[i],
    normalize = TRUE
  )
  
})

names(kernels) <- scenarios$scenario

# --------------------------------------------------
# 4. Convert matrices into spatial coordinates
# --------------------------------------------------

kernel_df <- map2_dfr(
  kernels,
  names(kernels),
  function(k, scenario_name) {
    
    nr <- nrow(k)
    nc <- ncol(k)
    
    # Matrix indices
    coords <- expand.grid(
      row = seq_len(nr),
      col = seq_len(nc)
    )
    
    # Convert indices into distances from source cell
    coords$x <- coords$col - (nc + 1) / 2
    coords$y <- coords$row - (nr + 1) / 2
    
    # Matrix values (column-major order)
    coords$probability <- as.vector(k)
    
    coords %>%
      mutate(
        distance = sqrt(x^2 + y^2),
        scenario = scenario_name
      )
  }
)

# --------------------------------------------------
# 5. Plot 2D kernels
# --------------------------------------------------

p1 <- ggplot(
  kernel_df,
  aes(x = x, y = y, fill = probability)
) +
  geom_tile() +
  facet_wrap(
    ~scenario,
    scales = "fixed"
  ) +
  scale_fill_viridis_c(
    option = "C",
    name = "Probability"
  ) +
  coord_equal() +
  labs(
    title = "Comparison of dispersal kernels",
    subtitle = "Annual mean dispersal distance = 1 cell",
    x = "Distance from source (cells)",
    y = "Distance from source (cells)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p1)

# --------------------------------------------------
# 6. Plot probability along horizontal axis
# --------------------------------------------------

profile_df <- kernel_df %>%
  filter(abs(y) < 1e-10)

p2 <- ggplot(
  profile_df,
  aes(
    x = x,
    y = probability,
    colour = scenario,
    group = scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Dispersal probability along one axis",
    x = "Distance from source (cells)",
    y = "Probability per cell",
    colour = "Scenario"
  ) +
  theme_minimal(base_size = 12)

print(p2)

# --------------------------------------------------
# 7. Calculate kernel statistics
# --------------------------------------------------

kernel_summary <- kernel_df %>%
  group_by(scenario) %>%
  summarise(
    total_probability = sum(probability),
    
    mean_radial_distance = sum(
      distance * probability
    ),
    
    mean_squared_distance = sum(
      distance^2 * probability
    ),
    
    radial_variance = sum(
      (distance - mean_radial_distance)^2 *
        probability
    ),
    
    max_radius = max(
      distance[probability > 0]
    ),
    
    .groups = "drop"
  )

print(kernel_summary)
