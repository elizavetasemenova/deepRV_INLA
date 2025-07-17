
ls10_dir <- '../RINLA_data/ls_10'

grid_dirs <- list.dirs(ls10_dir, full.names = TRUE, recursive = FALSE)
grid_dirs <- grid_dirs[grepl('grid_', basename(grid_dirs))]

read_grid_data <- function(grid_dir) {
  list(
    obs_mask = read.csv(file.path(grid_dir, 'obs_mask.csv'), header = FALSE)[[1]],
    s = read.csv(file.path(grid_dir, 's.csv'), header = FALSE)[[1]],
    y_obs = read.csv(file.path(grid_dir, 'y_obs.csv'), header = FALSE)[[1]]
  )
}


grid_data <- setNames(lapply(grid_dirs, read_grid_data), basename(grid_dirs))

# ---  only the smallest grid (grid_256) ---
smallest_grid <- 'grid_256'
data <- grid_data[[smallest_grid]]

df <- data.frame(
  y_obs = data$y_obs,
  s = data$s,
  obs_mask = data$obs_mask
)

df <- df[df$obs_mask == TRUE, ]
rownames(df) <- NULL

print(str(df))

# fit: y_obs ~ s (Poisson)
library(INLA)
result <- inla(y_obs ~ s, data = df, family = 'poisson')

print(summary(result))
