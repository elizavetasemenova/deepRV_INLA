#https://inlabru-org.github.io/inlabru/articles/svc.html

# libraries
library(maps)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra) # raster plotting
library(tidyr)
library(scales)
library(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
# Note: the 'splancs' package also needs to be installed,
# but doesn't need to be loaded

# set option
select <- dplyr::select
options(scipen = 99999)
options(max.print = 99999)
options(stringsAsFactors = FALSE)


# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

# make a base map
states <- maps::map("state", plot = FALSE, fill = TRUE) %>%
  sf::st_as_sf() %>%
  filter(ID %in% c(
    "texas", "oklahoma", "kansas", "missouri",
    "arkansas", "louisiana"
  )) %>%
  sf::st_make_valid() %>%
  sf::st_transform(epsg6703km) %>%
  sf::st_make_valid()

robins_subset <- read.csv(paste0(
  "https://raw.github.com/tmeeha/inlaSVCBC",
  "/master/code/modeling_data.csv"
)) %>%
  select(
    circle, bcr, state, year, std_yr, count, log_hrs,
    lon, lat, obs
  ) %>%
  mutate(year = year + 1899) %>%
  filter(
    state %in% c(
      "TEXAS", "OKLAHOMA", "KANSAS", "MISSOURI",
      "ARKANSAS", "LOUISIANA"
    ),
    year >= 1987
  )




data(robins_subset)
count_dat <- robins_subset %>%
  mutate(site_idx = as.numeric(factor(paste(circle, lon, lat)))) %>%
  group_by(site_idx) %>%
  mutate(n_years = n()) %>%
  filter(n_years >= 20) %>%
  ungroup() %>%
  mutate(
    std_yr = year - max(year),
    obs = seq_len(nrow(.)),
    site_idx = as.numeric(factor(paste(circle, lon, lat))),
    year_idx = as.numeric(factor(year)),
    site_year_idx = as.numeric(factor(paste(circle, lon, lat, year)))
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(epsg6703km) %>%
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) %>%
  arrange(circle, year)

# map it
ggplot() +
  geom_sf(
    data = count_dat %>% filter(year_idx %in% seq(1, 30, 3)),
    aes(col = log(count + 1))
  ) +
  geom_sf(data = states, fill = NA) +
  coord_sf(datum = NA) +
  facet_wrap(~year) +
  scale_color_distiller(palette = "Spectral") +
  theme_bw()


# make a set of distinct study sites for mapping
site_map <- count_dat %>%
  select(circle, easting, northing) %>%
  distinct() %>%
  select(circle, easting, northing)



# make a two extension hulls and mesh for spatial model
hull <- fm_extensions(
  count_dat,
  convex = c(200, 500),
  concave = c(350, 500)
)
mesh <- fm_mesh_2d_inla(
  fm_hexagon_lattice(hull[[1]], edge_len = 70),
  boundary = hull,
  max.edge = c(100, 600), # km inside and outside
  cutoff = 25,
  offset = c(100, 300),
  crs = fm_crs(count_dat)
) # cutoff is min edge

# plot it
ggplot() +
  gg(data = mesh) +
  geom_sf(data = site_map, col = "darkgreen", size = 1) +
  geom_sf(data = states, fill = NA) +
  theme_bw() +
  labs(x = "", y = "")

# make spde
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(500, 0.5),
  prior.sigma = c(1, 0.5)
)

# iid prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))

# components
svc_components <- ~ -1 +
  kappa(site_idx, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
  alpha(geometry, model = spde) +
  eps(geometry, weights = log_hrs, model = spde) +
  tau(geometry, weights = std_yr, model = spde)

# formula, with "." meaning "add all the model components":
svc_formula <- count ~ .

res <- bru(
  svc_components,
  bru_obs(
    formula = svc_formula,
    family = "nbinomial",
    data = count_dat
  ),
  options = list(
    control.compute = list(waic = TRUE, cpo = FALSE),
    control.inla = list(
      int.strategy = "eb"
    ),
    verbose=TRUE
  )
)

# view results
res$summary.hyperpar[-1, c(1, 2)]

summary(exp(res$summary.random$alp$"0.5quant")) # exp(alpha) posterior median

summary(res$summary.random$eps$mean) # epsilon

summary((exp(res$summary.random$tau$"0.5quant") - 1) * 100) # (exp(tau)-1)*100

# get easting and northing limits
bbox <- fm_bbox(hull[[1]])
grd_dims <- round(c(x = diff(bbox[[1]]), y = diff(bbox[[2]])) / 25)

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = bbox[[1]], ylim = bbox[[2]], dims = grd_dims
)

# pull data
kappa <- data.frame(
  median = exp(res$summary.random$kappa$"0.5quant"),
  range95 = exp(res$summary.random$kappa$"0.975quant") -
    exp(res$summary.random$kappa$"0.025quant")
)
alph <- data.frame(
  median = exp(res$summary.random$alpha$"0.5quant"),
  range95 = exp(res$summary.random$alpha$"0.975quant") -
    exp(res$summary.random$alpha$"0.025quant")
)
epsi <- data.frame(
  median = res$summary.random$eps$"0.5quant",
  range95 = (res$summary.random$eps$"0.975quant" -
               res$summary.random$eps$"0.025quant")
)
taus <- data.frame(
  median = (exp(res$summary.random$tau$"0.5quant") - 1) * 100,
  range95 = (exp(res$summary.random$tau$"0.975quant") -
               exp(res$summary.random$tau$"0.025quant")) * 100
)

# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(alpha = alph, epsilon = epsi, tau = taus),
  function(x) as.matrix(fm_evaluate(mesh_proj, x))
)

# make a terra raster stack with the posterior median and range95
out_stk <- rast()
for (j in 1:3) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j, crs = epsg6703km)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j, crs = epsg6703km)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk) <- out_j
}
names(out_stk) <- c(
  "alpha_median", "alpha_range95", "epsilon_median",
  "epsilon_range95", "tau_median", "tau_range95"
)
out_stk <- terra::mask(
  out_stk,
  terra::vect(sf::st_union(states)),
  updatevalue = NA,
  touches = FALSE
)

make_plot_field <- function(data_stk, scale_label) {
  ggplot(states) +
    geom_sf(fill = NA) +
    coord_sf(datum = NA) +
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_distiller(
      scale_label,
      palette = "Spectral",
      na.value = "transparent"
    ) +
    theme_bw() +
    geom_sf(fill = NA)
}
make_plot_site <- function(data, scale_label) {
  ggplot(states) +
    geom_sf() +
    coord_sf(datum = NA) +
    geom_sf(data = data, size = 1, mapping = aes(colour = value)) +
    scale_colour_distiller(scale_label, palette = "Spectral") +
    labs(x = "", y = "") +
    theme_bw() +
    geom_sf(fill = NA)
}

# medians
# fields alpha_s, epsilon_s, tau_s
pa <- make_plot_field(
  data_stk = out_stk[["alpha_median"]],
  scale_label = "posterior\nmedian\nexp(alpha_s)"
)
pe <- make_plot_field(
  data_stk = out_stk[["epsilon_median"]],
  scale_label = "posterior\nmedian\nepsilon_s"
)
pt <- make_plot_field(
  data_stk = out_stk[["tau_median"]],
  scale_label = "posterior\nmedian\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$median)),
  scale_label = "posterior\nmedian\nexp(kappa_s)"
)
# range95
# fields alpha_s, epsilon_s, tau_s
pa_range95 <- make_plot_field(
  data_stk = out_stk[["alpha_range95"]],
  scale_label = "posterior\nrange95\nexp(alpha_s)"
)
pe_range95 <- make_plot_field(
  data_stk = out_stk[["epsilon_range95"]],
  scale_label = "posterior\nrange95\nepsilon_s"
)
pt_range95 <- make_plot_field(
  data_stk = out_stk[["tau_range95"]],
  scale_label = "posterior\nrange95\n100(exp(tau_s)-1)"
)
# sites kappa_s
ps_range95 <- make_plot_site(
  data = cbind(site_map, data.frame(value = kappa$range95)),
  scale_label = "posterior\nrange95\nexp(kappa_s)"
)

# plot together
multiplot(ps, pa, pe, pt, cols = 2)


# plot together
multiplot(ps_range95, pa_range95, pe_range95, pt_range95, cols = 2)
