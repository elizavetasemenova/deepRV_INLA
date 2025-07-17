library(INLA)
library(ggplot2)
library(gridExtra)
library(viridis)
library(dplyr)

# data
obs_mask <- read.csv('../RINLA_data/2d/obs_mask.csv', header = TRUE)
s <- read.csv('../RINLA_data/2d/s.csv', header = TRUE)
y_obs <- read.csv('../RINLA_data/2d/y_obs.csv', header = TRUE)

df <- data.frame(
  obs_mask = as.numeric(obs_mask[[1]]),
  x1 = as.numeric(s[[1]]),
  x2 = as.numeric(s[[2]]),
  y_obs = as.numeric(y_obs[[1]])
)

# split data into observed and held-out
df_obs <- df[df$obs_mask == 1, ]
df_heldout <- df[df$obs_mask == 0, ]

cat("Data summary:\n")
cat("Total observations:", nrow(df), "\n")
cat("Observed points:", nrow(df_obs), "\n")
cat("Held-out points:", nrow(df_heldout), "\n")

# build mesh on all locations (observed + held-out) 
all_coords <- as.matrix(df[, c("x1", "x2")])
coords <- as.matrix(df_obs[, c("x1", "x2")])

data_range <- apply(all_coords, 2, range)
data_span <- max(data_range[2,] - data_range[1,])
convex_param <- max(0.1, data_span * 0.01)

mesh <- inla.mesh.2d(
  loc = coords,
  boundary = inla.nonconvex.hull(all_coords, convex = convex_param),
  max.edge = c(5, 20), 
  cutoff = 1
)
#TODO: checkout this warning

cat("Mesh created with", mesh$n, "nodes\n")
cat("Convex parameter used:", round(convex_param, 3), "\n")


p1 <- ggplot(df, aes(x = x1, y = x2)) +
  geom_point(aes(color = factor(obs_mask)), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("0" = "red", "1" = "blue"), 
                     labels = c("Held-out", "Observed"),
                     name = "Data") +
  labs(title = 'Data locations', x = 'x1', y = 'x2') +
  theme_minimal() +
  theme(legend.position = "bottom")

mesh_tri <- mesh$graph$tv
mesh_coords <- mesh$loc
mesh_lines <- do.call(rbind, lapply(1:nrow(mesh_tri), function(i) {
  tri <- mesh_tri[i,]
  data.frame(
    x = mesh_coords[c(tri, tri[1]), 1],
    y = mesh_coords[c(tri, tri[1]), 2],
    group = i
  )
}))

p2 <- ggplot() +
  geom_path(data = mesh_lines, aes(x = x, y = y, group = group), 
            color = 'black', alpha = 0.3, linewidth = 0.2) +
  geom_point(data = df_obs, aes(x = x1, y = x2), 
             color = 'blue', size = 1, alpha = 0.8) +
  labs(title = 'INLA mesh with observed data', x = 'x1', y = 'x2') +
  theme_minimal()

# plots before model fitting
grid.arrange(p1, p2, ncol = 2)


cat("\nFitting INLA model...\n")

# SPDE model (Matern 3/2 kernel, nu=1)
beta_shape1 <- 4
beta_shape2 <- 1
beta_median <- qbeta(0.5, beta_shape1, beta_shape2)
range_median <- exp(beta_median * log(100))
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(range_median, 0.5),
  prior.sigma = c(1, 0.99)
)

# projector matrix for observed
a_obs <- inla.spde.make.A(mesh, loc = coords)
s_index <- inla.spde.make.index("s", n.spde = spde$n.spde)
effects <- data.frame(s = s_index$s, intercept = 1)

stack_obs <- inla.stack(
  data = list(y = df_obs$y_obs),
  A = list(a_obs),
  effects = list(effects),
  tag = "obs"
)

formula <- y ~ 0 + intercept + f(s, model = spde)

# fit model with proper control settings
start_time <- Sys.time()
result <- inla(
  formula,
  data = inla.stack.data(stack_obs),
  family = "poisson",
  control.predictor = list(
    A = inla.stack.A(stack_obs), 
    compute = TRUE,
    link = 1  # TODO: log link for Poisson!!!
  ),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
)
inference_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("Inference time (seconds):", inference_time, "\n")

# --- evaluation on held-out data only ---
if (nrow(df_heldout) > 0) {
  cat("\nPredicting on held-out data...\n")
  
  coords_heldout <- as.matrix(df_heldout[, c("x1", "x2")])
  a_heldout <- inla.spde.make.A(mesh, loc = coords_heldout)
  
  # check for zero-weight rows (points outside mesh coverage)
  zero_rows <- which(rowSums(a_heldout) == 0)
  if (length(zero_rows) > 0) {
    cat("Warning:", length(zero_rows), "held-out points are outside mesh coverage and will be excluded\n")
    
    # remove points with zero weights
    valid_indices <- setdiff(1:nrow(df_heldout), zero_rows)
    df_heldout_valid <- df_heldout[valid_indices, ]
    coords_heldout_valid <- coords_heldout[valid_indices, , drop = FALSE]
    a_heldout_valid <- a_heldout[valid_indices, , drop = FALSE]
    
    if (nrow(df_heldout_valid) == 0) {
      cat("Error: No valid held-out points within mesh coverage\n")
    } else {
      cat("Using", nrow(df_heldout_valid), "valid held-out points for evaluation\n")
      
      stack_heldout <- inla.stack(
        data = list(y = NA),
        A = list(a_heldout_valid),
        effects = list(effects),
        tag = "heldout"
      )
      
      stack_full <- inla.stack(stack_obs, stack_heldout)
      
      result_full <- inla(
        formula,
        data = inla.stack.data(stack_full),
        family = "poisson",
        control.predictor = list(
          A = inla.stack.A(stack_full), 
          compute = TRUE,
          link = 1  # TODO: log link for Poisson!!!
        ),
        control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
      )
      
      df_heldout <- df_heldout_valid
    }
  } else {
    stack_heldout <- inla.stack(
      data = list(y = NA),
      A = list(a_heldout),
      effects = list(effects),
      tag = "heldout"
    )
    
    stack_full <- inla.stack(stack_obs, stack_heldout)
    
    result_full <- inla(
      formula,
      data = inla.stack.data(stack_full),
      family = "poisson",
      control.predictor = list(
        A = inla.stack.A(stack_full), 
        compute = TRUE,
        link = 1  # Use log link for Poisson
      ),
      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
    )
  }
  
  # held-out predictions
  idx_heldout <- inla.stack.index(stack_full, tag = "heldout")$data
  y_hat_heldout <- result_full$summary.fitted.values$mean[idx_heldout]
  y_hat_heldout_sd <- result_full$summary.fitted.values$sd[idx_heldout]
  
  # evaluation metrics on held-out data
  mse_heldout <- mean((df_heldout$y_obs - y_hat_heldout)^2)
  rmse_heldout <- sqrt(mse_heldout)
  mae_heldout <- mean(abs(df_heldout$y_obs - y_hat_heldout))
  
  # R-squared
  ss_res <- sum((df_heldout$y_obs - y_hat_heldout)^2)
  ss_tot <- sum((df_heldout$y_obs - mean(df_heldout$y_obs))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  
  cat("\n=== HELD-OUT DATA EVALUATION ===\n")
  cat("MSE:", round(mse_heldout, 4), "\n")
  cat("RMSE:", round(rmse_heldout, 4), "\n")
  cat("MAE:", round(mae_heldout, 4), "\n")
  cat("R-squared:", round(r_squared, 4), "\n")
  
  results_df <- data.frame(
    x1 = df_heldout$x1,
    x2 = df_heldout$x2,
    y_obs = df_heldout$y_obs,
    y_hat = y_hat_heldout,
    y_hat_sd = y_hat_heldout_sd,
    residual = df_heldout$y_obs - y_hat_heldout
  )
  
  # Save results
  write.csv(results_df, file = "../samples/inla_heldout_results.csv", row.names = FALSE)
  
  # --- visualise results  ---
  
  x1_range <- range(df$x1)
  x2_range <- range(df$x2)
  grid_res <- 50
  
  x1_seq <- seq(x1_range[1], x1_range[2], length.out = grid_res)
  x2_seq <- seq(x2_range[1], x2_range[2], length.out = grid_res)
  pred_grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)
  
  # projection matrix for grid
  coords_grid <- as.matrix(pred_grid)
  a_grid <- inla.spde.make.A(mesh, loc = coords_grid)
  
  # valid grid points
  valid_grid_indices <- which(rowSums(a_grid) > 0)
  if (length(valid_grid_indices) > 0) {
    pred_grid_valid <- pred_grid[valid_grid_indices, ]
    a_grid_valid <- a_grid[valid_grid_indices, , drop = FALSE]
    
    # stack for grid prediction
    stack_grid <- inla.stack(
      data = list(y = NA),
      A = list(a_grid_valid),
      effects = list(effects),
      tag = "grid"
    )
    
    # combine all stacks
    stack_all <- inla.stack(stack_obs, stack_heldout, stack_grid)
    
    # predict on full grid with proper link function
    result_grid <- inla(
      formula,
      data = inla.stack.data(stack_all),
      family = "poisson",
      control.predictor = list(
        A = inla.stack.A(stack_all), 
        compute = TRUE,
        link = 1  # TODO: use log link for Poisson!!!
      ),
      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
    )
    
    # grid predictions
    idx_grid <- inla.stack.index(stack_all, tag = "grid")$data
    grid_predictions <- result_grid$summary.fitted.values$mean[idx_grid]
    grid_sd <- result_grid$summary.fitted.values$sd[idx_grid]
    
    grid_df <- data.frame(
      x1 = pred_grid_valid$x1,
      x2 = pred_grid_valid$x2,
      prediction = grid_predictions,
      prediction_sd = grid_sd
    )
    
    # update held-out predictions from grid results
    idx_heldout_grid <- inla.stack.index(stack_all, tag = "heldout")$data
    y_hat_heldout <- result_grid$summary.fitted.values$mean[idx_heldout_grid]
    y_hat_heldout_sd <- result_grid$summary.fitted.values$sd[idx_heldout_grid]
    
    # recalculate metrics with updated predictions
    mse_heldout <- mean((df_heldout$y_obs - y_hat_heldout)^2)
    rmse_heldout <- sqrt(mse_heldout)
    mae_heldout <- mean(abs(df_heldout$y_obs - y_hat_heldout))
    
    ss_res <- sum((df_heldout$y_obs - y_hat_heldout)^2)
    ss_tot <- sum((df_heldout$y_obs - mean(df_heldout$y_obs))^2)
    r_squared <- 1 - (ss_res / ss_tot)
    
    cat("\n=== UPDATED HELD-OUT DATA EVALUATION (with grid predictions) ===\n")
    cat("MSE:", round(mse_heldout, 4), "\n")
    cat("RMSE:", round(rmse_heldout, 4), "\n")
    cat("MAE:", round(mae_heldout, 4), "\n")
    cat("R-squared:", round(r_squared, 4), "\n")
    
  } else {
    # fallback to original prediction method
    grid_df <- NULL
  }
  
  #  vis
  theme_results <- theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.position = "right"
    )
  
  r1 <- ggplot(results_df, aes(x = y_obs, y = y_hat)) +
    geom_point(alpha = 0.7, size = 2, color = "steelblue") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "darkblue", linewidth = 0.8) +
    labs(title = paste("Predicted vs observed (R² =", round(r_squared, 3), ")"),
         x = "observed", y = "predicted") +
    theme_results
  
  r2 <- ggplot(results_df, aes(x = y_hat, y = residual)) +
    geom_point(alpha = 0.7, size = 2, color = "steelblue") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    geom_smooth(method = "loess", se = FALSE, color = "darkblue", linewidth = 0.8) +
    labs(title = "Residuals vs fitted", x = "Fitted values", y = "Residuals") +
    theme_results
  
  # spatial field prediction 
  if (!is.null(grid_df)) {
    r3 <- ggplot() +
      geom_tile(data = grid_df, aes(x = x1, y = x2, fill = prediction), alpha = 0.8) +
      geom_point(data = df_obs, aes(x = x1, y = x2), color = "white", size = 1, alpha = 0.8) +
      geom_point(data = df_heldout, aes(x = x1, y = x2), color = "black", size = 1, alpha = 0.8) +
      scale_fill_viridis_c(name = "Predicted\nvalue") +
      labs(title = "spatial field predictions\n(white: observed, black: held-out)", 
           x = "x1", y = "x2") +
      theme_results +
      coord_equal()
  } else {
    r3 <- ggplot(results_df, aes(x = x1, y = x2, color = y_hat)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_viridis_c(name = "Predicted\nvalue") +
      labs(title = "held-out predictions", x = "x1", y = "x2") +
      theme_results +
      coord_equal()
  }
  

  if (!is.null(grid_df)) {
    r4 <- ggplot() +
      geom_tile(data = grid_df, aes(x = x1, y = x2, fill = prediction_sd), alpha = 0.8) +
      geom_point(data = df_obs, aes(x = x1, y = x2), color = "white", size = 1, alpha = 0.8) +
      geom_point(data = df_heldout, aes(x = x1, y = x2), color = "black", size = 1, alpha = 0.8) +
      scale_fill_viridis_c(name = "Prediction\nSD", option = "plasma") +
      labs(title = "prediction uncertainty\n(white: observed, black: held-out)", 
           x = "x1", y = "x2") +
      theme_results +
      coord_equal()
  } else {
    r4 <- ggplot(results_df, aes(x = x1, y = x2, color = y_hat_sd)) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_viridis_c(name = "Prediction\nSD", option = "plasma") +
      labs(title = "prediction uncertainty", x = "x1", y = "x2") +
      theme_results +
      coord_equal()
  }
  
  r5 <- ggplot(results_df, aes(x = x1, y = x2, color = residual)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, name = "Residual") +
    labs(title = "spatial residuals", x = "x1", y = "x2") +
    theme_results +
    coord_equal()
  
  r6 <- ggplot(results_df, aes(sample = residual)) +
    stat_qq(color = "steelblue", alpha = 0.7) +
    stat_qq_line(color = "red", linewidth = 1) +
    labs(title = "Q-Q plot of residuals", x = "theoretical quantiles", y = "sample quantiles") +
    theme_results
  
  
  cat("\n=== DISPLAYING RESULTS ===\n")
  
  grid.arrange(r1, r2, ncol = 2)
  
  if (!is.null(grid_df)) {
    grid.arrange(r3, r4, ncol = 2)
  } else {
    grid.arrange(r3, r4, ncol = 2)
  }
  
  grid.arrange(r5, r6, ncol = 2)
  
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("Mean observed value:", round(mean(df_heldout$y_obs), 3), "\n")
  cat("Mean predicted value:", round(mean(y_hat_heldout), 3), "\n")
  cat("Correlation coefficient:", round(cor(df_heldout$y_obs, y_hat_heldout), 3), "\n")
  cat("Mean prediction uncertainty (SD):", round(mean(y_hat_heldout_sd), 3), "\n")
  
} else {
  cat("No held-out data available for evaluation.\n")
}

cat("\n=== MODEL SUMMARY ===\n")
cat("DIC:", round(result$dic$dic, 2), "\n")
cat("WAIC:", round(result$waic$waic, 2), "\n")
cat("Number of mesh nodes:", mesh$n, "\n")
cat("Fixed effects:\n")
print(result$summary.fixed)

