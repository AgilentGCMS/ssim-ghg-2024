chol_solve <- function(R, b) {
  backsolve(R, backsolve(R, b, transpose = TRUE))
}

rmvnorm <- function(mu, S, chol_S, S_inv, chol_S_inv) {
  if (!missing(S) || !missing(chol_S)) {
    if (missing(chol_S)) {
      chol_S <- chol(S)
    }
    mu + as.vector(crossprod(chol_S, rnorm(nrow(chol_S))))
  } else {
    if (missing(chol_S_inv)) {
      chol_S_inv <- chol(S_inv)
    }
    mu + as.vector(backsolve(chol_S_inv, rnorm(nrow(chol_S_inv)), transpose = TRUE))
  }
}

solve_sym <- function(A, B) {
  R <- chol(A)
  if (missing(B)) {
    chol2inv(R)
  } else {
    chol_solve(R, B)
  }
}

state_catalog <- expand.grid(
  time_index = seq_len(n_months),
  region_index = seq_len(n_land_regions + n_ocean_regions)
) %>%
  mutate(state_index = row_number())

fn <- ncdf4::nc_open(file.path(data_dir, 'transcom/TRANSCOM_mask_GEOS_Chem_4x5.nc'))
base_mask <- ncdf4::ncvar_get(fn, 'Region00')
base_mask[] <- 0
for (region_number in seq_len(n_land_regions + n_ocean_regions)) {
  base_mask[
    ncdf4::ncvar_get(fn, sprintf('Region%02d', region_number)) == 1
  ] <- region_number
}
ncdf4::nc_close(fn)
fn <- ncdf4::nc_open(file.path(data_dir, 'priors/prior_SiB4.nc'))
prior_flux <- ncdf4::ncvar_get(fn, 'NEE')
ncdf4::nc_close(fn)
fn <- ncdf4::nc_open(file.path(data_dir, 'areas/area_5x4.nc'))
area <- ncdf4::ncvar_get(fn, 'grid_cell_area')
ncdf4::nc_close(fn)

region_names <- c(
  'North American Boreal',
  'North American Temperate',
  'Tropical South America',
  'South American Temperate',
  'Northern Africa',
  'Southern Africa',
  'Eurasia Boreal',
  'Eurasia Temperate',
  'Tropical Asia',
  'Australia',
  'Europe',
  'North Pacific Temperate',
  'West Pacific Tropical',
  'East Pacific Tropical',
  'South Pacific Temperate',
  'Northern Ocean',
  'North Atlantic Temperate',
  'Atlantic Tropical',
  'South Atlantic Temperate',
  'Southern Ocean',
  'Indian Tropical',
  'South Indian Temperate'
)

study_times <- seq(ymd_hms('2014-09-01 00:00:00'), ymd_hms('2016-08-01 00:00:00'), by = 'month')

prior_time_region <- data.frame(
  time_index = rep(seq_len(n_months), each = length(base_mask)),
  region_index = rep(as.vector(base_mask), n_months),
  area = rep(as.vector(area), n_months),
  flux = as.vector(prior_flux)
) %>%
  filter(region_index > 0) %>%
  group_by(time_index, region_index) %>%
  summarise(flux = 12 / 44 * 30.5 * 3600 * 24 * 1e3 * 1e-15 * sum(area * flux), .groups = 'drop') %>%
  left_join(
    state_catalog,
    by = c('time_index', 'region_index')
  ) %>%
  mutate(
    time = study_times[time_index],
    region_name = region_names[region_index]
  )

rm(fn)
rm(base_mask)
rm(prior_flux)
rm(area)

region_type_time <- expand.grid(
  time = study_times,
  region_type = factor(c('land', 'ocean'))
) %>%
  mutate(
    region_type_index = row_number()
  )

prior_region_type_time <- prior_time_region %>%
  mutate(
    region_type = factor(ifelse(region_index <= n_land_regions, 'land', 'ocean'))
  ) %>%
  group_by(time, region_type) %>%
  summarise(flux = sum(flux), .groups = 'drop') %>%
  left_join(
    region_type_time,
    by = c('time', 'region_type')
  ) %>%
  arrange(region_type_index)

region_type_time_state <- prior_time_region %>%
  mutate(
    region_type = factor(ifelse(region_index <= n_land_regions, 'land', 'ocean'))
  ) %>%
  group_by(time, region_type, state_index) %>%
  summarise(flux = sum(flux), .groups = 'drop') %>%
  left_join(
    region_type_time,
    by = c('time', 'region_type')
  )

X_state_to_region_month <- matrix(0, nrow = nrow(prior_time_region), ncol = nrow(state_catalog))
X_state_to_region_month[cbind(
  seq_len(nrow(prior_time_region)),
  prior_time_region$state_index
)] <- prior_time_region$flux

X_state_to_region_type_month <- matrix(0, nrow = 2 * n_months, ncol = nrow(state_catalog))
X_state_to_region_type_month[cbind(
  region_type_time_state$region_type_index,
  region_type_time_state$state_index
)] <- region_type_time_state$flux

x_to_land_ocean_flux <- function(x_mean, x_samples) {
  flux_mean <- as.vector(X_state_to_region_type_month %*% (1 + x_mean))
  if (!missing(x_samples)) {
    flux_samples <- X_state_to_region_type_month %*% t(1 + x_samples)
  }

  output <- prior_region_type_time
  output$flux <- flux_mean
  if (!missing(x_samples)) {
    output$flux_q025 <- matrixStats::rowQuantiles(flux_samples, probs = 0.025)
    output$flux_q975 <- matrixStats::rowQuantiles(flux_samples, probs = 0.975)
  }
  output
}

x_to_regional_flux <- function(x_mean, x_samples) {
  flux_mean <- as.vector(X_state_to_region_month %*% (1 + x_mean))
  if (!missing(x_samples)) {
    flux_samples <- X_state_to_region_month %*% t(1 + x_samples)
  }

  output <- prior_time_region
  output$flux <- flux_mean
  if (!missing(x_samples)) {
    output$flux_q025 <- matrixStats::rowQuantiles(flux_samples, probs = 0.025)
    output$flux_q975 <- matrixStats::rowQuantiles(flux_samples, probs = 0.975)
  }
  output
}
