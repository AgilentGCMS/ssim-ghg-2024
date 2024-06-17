### Load the Jacobian

load(file.path(data_dir, 'jacobians', 'trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda'))
H_all <- jacob * 12 / 44
rm(jacob)

### Load the data

load(file.path(data_dir, 'obs', 'obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda'))
obs_catalog$observation_mode <- with(obs_catalog, case_when(
  TYPE == 'IS' ~ 'IS',
  TYPE == 'TCCON' ~ 'TCCON',
  TYPE == 'OCO2' ~ case_when(
    substring(ID, 14) == '1' ~ 'LN',
    substring(ID, 14) == '2' ~ 'LG',
    substring(ID, 14) == '3' ~ 'LT',
    substring(ID, 14) == '4' ~ 'LTT',
    substring(ID, 14) == '6' ~ 'OG'
  )
))

### Subset the matrices

obs_to_include <- with(
  obs_catalog,
  observation_mode %in% c('IS', 'LN', 'LG')
    & DATE >= ymd_hms('2014-09-01 00:00:00')
    & DATE < ymd_hms('2015-09-01 00:00:00')
)
H <- H_all[obs_to_include, ]

obs_catalog <- obs_catalog[obs_to_include, ]

rm(H_all)
rm(obs_to_include)
