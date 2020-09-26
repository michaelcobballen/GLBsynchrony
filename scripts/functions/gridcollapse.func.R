#####################################
#### Function to prepare BBS data for Jags model using bbsBayes package
#####################################
gridcollapse = function(species, spcode, startyr, endyr, save.loc = "data/jags/"){
  
# stratify BBS data into 1 x 1 latitude / longitude grid cells

  stratified_data <- stratify(by = "latlong")
  
# re-classify data into 2 x 2 latitude / longitude grid cells

  # functions for identifying even and odd numbers
  is.even <- function(x) x %% 2 == 0
  is.odd <- function(x) x %% 2 != 0
  
  # pulling out the route_strat data and converting latlong into 2x2 grid ID
  rs = stratified_data$route_strat %>%
    mutate(grid_lat = ifelse(is.odd(trunc(Latitude)),trunc(Latitude)+1,trunc(Latitude))) %>%
    mutate(grid_lon = ifelse(is.odd(trunc(Longitude)),trunc(Longitude)-1,trunc(Longitude))) %>%
    mutate(strat_name = paste(grid_lat,grid_lon,sep="")) %>%
    dplyr::select(-grid_lat, -grid_lon)
  
  stratified_data[[2]] <- rs
  

# create csv jags data for a species

  # prepare data for Jags model
  jags_data <- prepare_jags_data(stratified_data, 
                                 species_to_run = species,
                                 min_n_routes = 2,
                                 min_max_route_years = 10,
                                 model = "gamye",
                                 heavy_tailed = T,
                                 min_year = startyr,
                                 max_year = endyr)
  
# collect relevant parts of jags_data into data frame

  data = data.frame(
    ncounts = c(jags_data[["ncounts"]], rep(NA,jags_data[["ncounts"]]-1)),
    nstrata = c(jags_data[["nstrata"]], rep(NA,jags_data[["ncounts"]]-1)),
    ymin = c(jags_data[["ymin"]], rep(NA,jags_data[["ncounts"]]-1)),
    ymax = c(jags_data[["ymax"]], rep(NA,jags_data[["ncounts"]]-1)),
    nknots = c(jags_data[["nknots"]], rep(NA,jags_data[["ncounts"]]-1)),
    stratify_by = c(jags_data[["stratify_by"]], rep(NA,jags_data[["ncounts"]]-1)),
    nonzeroweight = c(jags_data[["nonzeroweight"]],rep(NA,jags_data[["ncounts"]]-jags_data[["nstrata"]])),
    nobservers = c(jags_data[["nobservers"]], rep(NA,jags_data[["ncounts"]]-jags_data[["nstrata"]])),
    X.basis = c(jags_data[["X.basis"]], rep(NA, jags_data[["ncounts"]]-length(jags_data[["X.basis"]]))),
    count = jags_data[["count"]],
    strat_name = jags_data[["strat_name"]],
    strat = jags_data[["strat"]],
    obser = jags_data[["obser"]],
    year = jags_data[["year"]],
    firstyr = jags_data[["firstyr"]],
    month = jags_data[["month"]],
    day = jags_data[["day"]],
    r_year = jags_data[["r_year"]], 
    route = jags_data[["route"]]
  )
  
# export to csv for use in Amarel cluster script

  write.csv(data, paste0(save.loc, spcode, ".", substr(startyr,3,4), ".", substr(endyr,3,4), ".jags.csv"))
}

