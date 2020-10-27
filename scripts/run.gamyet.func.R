# Function to run the jagsUI model and jags_data exported from bbsBayes package
# (1x1 degree grid cell strata were first merged into 2x2 cells)
# (then jags_data was packaged into a csv to allow use of super computing cluster)
# (this was done as bbsBayes package wouldn't install on the cluster)

run.gamyet = function(spcode,
                      species,
                      startyr,
                      endyr,
                      chains = 3,
                      burn = 20000,
                      iter = 30000,
                      thin = 10,
                      par = T,
                      gamyet.loc = "scripts/gamyet.bug",
                      data.loc = "data/jags/",
                      save.posteriors = "data/index_posteriors/",
                      warnings = F) {
  
#####################################
#### # load libraries
#####################################
library(lattice)
library(jagsUI, warn.conflicts = FALSE)
library(dplyr)

#####################################
#### read in data and make into list form for jags
#####################################
data = read.csv(paste0(data.loc, spcode, ".", substr(startyr,3,4), ".", substr(endyr,3,4), ".jags.csv"))

jags_data = list(
  ncounts = data$ncounts[1],
  nstrata = data$nstrata[1],
  ymin = data$ymin[1],
  ymax = data$ymax[1],
  nonzeroweight = data$nonzeroweight[1:data$nstrata[1]],
  count = data$count,
  strat = data$strat,
  obser = data$obser,
  year = data$year,
  firstyr = data$firstyr,
  nobservers = data$nobservers[1:data$nstrata[1]],
  month = data$month,
  day = data$day,
  nknots = data$nknots[1],
  X.basis = matrix(c(data[is.na(data$X.basis)==F,]$X.basis), nrow = (data$ymax[1]-data$ymin[1]+1))
)

# note: intentionally left things out of the jags_data list that were made NULL in run-model.R

#####################################
#### Run jags model
#####################################
jags_job <- jagsUI::jags(data = jags_data,
                         inits = NULL,
                         parameters.to.save = c("n"),
                         model.file = gamyet.loc, # file path to model
                         n.chains = chains,
                         n.adapt = NULL,
                         n.burnin = burn,
                         n.iter = iter,
                         n.thin = thin,
                         parallel = T,
                         verbose = TRUE,
                         modules = NULL)

jags_job$strat_name <- data$strat_name
jags_job$stratify_by <- data$stratify_by[1]
jags_job$r_year <- data$r_year


#####################################
#### load the r_hat function
#####################################
r_hat <- function(jags_mod = NULL,
                  parameter_list = NULL,
                  threshold = NULL)
{
  rhat_list <- jags_mod$Rhat
  if (is.null(parameter_list))
  {
    parameter_list <- names(rhat_list)
  }
  to_return <- data.frame(Parameter = character(),
                          R_hat = double(),
                          stringsAsFactors = FALSE)
  
  for (p in parameter_list)
  {
    if (isFALSE(p %in% names(rhat_list)))
    {
      warning(paste("Parameter", p, "not found in supplied model."))
      next
    }
    current <- rhat_list[[p]]
    l <- length(dim(current))
    if (l == 0)
    {
      to_return <- rbind(to_return,
                         data.frame(Parameter = p,
                                    R_hat = current))
    }else if (l == 1)
    {
      for (i in 1:dim(current))
      {
        to_return <- rbind(to_return,
                           data.frame(Parameter = paste0(p, "[", i, "]"),
                                      R_hat = current[i]))
      }
    }else if (l == 2)
    {
      for (i in 1:dim(current)[1])
      {
        for (j in 1:dim(current)[2])
        {
          to_return <- rbind(to_return,
                             data.frame(Parameter = paste0(p, "[", i, "][", j, "]"),
                                        R_hat = current[i,j]))
        }
      }
    }else if (l == 3)
    {
      for (i in 1:dim(current)[1])
      {
        for (j in 1:dim(current)[2])
        {
          for (k in 1:dim(current)[3])
          {
            to_return <- rbind(to_return,
                               data.frame(Parameter = paste0(p, "[", i, "][", j, "][", k, "]"),
                                          R_hat = current[i,j,k]))
          }
        }
      }
    }
  }
  
  if (!is.null(threshold))
  {
    to_return <- to_return[which(to_return$R_hat >= threshold), ]
  }
  
  return(to_return)
  
}

#####################################
#### check the Rhat values for convergence failure
#### if there are failures, then throw a warning
#### alternatively, incoporate a logical option to automatically continue running the model until some minimum convergence criterion is met
#####################################
warning1 = "none"
warning2 = "none"

rhat_check = r_hat(jags_job,
                   threshold = 1.1)
if(nrow(rhat_check) > 1){
  failed = paste(rhat_check$Parameter,collapse = " ; ")
  nfail = nrow(rhat_check)
  warning1 = (paste("Warning",nfail,"parameters did not converge (Rhat > 1.1). Consider re-running with a longer burn-in and-or more posterior samples."))
  
  warning2 = (paste("Convergence failure on the following parameters:",failed))
}

#####################################
#### extract-index_data() function guts
#####################################

jags_mod = jags_job
alt_n = "n"

# functions for identifying even and odd numbers
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

  # Read area weights
  # (edited to collapse latlong grid cells into 2x2 format)
  stratify_by <- jags_mod$stratify_by
  all_area_weights <- utils::read.csv(paste0(data.loc,"db.csv")) %>%
    mutate(nchar = nchar(region)) %>%
    mutate(Latitude = as.numeric(substr(region, 1, 2))) %>%
    mutate(Longitude = as.numeric(case_when(nchar == 6 ~ substr(region, 4, 6),
                                            nchar == 7 ~ substr(region, 4, 7)))) %>%
    mutate(grid_lat = ifelse(is.odd(trunc(Latitude)),trunc(Latitude)+1,trunc(Latitude))) %>%
    mutate(grid_lon = ifelse(is.odd(trunc(Longitude)),trunc(Longitude)-1,trunc(Longitude))) %>%
    mutate(strat_name = paste(grid_lat,grid_lon,sep="")) %>%
    dplyr::select(-grid_lat, -grid_lon, -nchar, -Latitude, -Longitude) %>%
    group_by(strat_name) %>%
    summarize(area_sq_km = sum(area_sq_km)) %>%
    ungroup() %>%
    rename(region = strat_name)
  
  # Extract posterior data and other data from jags_mod
  n <- jags_mod$sims.list[["n"]]
  if (isTRUE(jags_mod$parallel))
  {
    bugs_data <- jags_mod$model$cluster1$data()
  }else
  {
    bugs_data <- jags_mod$model$data()
  }
  y_min <- bugs_data$ymin
  y_max <- bugs_data$ymax
  strat_list = unique(data.frame(strat_name = jags_mod$strat_name,
                                 strat = bugs_data$strat,
                                 stringsAsFactors = FALSE))
  strat_list = strat_list[order(strat_list$strat),]
  
  # Subset area weights based on strata used and ensure same order as JAGS
  strata_used <- strat_list$strat_name
  strata_num <- strat_list$strat
  area_weights <- all_area_weights[which(all_area_weights$region %in% strata_used), ]
  area_weights <- area_weights[ order(match(area_weights$region, strata_used)),]
  area_weights$num <- strata_num
  
  #####################
  #### Inserting get_prepared_data() function code from R script in bbsBayes package
  #####################
  get_prepared_data <- function(jags_data = NULL)
  {
    to_return <- data.frame(Year = data$r_year,
                            Year_Factored = jags_data$year,
                            Count = jags_data$count,
                            Stratum = data$strat_name,
                            Stratum_Factored = jags_data$strat,
                            Observer_Factored = jags_data$obser,
                            Route = data$route,
                            First_Year = jags_data$firstyr)
    
    return(to_return)
  }
  #####################
  
  if(!is.null(jags_data)){
    original_data = get_prepared_data(jags_data = jags_data)
  }else{
    original_data = NULL
  }
  
  data_list = list(n = n,
              area_weights = area_weights,
              y_min = y_min,
              y_max = y_max,
              r_year = jags_mod$r_year,
              bugs_data = bugs_data,
              original_data = original_data)

  
  # make list of grid cell names in the original order (c/o Adam Smith)
  df_dat <- get_prepared_data(jags_data)
  strat_table <- unique(df_dat[, c("Stratum", "Stratum_Factored")])
  strat_table <- strat_table[order(strat_table$Stratum_Factored),]
  
  # make a list of all grid cells (to match order of final index)
  npost = ((iter-burn)*3)/10
  gridlist = c()
  for (i in 1:dim(data_list$n)[2]) {
    gridlist[[i]] = rep(strat_table$Stratum[i], npost)
  }

  # format the posterior outputs
  dimnames(data_list$n) <- list(1:npost, strat_table$Stratum, startyr:endyr) 
  n = apply(data_list$n, 3L, c)
  
  # create summary of mean raw count data
  raw.data.warnings = data_list$original_data %>%
    group_by(Stratum, Year) %>%
    summarize(mean.count = mean(Count), n.routes = length(Count)) %>%
    ungroup() %>%
    arrange(Stratum, Year) %>%
    mutate(warning1 = warning1,
           warning2 = warning2)

# create final model output file (posterior distributions for each cell/year combo)
model.output = cbind.data.frame(
Stratum = c(t(do.call("rbind", gridlist))),
n.draw = rep(1:npost, dim(data_list$n)[2]),
n
) 

#####################################
#### save results to csv files
#####################################
  
# write index posterior draws to csv file
write.csv(
  model.output,
  paste0(
    save.posteriors,
    spcode,
    ".",
    substr(startyr, 3, 4),
    ".",
    substr(endyr, 3, 4),
    ".index.posteriors.csv"
  ),
  row.names = FALSE
)

# if "warnings == T", write csv of mean raw counts per grid & convergence warnings

if (warnings == T) {
  write.csv(
    raw.data.warnings,
    paste0(
      save.posteriors,
      spcode,
      ".",
      substr(startyr, 3, 4),
      ".",
      substr(endyr, 3, 4),
      ".rawcount.warnings.csv"
    ),
    row.names = FALSE
  )
}

}
