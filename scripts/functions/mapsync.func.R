# Note: this code may run quicker when not packaged into a function for some reason
# Note: parallel processing speeded it up, but stopped working correctly somewhere along the line (weird values)

# for testing
#spcode = "bais"
#startyr = 1994
#endyr = 2019
#n.iter = 3000
#n.post = 3000
#uncertainty = T
#maxzeros = 5
#other.metrics = T
#posterior.data.loc = "data/index_posteriors/"
#jags.data.loc = "data/jags/"
#save.to = "data/mapsync_iterations/" 

mapsync = function(spcode,
                   startyr,
                   endyr,
                   n.iter = 3000,
                   # number of data sets to generate
                   n.post = 3000,
                   # number of posterior values that exist for each grid cell / year
                   uncertainty = T,
                   # logical indicating whether to use the whole posterior or just the median value per grid cell / year
                   maxzeros = 5,
                   other.metrics = F, # logical indicating whether to calculate mean abundance, mean variance, mean trend in addition to mean synchrony
                   posterior.data.loc = "data/index_posteriors/",
                   jags.data.loc = "data/jags/",
                   save.to = "data/mapsync_iterations/") {
  

if(uncertainty == F){n.iter = 1}

# reading in summary of raw data within grid cells including average counts
# this allows excluding based on number of zero counts
zero.check = read.csv(paste0(
  jags.data.loc,
  spcode,
  ".",
  substr(startyr, 3, 4),
  ".",
  substr(endyr, 3, 4),
  ".jags.csv"
)) %>%
  group_by(strat_name, r_year) %>%
  summarise(mean.count = mean(count)) %>% # average route-level counts by grid/year
  ungroup() %>%
  group_by(strat_name) %>%
  summarise(zeros = sum(mean.count == 0)) %>% # count how many "zero years" per grid
  ungroup() %>%
  mutate(exclude = case_when(zeros > maxzeros ~ 1,
                             TRUE ~ 0),
         grid = strat_name) %>%
  dplyr::select(grid, exclude)

# read in posterior distributions for each grid cell / year 
# exclude grid cells with > "maxzeros" number of zero values
post = read.csv(paste0(
  posterior.data.loc,
  spcode,
  ".",
  substr(startyr, 3, 4),
  ".",
  substr(endyr, 3, 4),
  ".index.posteriors.csv"
)) %>%
  rename(grid = Stratum) %>% arrange(grid) %>%
  left_join(zero.check, by = "grid") %>%
  filter(exclude == 0) %>%
  dplyr::select(-exclude)

# split the posterior data into a list of data for each grid cell
post.split = split(post, f = post$grid)

# read number of grids and years
n.grid = length(unique(post$grid)) # number of grid cells for species
n.year = endyr - startyr + 1 # 

# functions to create datasets from random draws of the posteriors
# skipped if uncertainty = F
if(uncertainty == T){
chooserand = function(colnum) {
    lapply(post.split, function(x) x[round(runif(1, 1, n.post)), colnum])
  }

makelots = function(){  
  makeone.list = lapply(4:(n.year+3), chooserand)  # n.year number of years; chooses random column value for each year
  makeone = data.frame(do.call("cbind", makeone.list))
  return(makeone)
}
}

# function to calculate median posterior values if uncertainty == F
if(uncertainty == F){
  getmed = function(x) {apply(x[,4:ncol(x)], 2L, median)}
  makelots = function(){
    makeone.list = lapply(post.split, getmed)
    makeone = data.frame(do.call("cbind", makeone.list))
    return(makeone)
}
}

# creating lots of data sets from random draws of the posteriors
lots.list <- replicate(n.iter, makelots())

if(uncertainty==T){
lots.df = data.frame(do.call("rbind", lots.list))
}

if(uncertainty==F){
lots.df = data.frame(do.call("cbind", lots.list))
}

iteration.list = lapply(1:n.iter, function(x) rep(x, n.year)) 
lots.df$iteration = do.call(c, iteration.list)
lots.df$Year = rep(1:n.year, n.iter)

# split up data.frame to perform transformations on each data set iteration
lots.df.split = split(lots.df, lots.df$iteration)

# a function to linear detrend each time series
lin.func = function(x) {
  lots.c = as.numeric(do.call("rbind", x[1:(ncol(x)-2)])) 
  lots.mat = matrix(lots.c, ncol = (ncol(x)-2), byrow = T) 
  lots.lin = apply(lots.mat, 2L, function(y) resid(lm(y~c(1:26))))
  return(lots.lin)
}

lin.list <- lapply(lots.df.split, lin.func)

if(other.metrics == T) {
  # a function to take the mean of each time series (for when "other.metrics" = T)
  mean.func = function(x) {
    lots.c = as.numeric(do.call("rbind", x[1:(ncol(x) - 2)]))
    lots.mat = matrix(lots.c, ncol = (ncol(x) - 2), byrow = T)
    lots.mean = colMeans(lots.mat, 1L)
    return(matrix(lots.mean, 1, 15))
  }
  
  mean.list <- lapply(lots.df.split, mean.func)
  
  # a function to calculate variance of each detrended time series (for when "other.metrics" = T)
  var.func = function(x) {
    lots.var = apply(x, 2L, var)
    return(matrix(lots.var, 1, 15))
  }
  
  var.list <- lapply(lin.list, var.func)
  
  # a function to extract the linear trend coefficient from each time series
  trend.func = function(x) {
    lots.c = as.numeric(do.call("rbind", x[1:(ncol(x) - 2)]))
    lots.mat = matrix(lots.c, ncol = (ncol(x) - 2), byrow = T)
    lots.trend = apply(lots.mat, 2L, function(y)
      coef(lm(y ~ c(1:26)))[2])
    return(matrix(lots.trend, 1, 15))
  }
  
  trend.list <- lapply(lots.df.split, trend.func)
}

# create correlation matrix among all grid cells (of linear detrended time series)
lin.cor = lapply(lin.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
lin.cor.n = lapply(lin.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$n)
if(min(as.vector(do.call("rbind",lin.cor.n)))<20){warning("At least one grid pair had n < 20 for correlations")}

# combine all iterations of the linear detrended correlation matrices
lin.cor.mat = do.call("rbind", lin.cor)

if(other.metrics == T){
  mean.df = as.data.frame(do.call("rbind", mean.list)) %>%
    mutate(iteration = 1:n.iter) %>%
    pivot_longer(cols = 1:15) %>%
    mutate(name = rep(unique(post$grid),n.iter)) %>%
    rename(mean = value, grid = name)

  var.df = as.data.frame(do.call("rbind", var.list)) %>%
    mutate(iteration = 1:n.iter) %>%
    pivot_longer(cols = 1:15) %>%
    mutate(name = rep(unique(post$grid),n.iter)) %>%
    rename(var = value, grid = name)

  trend.df = as.data.frame(do.call("rbind", trend.list)) %>%
    mutate(iteration = 1:n.iter) %>%
    pivot_longer(cols = 1:15) %>%
    mutate(name = rep(unique(post$grid),n.iter)) %>%
    rename(trend = value, grid = name)
}

# combine all matrices into a "long" data.frame with grid labels
gridA = rep(unique(post$grid), n.grid)
gridB.list = lapply(1:n.grid, function(x)
  rep(unique(post$grid)[x], n.grid))
gridB = c(do.call("cbind", gridB.list))
iter.list = lapply(1:n.iter, function(x) rep(x, n.grid*n.grid))
grid.cor = cbind.data.frame(
  gridA = rep(gridA, n.iter),
  gridB = rep(gridB, n.iter),
  lin.cor = as.vector(t(lin.cor.mat))
) %>%
  mutate(upper.tri = rep(as.vector(t(
    upper.tri(lin.cor.mat[1:n.grid, 1:n.grid], diag = TRUE)
  )), n.iter)) %>%
  mutate(iteration = c(do.call("cbind", iter.list))) %>%
  filter(upper.tri == T) %>% # this excludes duplicate and diagonal correlations from the matrix
  dplyr::select(-upper.tri) %>%
  arrange(iteration, gridA, gridB) 

# calculate inter-grid distances for 2 X 2 degree lat/long grid cells
gridlat = as.numeric(substr(as.character(unique(post$grid)), 1, 2)) - 0.5
gridlong = as.numeric(substr(as.character(unique(post$grid)), 3, 10)) + 0.5
grid.spatial = sp::SpatialPointsDataFrame(
  cbind(gridlong, gridlat),
  cbind.data.frame(grid = as.character(unique(post$grid)), as.numeric(gridlong), as.numeric(gridlat))
)
sp::proj4string(grid.spatial) = sp::CRS("+init=epsg:4326")

df.grid.spatial = data.frame(
  grid = grid.spatial$grid,
  lat = sp::coordinates(grid.spatial)[, 2],
  lon = sp::coordinates(grid.spatial)[, 1]
)

#calculate Haversine distances between grid cells
dmat.grid = round(geosphere::distm(df.grid.spatial[, c("lon", "lat")]) /
                    1000)

dgrid = data.frame(
  gridA = gridB,
  gridB = gridA,
  dist = as.vector(dmat.grid),
  upper.tri = as.vector(upper.tri(dmat.grid, diag = TRUE))
) %>%
  filter(upper.tri == T) %>%
  dplyr::select(-upper.tri) %>%
  arrange(gridA, gridB)
  
# adding in the distance between grid cells again
grid.cordist = grid.cor %>%
  left_join(dgrid, by = c('gridA', 'gridB')) %>%
  dplyr::select(gridA,
                gridB,
                lin.cor,
                dist,
                iteration) %>%
  arrange(iteration, gridA, gridB) %>%
  mutate(lin.cor = case_when(dist != 0 ~ lin.cor,
                             TRUE ~ as.numeric(NA)))

if (other.metrics == T) {
  grid.cordist = grid.cordist %>%
    left_join(mean.df, by = c("gridA" = "grid", "iteration")) %>%
    rename(meanA = mean) %>%
    left_join(mean.df, by = c("gridB" = "grid", "iteration")) %>%
    rename(meanB = mean) %>%
    left_join(var.df, by = c("gridA" = "grid", "iteration")) %>%
    rename(varA = var) %>%
    left_join(var.df, by = c("gridB" = "grid", "iteration")) %>%
    rename(varB = var) %>%
    left_join(trend.df, by = c("gridA" = "grid", "iteration")) %>%
    rename(trendA = trend) %>%
    left_join(trend.df, by = c("gridB" = "grid", "iteration")) %>%
    rename(trendB = trend)
}

# create sets of adjacent cells within 400 km ("grid" being the focal cell)
mapsync_step1 =
  rbind(
    data.frame(grid.cordist, grid = grid.cordist$gridA),
    data.frame(
      filter(grid.cordist, dist != 0),
      grid = filter(grid.cordist, dist != 0)$gridB
    )
  ) %>%
  filter(dist < 400)

if (other.metrics == T) {
  mapsync_step1 = mapsync_step1 %>%
    mutate(
      mean = case_when(grid == gridA ~ meanB,
                       TRUE ~ meanA),
      var = case_when(grid == gridA ~ varB,
                      TRUE ~ varA),
      trend = case_when(grid == gridA ~ trendB,
                        TRUE ~ trendA)
    ) %>%
    dplyr::select(-meanA, -meanB, -varA, -varB, -trendA, -trendB)
}

# create 400 km moving window average of synchrony, abundance, variance, trend
if (other.metrics == T) {
  mapsync = mapsync_step1 %>%
    group_by(iteration, grid) %>%
    summarise(
      meancor.lin = mean(lin.cor, na.rm = T),
      mean.mean = mean(mean),
      mean.var = mean(var),
      mean.trend = mean(trend),
      n.grids.per.mean = length(grid) # for mean correlations (others include center cell, so = + 1)
    ) %>%
    ungroup() %>%
    mutate(grid = as.character(grid)) %>%
    filter(n.grids.per.mean > 1) %>%
    dplyr::select(-n.grids.per.mean) %>%
    arrange(iteration, grid)
} else{
  mapsync = mapsync_step1 %>%
    group_by(iteration, grid) %>%
    summarise(meancor.lin = mean(lin.cor, na.rm = T),
              n.grids.per.mean = length(grid)) %>% # for mean correlations (others include center cell, so = + 1)) %>%
    ungroup() %>%
    mutate(grid = as.character(grid)) %>%
    filter(n.grids.per.mean > 1) %>%
    dplyr::select(-n.grids.per.mean) %>%
    arrange(iteration, grid)
}

# write csv for iterations method (uncertainty == T)
if(uncertainty==T){
write.csv(mapsync, paste0(save.to, spcode,".",substr(startyr,3,4),".", substr(endyr,3,4),".mapsync.iterations.csv"), row.names=FALSE)
}

# write csv for "median posterior" method (uncertainty == F)
if(uncertainty==F){
write.csv(mapsync, paste0(save.to, spcode,".",substr(startyr,3,4),".", substr(endyr,3,4),".mapsync.medians.csv"), row.names=FALSE)
}
}
