# Note: this code may run quicker when not packaged into a function for some reason
# Note: parallel processing speeded it up, but stopped working correctly somewhere along the line (weird values)

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
                   jags.data.loc = "data/raw_counts/",
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
  summarise(mean.count = mean(count)) %>%
  ungroup() %>%
  group_by(strat_name) %>%
  summarise(zeros = sum(mean.count == 0)) %>%
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
  lots.c = as.numeric(do.call("rbind", x[1:(ncol(x)-2)])) # n.grid = number of grid cells
  lots.mat = matrix(lots.c, ncol = (ncol(x)-2), byrow = T) # n.grid = number of grid cells
  lots.lin = apply(lots.mat, 2L, function(y) resid(lm(y~c(1:26))))
  return(lots.lin)
}

lin.list <- lapply(lots.df.split, lin.func)

# create correlation matrix among all grid cells (linear detrend, and ar1 and diff methods if called)
lin.cor = lapply(lin.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
lin.cor.n = lapply(lin.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$n)
if(min(as.vector(do.call("rbind",lin.cor.n)))<20){warning("At least one grid pair had n < 20 for correlations")}
if(all.methods == T){
lin.l.cor = lapply(lin.l.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
ar1.cor <- lapply(ar1.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
ar1.l.cor <- lapply(ar1.l.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
dif.cor = lapply(dif.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
dif.l.cor = lapply(dif.l.list, function(x) psych::corr.test(x, adjust = "none", ci = F)$r)
}

# combine all iterations of the linear and (if called for) AR1 and differenced detrended correlation matrices
lin.cor.mat = do.call("rbind", lin.cor)

# combine all matrices into a "long" data.frame with grid labels
gridA = rep(unique(post$grid), n.grid)
gridB.list = lapply(1:n.grid, function(x)
  rep(unique(post$grid)[x], n.grid))
gridB = c(do.call("cbind", gridB.list))
iter.list = lapply(1:n.iter, function(x) rep(x, n.grid*n.grid))
if(all.methods == F){
grid.cor = cbind.data.frame(
  gridA = rep(gridA, n.iter),
  gridB = rep(gridB, n.iter),
  lin.cor = as.vector(t(lin.cor.mat))
) %>%
  mutate(upper.tri = rep(as.vector(t(
    upper.tri(lin.cor.mat[1:n.grid, 1:n.grid])
  )), n.iter)) %>%
  mutate(iteration = c(do.call("cbind", iter.list))) %>%
  filter(upper.tri == T) %>% # this excludes duplicate and diagonal correlations from the matrix
  dplyr::select(-upper.tri) %>%
  arrange(iteration, gridA, gridB)
}

if (all.methods == T) {
  grid.cor = cbind.data.frame(
    gridA = rep(gridA, n.iter),
    gridB = rep(gridB, n.iter),
    lin.cor = as.vector(t(lin.cor.mat)),
    lin.l.cor = as.vector(t(lin.l.cor.mat)),
    ar1.cor = as.vector(t(ar1.cor.mat)),
    ar1.l.cor = as.vector(t(ar1.l.cor.mat)),
    dif.cor = as.vector(t(dif.cor.mat)),
    dif.l.cor = as.vector(t(dif.l.cor.mat))
  ) %>%
    mutate(upper.tri = rep(as.vector(t(
      upper.tri(lin.cor.mat[1:n.grid, 1:n.grid])
    )), n.iter)) %>%
    mutate(iteration = c(do.call("cbind", iter.list))) %>%
    filter(upper.tri == T) %>% # this excludes duplicate and diagonal correlations from the matrix
    dplyr::select(-upper.tri) %>%
    arrange(iteration, gridA, gridB)
}

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
  upper.tri = as.vector(upper.tri(dmat.grid))
) %>%
  filter(upper.tri == T) %>%
  dplyr::select(-upper.tri) %>%
  arrange(gridA, gridB)
  
# adding in the distance between grid cells again
if(all.methods == F){
grid.cordist = grid.cor %>%
  left_join(dgrid, by = c('gridA', 'gridB')) %>%
  dplyr::select(gridA,
                gridB,
                lin.cor,
                dist,
                iteration) %>%
  arrange(iteration, gridA, gridB)
}

if (all.methods == T) {
  grid.cordist = grid.cor %>%
    left_join(dgrid, by = c('gridA', 'gridB')) %>%
    dplyr::select(
      gridA,
      gridB,
      lin.cor,
      lin.l.cor,
      ar1.cor,
      ar1.l.cor,
      dif.cor,
      dif.l.cor,
      dist,
      iteration
    ) %>%
    arrange(iteration, gridA, gridB)
}

# create moving average of all adjacent cells within 400 km
if(all.methods == F){
mapsync =
  rbind(
    data.frame(grid.cordist, grid = grid.cordist$gridA),
    data.frame(grid.cordist, grid = grid.cordist$gridB)
  ) %>%
  filter(dist < 400) %>%
  group_by(iteration, grid) %>%
  summarise(
    meancor.lin = mean(lin.cor),
    n.grids.per.mean = length(grid)
  ) %>%
  ungroup() %>%
  mutate(grid = as.character(grid)) %>%
  filter(n.grids.per.mean > 1) %>%
  dplyr::select(-n.grids.per.mean) %>%
  arrange(iteration, grid)
}

if(all.methods == T){
mapsync =
  rbind(
    data.frame(grid.cordist, grid = grid.cordist$gridA),
    data.frame(grid.cordist, grid = grid.cordist$gridB)
  ) %>%
  filter(dist < 400) %>%
  group_by(iteration, grid) %>%
  summarise(
    meancor.lin = mean(lin.cor),
    meancor.lin.l = mean(lin.l.cor),
    meancor.ar1 = mean(ar1.cor),
    meancor.ar1.l = mean(ar1.l.cor),
    meancor.dif = mean(dif.cor),
    meancor.dif.l = mean(dif.l.cor),
    n.grids.per.mean = length(grid)
  ) %>%
  ungroup() %>%
  mutate(grid = as.character(grid)) %>%
  filter(n.grids.per.mean > 1) %>%
  dplyr::select(-n.grids.per.mean) %>%
  arrange(iteration, grid)
}

# write csv for iterations method (uncertainty == T)
if(uncertainty==T){
write.csv(mapsync, paste0(save.to, spcode,".",substr(startyr,3,4),".", substr(endyr,3,4),".mapsync.iterations.csv"))
}

# write csv for "median posterior" method (uncertainty == F)
if(uncertainty==F){
write.csv(mapsync, paste0(save.to, spcode,".",substr(startyr,3,4),".", substr(endyr,3,4),".mapsync.medians.csv"))
}
}