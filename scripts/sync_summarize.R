
# make a list of files containing the mean synchrony data sets
mapsync.files = sprintf(
  "data//mapsync_iterations//%s.mapsync.iterations.csv",
  paste0(
    species$spcode,
    ".",
    substr(species$startyr, 3, 4),
    ".",
    substr(species$endyr, 3, 4)
  )
)

# read in all mapsync data into list form (takes a few minutes)
mapsync.iterations = lapply(mapsync.files, function(x)
  read.csv(x))
names(mapsync.iterations) <-   paste0(species$spcode,
                                      ".",
                                      substr(species$startyr, 3, 4),
                                      ".",
                                      substr(species$endyr, 3, 4))

# get a list of the grid cells in common between periods for each species (for filtering)
numsp = length(unique(species$spcode))
gridcom = function(x) {
  Reduce(intersect, list(
    unique(mapsync.iterations[[x]]$grid),
    unique(mapsync.iterations[[x + numsp]]$grid)
  ))
}

grid.common = lapply(1:numsp, gridcom)
names(grid.common) = unique(species$spcode)

grid.common2 = c(grid.common, grid.common)

# subset only the grid cells in common between the two time periods
mapsync.com = list()
for (i in 1:nrow(species)) {
  mapsync.com[[i]] = filter(mapsync.iterations[[i]], grid %in% grid.common2[[i]])
}
names(mapsync.com) <- species$spcode

# separate out period one, period two, and subtract the two to get change in synchrony
mapsync.com1 = mapsync.com[1:numsp]
mapsync.com2 = mapsync.com[(numsp+1):nrow(species)]
mapsync.dif1 = Map(function(listper2, listper1)
  listper2[, 4:ncol(listper2)] - listper1[, 4:ncol(listper2)],
  mapsync.com2,
  mapsync.com1)

# add in iteration and grid variables back to synchrony change
mapsync.dif = Map(function(mainlist, listnames)
  cbind(listnames[, 2:3], mainlist[, 1:ncol(mainlist)]),
  mapsync.dif1,
  mapsync.com1)

# add in species names
mapsync.dif.names = list()
for (i in 1:numsp) {
  mapsync.dif.names[[i]] = cbind.data.frame(species = species$spcode[i], mapsync.dif[[i]])
}

# create giant data frame with all iterations of dif, per1, and per2 for all spp
syncdif.all = cbind(
  do.call("rbind", mapsync.dif.names),
  # difference
  do.call("rbind", mapsync.com1)[, c(4:ncol(do.call("rbind", mapsync.com1)))],
  # period 1
  do.call("rbind", mapsync.com2)[, c(4:ncol(do.call("rbind", mapsync.com2)))] # period 2
)

# Name columns appropriately for the number of mapsync data columns (i.e., made using "all.methods = T" or "F")
if(ncol(mapsync.dif[[1]]) == 8){
colnames(syncdif.all) = c(
  "spcode",
  "iteration",
  "grid",
  "ar1.dif",
  "lin.dif",
  "diff.dif",
  "ar1.l.dif",
  "lin.l.dif",
  "diff.l.dif",
  "ar1.one",
  "lin.one",
  "diff.one",
  "ar1.l.one",
  "lin.l.one",
  "diff.l.one",
  "ar1.two",
  "lin.two",
  "diff.two",
  "ar1.l.two",
  "lin.l.two",
  "diff.l.two")} else{
    colnames(syncdif.all) = c(
      "spcode",
      "iteration",
      "grid",
      "lin.dif",
      "lin.one",
      "lin.two")
}

syncdif = syncdif.all %>%
  group_by(iteration, grid) %>%
  mutate(spp.per.grid = length(spcode)) %>%
  summarise_at(vars(-spcode), mean
  ) %>%
  ungroup() %>%
  filter(spp.per.grid >= 3) %>%
  dplyr::select(-spp.per.grid)

# summarize synchrony change iterations by quantiles (all species)
syncdif.split = split(syncdif, f = syncdif$grid)

syncdif.sum.list = lapply(syncdif.split, function(x)
  apply(x[, 3:ncol(syncdif)], 2L, FUN = quantile, probs = c(0.05, 0.5, 0.95)))

syncdif.sum = data.frame(do.call("rbind", syncdif.sum.list)) %>%
  mutate(quant = rep(c("p5", "med", "p95"), length(unique(syncdif$grid))),
         grid = c(mapply(
           unique(syncdif$grid),
           FUN = function(x)
             rep(x, 3)
         ))) %>%
  mutate(
    lat = as.numeric(substr(grid, 1, 2)) - 0.5,
    lon = as.numeric(substr(grid, 3, 10)) + 0.5
  )

# join averaged species syncdif data to individual species data
syncdif.spp.prep = syncdif.all %>%
  bind_rows(syncdif) %>%
  mutate(spcode = as.factor(case_when(is.na(spcode) ~ "All species",
                                      TRUE ~ spcode)))

# synchrony change of each grid cell (median of 3000 iterations) for each species
syncdif.spp.map = cbind(aggregate(
  syncdif.all[,c(4:ncol(syncdif.all))],
  by = list(syncdif.all$spcode, syncdif.all$grid),
  median
),
aggregate(
  syncdif.all[,c(4:ncol(syncdif.all))],
  by = list(syncdif.all$spcode, syncdif.all$grid),
  FUN = function(x) quantile(x, 0.05)
),
aggregate(
  syncdif.all[,c(4:ncol(syncdif.all))],
  by = list(syncdif.all$spcode, syncdif.all$grid),
  FUN = function(x) quantile(x, 0.95)
)
)
if(ncol(mapsync.dif[[1]]) == 8) {
  colnames(syncdif.spp.map) = c(
    "spcode",
    "grid",
    "ar1.dif_med",
    "lin.dif_med",
    "diff.dif_med",
    "ar1.l.dif_med",
    "lin.l.dif_med",
    "diff.l.dif_med",
    "ar1.one_med",
    "lin.one_med",
    "diff.one_med",
    "ar1.l.one_med",
    "lin.l.one_med",
    "diff.l.one_med",
    "ar1.two_med",
    "lin.two_med",
    "diff.two_med",
    "ar1.l.two_med",
    "lin.l.two_med",
    "diff.l.two_med",
    "w",
    "x",
    "ar1.dif_p5",
    "lin.dif_p5",
    "diff.dif_p5",
    "ar1.l.dif_p5",
    "lin.l.dif_p5",
    "diff.l.dif_p5",
    "ar1.one_p5",
    "lin.one_p5",
    "diff.one_p5",
    "ar1.l.one_p5",
    "lin.l.one_p5",
    "diff.l.one_p5",
    "ar1.two_p5",
    "lin.two_p5",
    "diff.two_p5",
    "ar1.l.two_p5",
    "lin.l.two_p5",
    "diff.l.two_p5",
    "y",
    "z",
    "ar1.dif_p95",
    "lin.dif_p95",
    "diff.dif_p95",
    "ar1.l.dif_p95",
    "lin.l.dif_p95",
    "diff.l.dif_p95",
    "ar1.one_p95",
    "lin.one_p95",
    "diff.one_p95",
    "ar1.l.one_p95",
    "lin.l.one_p95",
    "diff.l.one_p95",
    "ar1.two_p95",
    "lin.two_p95",
    "diff.two_p95",
    "ar1.l.two_p95",
    "lin.l.two_p95",
    "diff.l.two_p95"
  )
} else{
  colnames(syncdif.spp.map) = c(
    "spcode",
    "grid",
    "lin.dif_med",
    "lin.one_med",
    "lin.two_med",
    "w",
    "x",
    "lin.dif_p5",
    "lin.one_p5",
    "lin.two_p5",
    "y",
    "z",
    "lin.dif_p95",
    "lin.one_p95",
    "lin.two_p95"
  )
}
syncdif.spp.map$lat = as.numeric(substr(syncdif.spp.map$grid,1,2))-0.5
syncdif.spp.map$lon = as.numeric(substr(syncdif.spp.map$grid,3,10))+0.5
syncdif.spp.map = syncdif.spp.map %>%
  dplyr::select(-w, -x, -y, -z) %>%
  arrange(spcode, grid)

# average synchrony change of grid cells within species and iterations
syncdif.spp = cbind(aggregate(
  syncdif.spp.prep[,c(4:ncol(syncdif.spp.prep))],
  by = list(syncdif.spp.prep$iteration, syncdif.spp.prep$spcode),
  mean
),
aggregate(
  syncdif.spp.prep[,c("lin.dif")],
  by = list(syncdif.spp.prep$iteration, syncdif.spp.prep$spcode),
  FUN = function(x) length(x)
)
)

# get the median and quantiles of the mean synchrony change iterations 
# for plotting, individual species and all species averaged
syncdif.spp.sum = cbind(aggregate(
  syncdif.spp[,c(3:(ncol(syncdif.spp)-3))],
  by = list(syncdif.spp$Group.2),
  FUN = function(x) quantile(x, 0.5)),
  aggregate(
    syncdif.spp[,c(3:(ncol(syncdif.spp)-3))],
    by = list(syncdif.spp$Group.2),
    FUN = function(x) quantile(x, 0.05)),
  aggregate(
    syncdif.spp[,c(3:(ncol(syncdif.spp)-3))],
    by = list(syncdif.spp$Group.2),
    FUN = function(x) quantile(x, 0.95)),
  aggregate(
    syncdif.spp[,c("x")],
    by = list(syncdif.spp$Group.2),
    mean)
)


############## progress


colnames(syncdif.spp.sum) = c("spcode", "lin.dif_med", "ar1.l.dif_med","lin.one_med", "ar1.l.one_med","lin.two_med", "ar1.l.two_med", "lin.dif_p5", "ar1.l.dif_p5", "lin.one_p5", "ar1.l.one_p5", "lin.two_p5", "ar1.l.two_p5", "lin.dif_p95", "ar1.l.dif_p95", "lin.one_p95", "ar1.l.one_p95", "lin.two_p95", "ar1.l.two_p95", "n")

syncdif.spp.sum = syncdif.spp.sum %>%
  mutate(
    spcode = as.character(spcode),
    fullsp = case_when(
      spcode == "eame" ~ "E. Meadowlark",
      spcode == "sppi" ~ "Sprague's Pipit",
      spcode == "noha" ~ "N. Harrier",
      spcode == "cclo" ~ "Chestnut-col. Longspur",
      spcode == "grsp" ~ "Grasshopper Sparrow",
      spcode == "upsa" ~ "Upland Sandpiper",
      spcode == "rnep" ~ "Ring-necked Pheasant",
      spcode == "vesp" ~ "Vesper Sparrow",
      spcode == "bobo" ~ "Bobolink",
      spcode == "lcsp" ~ "LeConte's Sparrow",
      spcode == "hola" ~ "Horned Lark",
      spcode == "savs" ~ "Savannah Sparrow",
      spcode == "weme" ~ "W. Meadowlark",
      spcode == "larb" ~ "Lark Bunting",
      spcode == "dick" ~ "Dickcissel",
      spcode == "bais" ~ "Baird's Sparrow",
      spcode == "sewr" ~ "Sedge Wren",
      spcode == "casp" ~ "Cassin's Sparrow",
      spcode == "lbcu" ~ "Long-billed Curlew",
      spcode == "mclo" ~ "McCown's Longspur",
      spcode == "All species" ~ "All species"
    )
  ) %>%
  mutate(
    order = case_when(spcode == "All species" ~ -0.5,
                      TRUE ~ lin.dif_med),
    fullsp = forcats::fct_reorder(fullsp, order, .desc = T)
  ) 

syncdif.spp.map = syncdif.spp.map %>%
  left_join(syncdif.spp.sum[,c("spcode", "fullsp")], by = "spcode") 

# Remove unneeded objects
rm(grid.common)
rm(grid.common2)
rm(mapsync.com)
rm(mapsync.com1)
rm(mapsync.com2)
rm(mapsync.dif.names)
rm(mapsync.dif)
rm(mapsync.dif1)
rm(mapsync.iterations)
rm(syncdif.all)
rm(syncdif.spp.prep)
rm(syncdif.spp)
rm(syncdif.split)
rm(syncdif)
rm(syncdif.sum.list)











else{
  # make a list of files containing the mean synchrony data sets
  mapsync.postmeds.files = sprintf(
    "scripts//amarel//%s.mapsync.postmeds.csv",
    paste0(
      species$spcode,
      ".",
      substr(species$startyr, 3, 4),
      ".",
      substr(species$endyr, 3, 4)
    )
  )
  
  # read in all mapsync data into list form
  mapsync.postmeds = lapply(mapsync.postmeds.files, function(x) read.csv(x))
  names(mapsync.postmeds) <-   paste0(species$spcode,
                                      ".",
                                      substr(species$startyr, 3, 4),
                                      ".",
                                      substr(species$endyr, 3, 4))
  
  # get a list of the grid cells in common between periods for each species (for filtering)
  gridcom.postmeds = function(x) {Reduce(intersect, list(
    unique(mapsync.postmeds[[x]]$grid),
    unique(mapsync.postmeds[[x + 19]]$grid)
  ))}
  
  grid.common.postmeds = lapply(1:19, gridcom.postmeds)
  names(grid.common.postmeds) = species[1:19,]$spcode
  
  grid.common.postmeds2 = c(grid.common.postmeds, grid.common.postmeds)
  
  # subset only the grid cells in common between the two time periods
  mapsync.postmeds.com = list()
  for(i in 1:38){
    mapsync.postmeds.com[[i]] = filter(mapsync.postmeds[[i]], grid %in% grid.common.postmeds2[[i]])
  }
  names(mapsync.postmeds.com) <- species$spcode
  
  # separate out period one, period two, and subtract the two to get change in synchrony
  mapsync.postmeds.com1 = mapsync.postmeds.com[1:19]
  mapsync.postmeds.com2 = mapsync.postmeds.com[20:38]
  mapsync.postmeds.dif1 = Map(function(listper2, listper1)
    listper2[, 4:9] - listper1[, 4:9],
    mapsync.postmeds.com2,
    mapsync.postmeds.com1)
  
  # add in iteration and grid variables back to synchrony change
  mapsync.postmeds.dif = Map(function(mainlist, listnames) cbind(listnames[,2:3], mainlist[,1:6]), 
                             mapsync.postmeds.dif1, mapsync.postmeds.com1)
  
  mapsync.postmeds.dif.names = list()
  for(i in 1:19){
    mapsync.postmeds.dif.names[[i]] = cbind.data.frame(species = species$spcode[i], mapsync.postmeds.dif[[i]])  
  }
  
  syncdif.postmeds.all = do.call("rbind", mapsync.postmeds.dif.names)
  
  syncdif.postmeds.sum = syncdif.postmeds.all %>%
    group_by(grid) %>%
    summarise(
      ar1 = mean(meancor.ar1),
      lin = mean(meancor.lin),
      dif = mean(meancor.lin),
      ar1.l = mean(meancor.ar1.l),
      lin.l = mean(meancor.lin.l),
      dif.l = mean(meancor.dif.l),
      spp.per.grid = length(species)
    ) %>%
    ungroup() %>%
    filter(spp.per.grid >= 3) %>%
    dplyr::select(-spp.per.grid) %>%
    left_join(mean.routes.wide, by = "grid") %>%
    mutate(lat = as.numeric(substr(grid,1,2))-0.5,
           lon = as.numeric(substr(grid,3,10))+0.5) 
}