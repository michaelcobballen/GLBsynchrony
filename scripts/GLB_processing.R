#### allsync: Combines function results (pairwise gridcell synchrony) for all grassland species. ####
allsync1 = rbind(
  grsp.66.91 = cbind.data.frame(bbsync("5460",1966,1991,"AR","sync.dist"),sp="GRSP",per="1966-1991"),
  grsp.92.17 = cbind.data.frame(bbsync("5460",1992,2017,"AR","sync.dist"),sp="GRSP",per="1992-2017"),
  bobo.66.91 = cbind.data.frame(bbsync("4940",1966,1991,"AR","sync.dist"),sp="BOBO",per="1966-1991"),
  bobo.92.17 = cbind.data.frame(bbsync("4940",1992,2017,"AR","sync.dist"),sp="BOBO",per="1992-2017"),
  eame.66.91 = cbind.data.frame(bbsync("5010",1966,1991,"AR","sync.dist"),sp="EAME",per="1966-1991"),
  eame.92.17 = cbind.data.frame(bbsync("5010",1992,2017,"AR","sync.dist"),sp="EAME",per="1992-2017")
)

#.rs.restartR()
# source("scripts/GLBsynchrony_libraries.R") # reloads all required libraries after restart
# putting this between batches might make it faster

# Batch 2
allsync2 = rbind(
  savs.66.91 = cbind.data.frame(bbsync("5420",1966,1991,"AR","sync.dist"),sp="SAVS",per="1966-1991"),
  savs.92.17 = cbind.data.frame(bbsync("5420",1992,2017,"AR","sync.dist"),sp="SAVS",per="1992-2017"),
  noha.66.91 = cbind.data.frame(bbsync("3310",1966,1991,"AR","sync.dist"),sp="NOHA",per="1966-1991"),
  noha.92.17 = cbind.data.frame(bbsync("3310",1992,2017,"AR","sync.dist"),sp="NOHA",per="1992-2017"),
  upsa.66.91 = cbind.data.frame(bbsync("2610",1966,1991,"AR","sync.dist"),sp="UPSA",per="1966-1991"),
  upsa.92.17 = cbind.data.frame(bbsync("2610",1992,2017,"AR","sync.dist"),sp="UPSA",per="1992-2017")
)

# Batch 3
allsync3 = rbind(
  lbcu.66.91 = cbind.data.frame(bbsync("2640",1966,1991,"AR","sync.dist"),sp="LBCU",per="1966-1991"),
  lbcu.92.17 = cbind.data.frame(bbsync("2640",1992,2017,"AR","sync.dist"),sp="LBCU",per="1992-2017"),
  hola.66.91 = cbind.data.frame(bbsync("4740",1966,1991,"AR","sync.dist"),sp="HOLA",per="1966-1991"),
  hola.92.17 = cbind.data.frame(bbsync("4740",1992,2017,"AR","sync.dist"),sp="HOLA",per="1992-2017"),
  hesp.66.91 = cbind.data.frame(bbsync("5470",1966,1991,"AR","sync.dist"),sp="HESP",per="1966-1991"),
  hesp.92.17 = cbind.data.frame(bbsync("5470",1992,2017,"AR","sync.dist"),sp="HESP",per="1992-2017")
)

# Batch 4
allsync4 = rbind(
  sewr.66.91 = cbind.data.frame(bbsync("7240",1966,1991,"AR","sync.dist"),sp="SEWR",per="1966-1991"),
  sewr.92.17 = cbind.data.frame(bbsync("7240",1992,2017,"AR","sync.dist"),sp="SEWR",per="1992-2017"),
  sppi.66.91 = cbind.data.frame(bbsync("7000",1966,1991,"AR","sync.dist"),sp="SPPI",per="1966-1991"),
  sppi.92.17 = cbind.data.frame(bbsync("7000",1992,2017,"AR","sync.dist"),sp="SPPI",per="1992-2017"),
  dick.66.91 = cbind.data.frame(bbsync("6040",1966,1991,"AR","sync.dist"),sp="DICK",per="1966-1991"),
  dick.92.17 = cbind.data.frame(bbsync("6040",1992,2017,"AR","sync.dist"),sp="DICK",per="1992-2017")
)

# Batch 5
allsync5 = rbind(
  vesp.66.91 = cbind.data.frame(bbsync("5400",1966,1991,"AR","sync.dist"),sp="VESP",per="1966-1991"),
  vesp.92.17 = cbind.data.frame(bbsync("5400",1992,2017,"AR","sync.dist"),sp="VESP",per="1992-2017"),
  larb.66.91 = cbind.data.frame(bbsync("6050",1966,1991,"AR","sync.dist"),sp="LARB",per="1966-1991"),
  larb.92.17 = cbind.data.frame(bbsync("6050",1992,2017,"AR","sync.dist"),sp="LARB",per="1992-2017"),
  bais.66.91 = cbind.data.frame(bbsync("5450",1966,1991,"AR","sync.dist"),sp="BAIS",per="1966-1991"),
  bais.92.17 = cbind.data.frame(bbsync("5450",1992,2017,"AR","sync.dist"),sp="BAIS",per="1992-2017")
)

# Batch 6
allsync6 = rbind(
  casp.66.91 = cbind.data.frame(bbsync("5780",1966,1991,"AR","sync.dist"),sp="CASP",per="1966-1991"),
  casp.92.17 = cbind.data.frame(bbsync("5780",1992,2017,"AR","sync.dist"),sp="CASP",per="1992-2017"),
  lesp.66.91 = cbind.data.frame(bbsync("5480",1966,1991,"AR","sync.dist"),sp="LESP",per="1966-1991"),
  lesp.92.17 = cbind.data.frame(bbsync("5480",1992,2017,"AR","sync.dist"),sp="LESP",per="1992-2017"),
  weme.66.91 = cbind.data.frame(bbsync("5011",1966,1991,"AR","sync.dist"),sp="WEME",per="1966-1991"),
  weme.92.17 = cbind.data.frame(bbsync("5011",1992,2017,"AR","sync.dist"),sp="WEME",per="1992-2017")
)

# Batch 7
allsync7 = rbind(
  rnep.66.91 = cbind.data.frame(bbsync("3091",1966,1991,"AR","sync.dist"),sp="RNEP",per="1966-1991"),
  rnep.92.17 = cbind.data.frame(bbsync("3091",1992,2017,"AR","sync.dist"),sp="RNEP",per="1992-2017"),
  mclo.66.91 = cbind.data.frame(bbsync("5390",1966,1991,"AR","sync.dist"),sp="MCLO",per="1966-1991"),
  mclo.92.17 = cbind.data.frame(bbsync("5390",1992,2017,"AR","sync.dist"),sp="MCLO",per="1992-2017"),
  cclo.66.91 = cbind.data.frame(bbsync("5380",1966,1991,"AR","sync.dist"),sp="CCLO",per="1966-1991"),
  cclo.92.17 = cbind.data.frame(bbsync("5380",1992,2017,"AR","sync.dist"),sp="CCLO",per="1992-2017")
)

# Batch 8
allsync8 = rbind(
  stgr.66.91 = cbind.data.frame(bbsync("3080",1966,1991,"AR","sync.dist"),sp="STGR",per="1966-1991"),
  stgr.92.17 = cbind.data.frame(bbsync("3080",1992,2017,"AR","sync.dist"),sp="STGR",per="1992-2017"),
  mopl.66.91 = cbind.data.frame(bbsync("2810",1966,1991,"AR","sync.dist"),sp="MOPL",per="1966-1991"),
  mopl.92.17 = cbind.data.frame(bbsync("2810",1992,2017,"AR","sync.dist"),sp="MOPL",per="1992-2017")
)

# FEHA, GRPC, BNOW, and SEOW didn't work - no qualifying cells

allsync = rbind(allsync1,allsync2,allsync3,allsync4,allsync5,allsync6,allsync7,allsync8) 

allsync = allsync %>%
  group_by(sp, per) %>%
  mutate(pairs = length(sp)) %>%
  ungroup()

#### mapsync: Creates moving average of all adjacent cells within 400 km ####
mapsync = 
  rbind(data.frame(allsync, grid = allsync$grid4A), data.frame(allsync, grid = allsync$grid4B)) %>%
  filter(dist<400) %>%
  group_by(grid, per, sp) %>%
  summarise(meancor = mean(cor), n = length(cor)) %>%
  ungroup() %>%
  mutate(grid_lat = as.numeric(substr(grid,1,2)), grid_lon = as.numeric(paste("-",substr(grid,4,10),sep=""))) %>%
  mutate(corcat = cut(meancor, breaks = c(-1,0,0.1,.2,.3,.4,.5,.6,.7))) %>%
  group_by(sp, per) %>%
  mutate(cells = length(sp)) %>%
  ungroup()

#### syncdiff: Calculate change in synchrony between periods for each species/grid cell. ####
mapsync.early = mapsync %>% filter(per=="1966-1991") %>% mutate(grid.sp = paste(grid,sp,sep="."))
mapsync.late = mapsync %>% filter(per=="1992-2017") %>% mutate(grid.sp = paste(grid,sp,sep="."))
syncdiff = full_join(mapsync.early, mapsync.late, by = "grid.sp") %>%
  filter(is.na(meancor.x)==F, is.na(meancor.y)==F) %>%
  mutate(syncdiff = meancor.y - meancor.x) %>%
  mutate(diffcat = cut(syncdiff, breaks = c(-.6,-.5,-.4,-.3,-.2,-.1,0,.1,.2,.3,.4,.5,.6))) %>%
  mutate(incdec = ifelse(syncdiff>0,"Increasing","Decreasing")) %>%
  dplyr::select(-grid.y,-grid_lat.y, -grid_lat.y)


#### meansyncdiff: Calculate mean change in synchrony across all GLB species by grid cell. ####
meansyncdiff = syncdiff %>%
  group_by(grid.x, grid_lat.x, grid_lon.x) %>%
  summarise(meansyncdiff = mean(syncdiff), n = length(sp.x)) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  mutate(meandiffcat = cut(meansyncdiff, breaks = c(-.3,-.25,-.2,-.15,-.1,-.05,0,.05,.1,.15,.2,.25,.3))) %>%
  mutate(meanincdec = ifelse(meansyncdiff>0,"Increasing","Decreasing")) %>%
  mutate(species = "All species")