###############################################################################################
###
### This part is an attempt to generate the BBS indices myself
### For 1-degree-grid-cell strata via methods in Michel et al. Ecography (2016)
### i.e., the conditional Autoregressive spatial model described with code in Smith et al. PlosOne (2015)
###
###############################################################################################

# Links:
# BBS Route-level data fields: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/Summary.txt
# BBS RouteID info: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/RegionCodes.txt
# BBS phys strata map: https://www.pwrc.usgs.gov/bbs/StrataNames/strata_map_small.htm
# BBS Route info txt with BCR and strata codes: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/RouteInf.txt
# BBS BCR shapefile: http://www.pwrc.usgs.gov/bba/index.cfm?fa=bba.getdata
# Smith et al. sup data for CAR model and indices of aerial insectivores: 
#       https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130768#sec016
#       https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0130768.s005&type=supplementary
# Michel et al. README file: https://datadryad.org/bitstream/handle/10255/dryad.96773/README.txt?sequence=2
# Michel et al. paper with link to sup code files: https://onlinelibrary.wiley.com/doi/abs/10.1111/ecog.01798

# Michel et al. keep all strata with "two or more routes where the species has been observed 
# and at least one of which with >= 10 years of data


########################################################################################

library(dplyr)
library(stringr)
library(sp)
library(rgdal)
library(maptools)

### data prep for formatting route-level data for WinBugs
# BBS fields: ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/NABBS_DataFieldDefinitions_1996-2017.xml
#  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS

# read and combine all route-level data stored in separate files by state

path <- "data//states//"
files <- list.files(path=path, pattern="*.csv")
for(file in files)
{
  perpos <- which(strsplit(file, "")[[1]]==".")
  assign(
    gsub(" ","",substr(file, 1, perpos-1)), 
    read.csv(paste(path,file,sep=""),header=T,skip=0))
}


b = rbind.data.frame(Alabama,Alaska,Alberta,Arizona,Arkansa,BritCol,Califor,Colorad,Connect,Delawar,
                   Florida,Georgia,Idaho,Illinoi,Indiana,Iowa,Kansas,Kentuck, Louisia,Maine,
                   Manitob,Marylan,Massach,Michiga,Minneso,Mississ,Missour,Montana,NBrunsw,
                   NCaroli,NDakota,Nebrask,Nevada,Newfoun,NHampsh,NJersey,NMexico,NovaSco,Nunavut,
                   NWTerri,NYork,Ohio,Oklahom,Ontario,Oregon,PEI,Pennsyl,Quebec,RhodeIs,
                   Saskatc,SCaroli,SDakota,Tenness,Texas,Utah,Vermont,Virgini,W_Virgi,Washing,
                   Wiscons,Wyoming,Yukon) %>% 
        mutate(rid = paste(StateNum,".",Route,sep=""))
r = read.csv("data/routes.csv") %>% 
        mutate(statenum=as.character(statenum),rid = paste(statenum,".",Route,sep=""))
w = read.csv("data/weather.csv") %>% 
        mutate(rid = paste(StateNum,".",Route,sep=""))

brw = b %>%
  left_join(r,by = 'rid') %>%
  left_join(w,by = 'RouteDataID') %>%
  filter(RunType==1,RPID.x ==101) %>%   
  select(RouteDataID, CountryNum = CountryNum.x, StateNum = StateNum.x, Route = Route.x,
         Year = Year.x, AOU, StopTotal, rid = rid.x, RouteName, lat = Latitude, lon = Longitude, 
         Stratum, BCR, ObsN) %>%
  mutate(sp.number = as.numeric(ifelse(str_detect(AOU,"5460"),StopTotal,0))) %>%
  mutate(grid = ifelse(nchar(trunc(lon))==3,paste(substr(lat,1,2),substr(lon,1,3),sep=""),
                       paste(substr(lat,1,2),substr(lon,1,4),sep=""))) %>%
  mutate(grid_lat = as.numeric(ifelse(nchar(trunc(lat))==3,paste(substr(lat,1,2),".5",sep=""),
                                      paste(substr(lat,1,2),".5",sep="")))) %>%
  mutate(grid_lon = as.numeric(ifelse(nchar(trunc(lon))==3,paste(substr(lon,1,3),".5",sep=""),
                                      paste(substr(lon,1,4),".5",sep=""))))


brw.route.year = brw %>%
  group_by(RouteDataID) %>%
  summarise(sp.number = sum(sp.number)) %>%
  left_join(distinct(dplyr::select(brw,-AOU,-sp.number,-StopTotal)),by="RouteDataID") %>%
  mutate(ObsN.rid = paste(ObsN,".",rid,sep="")) %>%
  group_by(ObsN.rid) %>%
  arrange(ObsN.rid,Year) %>%
  mutate(ObsN.rid.year = row_number(), yr = Year-1965) %>%  
  mutate(ObsN.rid.first = as.numeric(ifelse(ObsN.rid.year==1,1,0))) %>%
  ungroup() %>%
  filter(is.na(lat)==F) %>%
  group_by(grid) %>%
  mutate(ObsN.rid.count = as.numeric(as.factor(ObsN.rid))) %>%
  ungroup()


# the filter command in line 90 excludes 20 routes with no associated route time/location info in "r" file
# use this next line to test they are gone
# test = brw.route[is.na(brw.route$lat),]

brw.route = brw.route.year %>%
  group_by(rid) %>%
  summarise(sp.yrs = sum(sp.number>0), yrs.data = length(sp.number), sp.mean = mean(sp.number),
            nobservers = max(ObsN.rid.count)) %>%
  left_join(distinct(select(brw.route.year,-sp.number,-Year,-RouteDataID,-ObsN,-ObsN.rid,-ObsN.rid.year,-yr,-ObsN.rid.first,-ObsN.rid.count)),by="rid") %>%
  mutate(sp.present = as.numeric(ifelse(sp.yrs>0,1,0))) %>%
  arrange(lat) %>%
  ungroup() 


brw.route.spatial = sp::SpatialPointsDataFrame(brw.route[,c("lon","lat")],brw.route)
proj4string(brw.route.spatial) = CRS("+init=epsg:4326")
plot(brw.route.spatial)
# writeOGR(brw.route.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_routes1", driver="ESRI Shapefile")

grid.summary = brw.route %>%
  group_by(grid) %>%
  summarise(grid_lat=mean(grid_lat), grid_lon=mean(grid_lon), num.routes = length(sp.present), 
            sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
            ten_yr_rtes = sum(yrs.data>9), nobservers = max(nobservers)) %>%
  mutate(sp.two.rtes.yes = as.numeric(ifelse(num_sp_rtes > 1,1,0))) %>%
  mutate(ten.years.yes = as.numeric(ifelse(ten_yr_rtes > 0,1,0))) %>%
  mutate(include_strat = as.numeric(ifelse(sp.two.rtes.yes + ten.years.yes == 2,1,0))) %>%
  mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
  select(-sp.two.rtes.yes,-ten.years.yes)

grid.summary.spatial = sp::SpatialPointsDataFrame(grid.summary[,c("grid_lon","grid_lat")],grid.summary)
proj4string(grid.summary.spatial) = CRS("+init=epsg:4326")
plot(grid.summary.spatial)
# writeOGR(grid.summary.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_grid", driver="ESRI Shapefile")

route.year.grid = brw.route.year %>%
  left_join(select(grid.summary, grid, include_strat), by = "grid") %>%
  filter(include_strat==1) %>%
  select(count = sp.number, grid, obser = ObsN.rid, firstyr = ObsN.rid.first, 
         year = yr, -include_strat, grid_lat, grid_lon) %>%
  mutate(grid = gsub("-","_",grid))
  
route.year.grid.spatial = SpatialPointsDataFrame(route.year.grid[,c("grid_lon","grid_lat")],route.year.grid)
proj4string(route.year.grid.spatial) = CRS("+init=epsg:4326")
plot(route.year.grid.spatial)
# writeOGR(route.year.grid.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_route_yrs", driver="ESRI Shapefile")

# just the grid cells included in the analysis (Michel et al. criteria - Ecography 2016)
grid.grsp.spatial = SpatialPointsDataFrame(grid.summary[grid.summary$include_strat==1,c("grid_lon","grid_lat")],grid.summary[grid.summary$include_strat==1,])
proj4string(grid.grsp.spatial) = CRS("+init=epsg:4326")
plot(grid.grsp.spatial)
# writeOGR(grid.grsp.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid", driver="ESRI Shapefile")

#### attempting to make a grid of spatial polygons that match up with data
rr <- raster::raster(ext = extent(-177, -52, 24, 69), res=c(1,1))
values(rr) <- 1:ncell(rr)
p <- rasterToPolygons(rr)
attributes = dplyr::tibble(grid_lon = coordinates(p)[,1], grid_lat = coordinates(p)[,2]) %>%
  mutate(join = 1, grid = paste(trunc(grid_lat),trunc(grid_lon),sep="")) %>%
  left_join(filter(grid.summary,include_strat==1), by="grid") 

p$grid = attributes$grid
p$include_strat = attributes$include_strat
p$nonzeroweight = attributes$nonzeroweight
p$nobservers = attributes$nobservers
p$strata = log(attributes$sp.mean)


grid.grsp = p[p$include_strat==1 & is.na(p$include_strat)==F,]
grid.grsp@data = grid.grsp@data %>%
  select(layer, grid, include_strat, nonzeroweight, nobservers,strata) %>%
  mutate(grid = gsub("-","_",grid), strat = 1:558)
head(grid.grsp@data)

plot(grid.grsp)
# writeOGR(grid.grsp, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid_poly", driver="ESRI Shapefile")

proj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=10000000 +y_0=10000000 +datum=NAD83 +units=m +no_defs"
grid.grsp.lamb0 = spTransform(grid.grsp, proj)
plot(grid.grsp.lamb0)
head(grid.grsp.lamb0); nrow(grid.grsp.lamb0)
# writeOGR(grid.grsp.lamb0, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid_poly_lamb6", driver="ESRI Shapefile")


# changing names of the ID variable of the polygons to grid labels so names are visible in geobugs map
# https://eco-data-science.github.io/spatial_analysis2_R/
sapply(slot(grid.grsp.lamb0, "polygons"), function(x) slot(x, "ID"))

for (i in 1:length(slot(grid.grsp.lamb0, "polygons"))){
  slot(slot(grid.grsp.lamb0, "polygons")[[i]], "ID") = grid.grsp.lamb0$grid[i]
}
sapply(slot(grid.grsp.lamb0, "polygons"), function(x) slot(x, "ID"))


### converting spatial data to geobugs format if needed

#sp2WB(grid.grsp.lamb0, "grid.grsp.lamb3.txt")


#### TO SUMMARIZE DATA BY YEAR AND GRID
route.year.grid = route.year.grid %>%
  left_join(grid.grsp.lamb0@data, by = "grid")

grid.year = route.year.grid %>%
  group_by(grid, year, grid_lat,grid_lon) %>%
  summarise(mean.count = mean(count)) %>%
  ungroup() %>%
  mutate(year = year+1965, grid.year = paste(grid,year,sep="."))


