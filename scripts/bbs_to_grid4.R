###############################################################################################
###
### An attempt to generate the BBS indices myself
### for 2x2 degree grid-cell strata via methods in Michel et al. Ecography (2016)
### to use in the conditional autoregressive (CAR) spatial model described with code in Smith et al. PlosOne (2015)
###
###############################################################################################

# see Links at bottom for data sources, papers, etc.

########################################################################################

library(dplyr)
library(stringr)
library(sp)
library(rgdal)
library(maptools)

### data prep for formatting route-level data into 2x2 degree lat/long grid for WinBugs
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

# functions for identifying even and odd numbers
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

brw = b %>%
  left_join(r,by = 'rid') %>%
  left_join(w,by = 'RouteDataID') %>%
  filter(RunType==1,RPID.x ==101) %>%   
  select(RouteDataID, CountryNum = CountryNum.x, StateNum = StateNum.x, Route = Route.x,
         Year = Year.x, AOU, StopTotal, rid = rid.x, RouteName, lat = Latitude, lon = Longitude, 
         Stratum, BCR, ObsN) %>%
  mutate(sp.number = as.numeric(ifelse(str_detect(AOU,"5460"),StopTotal,0))) %>%
  mutate(grid4_lat = ifelse(is.odd(trunc(lat)),trunc(lat)+1,trunc(lat))) %>%
  mutate(grid4_lon = ifelse(is.odd(trunc(lon)),trunc(lon)-1,trunc(lon))) %>%
  mutate(grid4 = paste(grid4_lat,grid4_lon,sep=""))

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
  filter(is.na(lat)==F) %>% # this excludes ~ 20 routes without attribute data
  group_by(grid4) %>%
  mutate(ObsN.rid.count = as.numeric(as.factor(ObsN.rid))) %>%
  ungroup()

brw.route.year.92.17 = brw %>%
  filter(Year>1991) %>%
  group_by(RouteDataID) %>%
  summarise(sp.number = sum(sp.number)) %>%
  left_join(distinct(dplyr::select(brw,-AOU,-sp.number,-StopTotal)),by="RouteDataID") %>%
  mutate(ObsN.rid = paste(ObsN,".",rid,sep="")) %>%
  group_by(ObsN.rid) %>%
  arrange(ObsN.rid,Year) %>%
  mutate(ObsN.rid.year = row_number(), yr = Year-1965) %>%  
  mutate(ObsN.rid.first = as.numeric(ifelse(ObsN.rid.year==1,1,0))) %>%
  ungroup() %>%
  filter(is.na(lat)==F) %>% # this excludes ~ 20 routes without attribute data
  group_by(grid4) %>%
  mutate(ObsN.rid.count = as.numeric(as.factor(ObsN.rid))) %>%
  ungroup()

brw.route = brw.route.year %>%
  group_by(rid) %>%
  summarise(sp.yrs = sum(sp.number>0), yrs.data = length(sp.number), sp.mean = mean(sp.number),
            nobservers = max(ObsN.rid.count)) %>%
  left_join(distinct(select(brw.route.year,-sp.number,-Year,-RouteDataID,-ObsN,-ObsN.rid,-ObsN.rid.year,-yr,-ObsN.rid.first,-ObsN.rid.count)),by="rid") %>%
  mutate(sp.present = as.numeric(ifelse(sp.yrs>0,1,0))) %>%
  arrange(lat) %>%
  ungroup() 

brw.route.92.17 = brw.route.year %>%
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

grid4.summary = brw.route %>%
  group_by(grid4) %>%
  summarise(grid4_lat=mean(grid4_lat), grid4_lon=mean(grid4_lon), num.routes = length(sp.present), 
            sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
            ten_yr_rtes = sum(yrs.data>9), nobservers = max(nobservers)) %>%
  mutate(sp.two.rtes.yes = as.numeric(ifelse(num_sp_rtes > 1,1,0))) %>%
  mutate(ten.years.yes = as.numeric(ifelse(ten_yr_rtes > 0,1,0))) %>%
  mutate(include_strat = as.numeric(ifelse(sp.two.rtes.yes + ten.years.yes == 2,1,0))) %>%
  mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
  select(-sp.two.rtes.yes,-ten.years.yes)

grid4.summary.92.17 = brw.route.92.17 %>%
  group_by(grid4) %>%
  summarise(grid4_lat=mean(grid4_lat), grid4_lon=mean(grid4_lon), num.routes = length(sp.present), 
            sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
            ten_yr_rtes = sum(yrs.data>9), nobservers = max(nobservers)) %>%
  mutate(sp.two.rtes.yes = as.numeric(ifelse(num_sp_rtes > 1,1,0))) %>%
  mutate(ten.years.yes = as.numeric(ifelse(ten_yr_rtes > 0,1,0))) %>%
  mutate(include_strat = as.numeric(ifelse(sp.two.rtes.yes + ten.years.yes == 2,1,0))) %>%
  mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
  select(-sp.two.rtes.yes,-ten.years.yes)

grid4.summary.spatial = sp::SpatialPointsDataFrame(grid4.summary[,c("grid4_lon","grid4_lat")],grid4.summary)
proj4string(grid4.summary.spatial) = CRS("+init=epsg:4326")
plot(grid4.summary.spatial)
# writeOGR(grid4.summary.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_grid4", driver="ESRI Shapefile")

route.year.grid4 = brw.route.year %>%
  left_join(select(grid4.summary, grid4, include_strat), by = "grid4") %>%
  filter(include_strat==1) %>%
  select(count = sp.number, grid4, obser = ObsN.rid, firstyr = ObsN.rid.first, 
         year = yr, -include_strat, grid4_lat, grid4_lon) %>%
  mutate(grid4 = gsub("-","_",grid4))

route.year.grid4.92.17 = brw.route.year.92.17 %>%
  left_join(select(grid4.summary, grid4, include_strat), by = "grid4") %>%
  filter(include_strat==1) %>%
  select(count = sp.number, grid4, obser = ObsN.rid, firstyr = ObsN.rid.first, 
         year = yr, -include_strat, grid4_lat, grid4_lon) %>%
  mutate(grid4 = gsub("-","_",grid4))
  
route.year.grid4.spatial = SpatialPointsDataFrame(route.year.grid4[,c("grid4_lon","grid4_lat")],route.year.grid4)
proj4string(route.year.grid4.spatial) = CRS("+init=epsg:4326")
plot(route.year.grid4.spatial)
# writeOGR(route.year.grid.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_route_yrs", driver="ESRI Shapefile")

route.year.grid4.92.17.spatial = SpatialPointsDataFrame(route.year.grid4.92.17[,c("grid4_lon","grid4_lat")],route.year.grid4.92.17)
proj4string(route.year.grid4.92.17.spatial) = CRS("+init=epsg:4326")
plot(route.year.grid4.92.17.spatial)
# writeOGR(route.year.grid.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_route_yrs", driver="ESRI Shapefile")


#### TO SUMMARIZE DATA BY YEAR AND GRID

grid4.year = route.year.grid4 %>%
  group_by(grid4, year, grid4_lat, grid4_lon) %>%
  summarise(mean.count = mean(count)) %>%
  ungroup() %>%
  mutate(year = year+1965, grid4.year = paste(grid4,year,sep="."))

grid4.year.92.17 = route.year.grid4.92.17 %>%
  group_by(grid4, year, grid4_lat, grid4_lon) %>%
  summarise(mean.count = mean(count)) %>%
  ungroup() %>%
  mutate(year = year+1965, grid4.year = paste(grid4,year,sep="."))



### create criteria for including only grid cells with < XX years of zeros

grid4.include2 = grid4.year %>%
  group_by(grid4) %>%
  summarise(num.zeros = sum(mean.count==0)) %>%
  mutate(include2 = ifelse(num.zeros >4,0,1)) %>%
  ungroup()

grid4.include2.92.17 = grid4.year.92.17 %>%
  group_by(grid4) %>%
  summarise(num.zeros = sum(mean.count==0)) %>%
  mutate(include2 = ifelse(num.zeros >4,0,1)) %>%
  ungroup()

# adding the second inclusion criteria to grid.summary
# subsetting just the included grid cells based on Michel et al. criteria*
## *Michel et al. criteria (Ecography 2016): at least 2 routes with species, and at least one route with 10+ years of data
## additional restrictions are based on number of years with zeros

grid4.summary2 = grid4.summary %>%
  mutate(grid4 = gsub("-","_",grid4)) %>%
  left_join(grid4.include2, by = "grid4") %>%
  filter(include_strat==1 & include2 == 1)

grid4.summary2.92.17 = grid4.summary.92.17 %>%
  mutate(grid4 = gsub("-","_",grid4)) %>%
  left_join(grid4.include2.92.17, by = "grid4") %>%
  filter(include_strat==1 & include2 == 1)

# make into spatial object

grid4.grsp.spatial = SpatialPointsDataFrame(grid4.summary2[,c("grid4_lon","grid4_lat")],grid4.summary2)
proj4string(grid4.grsp.spatial) = CRS("+init=epsg:4326")
plot(grid4.grsp.spatial)
# writeOGR(grid.grsp.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid", driver="ESRI Shapefile")

grid4.grsp.spatial.92.17 = SpatialPointsDataFrame(grid4.summary2.92.17[,c("grid4_lon","grid4_lat")],grid4.summary2.92.17)
proj4string(grid4.grsp.spatial.92.17) = CRS("+init=epsg:4326")
plot(grid4.grsp.spatial.92.17)
# writeOGR(grid.grsp.spatial.92.17, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid", driver="ESRI Shapefile")


#### attempting to make a grid of spatial polygons that match up with data
rr <- raster::raster(ext = extent(-177, -52, 24, 69), res=c(2,2))
values(rr) <- 1:ncell(rr)
p <- rasterToPolygons(rr)
attributes = tibble(grid4_lon = coordinates(p)[,1], grid4_lat = coordinates(p)[,2]) %>%
  mutate(join = 1, grid4 = paste(trunc(grid4_lat),trunc(grid4_lon),sep="")) %>%
  mutate(grid4 = gsub("-","_",grid4)) %>%
  left_join(grid4.summary2.92.17, by="grid4") 

p$grid4 = attributes$grid4
p$include_strat = attributes$include_strat
p$include2 = attributes$include2
p$nonzeroweight = attributes$nonzeroweight
p$nobservers = attributes$nobservers
p$strata = log(attributes$sp.mean)


grid4.grsp.92.17 = p[p$include_strat==1 & is.na(p$include_strat)==F,]
grid4.grsp.92.17@data = grid4.grsp@data %>%
  select(layer, grid4, include_strat, include2, nonzeroweight, nobservers, strata) %>%
  mutate(strat = 1:nrow(grid4.grsp.92.17@data))
head(grid4.grsp.92.17@data); nrow(grid4.grsp.92.17@data)

plot(grid4.grsp.92.17)
# writeOGR(grid4.grsp, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_poly", driver="ESRI Shapefile")

proj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=10000000 +y_0=10000000 +datum=NAD83 +units=m +no_defs"
grid4.grsp.92.17.lamb0 = spTransform(grid4.grsp.92.17, proj)
plot(grid4.grsp.92.17.lamb0, col="gray")
head(grid4.grsp.92.17.lamb0); nrow(grid4.grsp.92.17.lamb0)
# writeOGR(grid4.grsp.92.17.lamb0, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_92_17_poly_lamb1", driver="ESRI Shapefile")


# changing names of the ID variable of the polygons to grid labels so names are visible in geobugs map
# https://eco-data-science.github.io/spatial_analysis2_R/
sapply(slot(grid4.grsp.92.17.lamb0, "polygons"), function(x) slot(x, "ID"))

for (i in 1:length(slot(grid4.grsp.92.17.lamb0, "polygons"))){
  slot(slot(grid4.grsp.92.17.lamb0, "polygons")[[i]], "ID") = grid4.grsp.92.17.lamb0$grid4[i]
}
sapply(slot(grid4.grsp.92.17.lamb0, "polygons"), function(x) slot(x, "ID"))


### converting spatial data to geobugs format if needed

#sp2WB(grid.grsp.lamb0, "grid.grsp.lamb3.txt")



### attaching the include2 criteria back onto the summary by grid4 and year

grid4.year = grid4.year %>%
  left_join(grid4.include2,by = "grid4") %>%
  filter(include2 == 1)

grid4.year.92.17 = grid4.year.92.17 %>%
  left_join(grid4.include2.92.17,by = "grid4") %>%
  filter(include2 == 1)




### Links to data sources, code, papers, etc.
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
