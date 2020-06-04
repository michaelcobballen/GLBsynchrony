bbsync = function(spcode, yearstart, yearend, method, output){
  
  # if running these as part of the function, then put "#" before them; if not, then remove the "#"
  #spcode = "5460" # AOU species number; #  spcode = s05460 for GRSP, s05010 for EAME, s04940 for BOBO, s05420 for SAVS
  #yearstart = 1992
  #yearend = 2017
  #method = "AR" # "AR" = residuals from 1st order autoregressive model, "change" = first difference of logged (count + 1)
  minyearsdata = 10 # each grid cell contains at least one route with at least this many years of data (10 is based on criteria of Michel et al. Ecography 2016)
  minsproutes = 2 # each grid cell contains at least this many routes at which the species was detected (2 is based on criteria of Michel et al. Ecography 2016)
  maxzeros = 5 # maximum number of years allowable where the mean count for the species is zero (I chose five)
  minyearscor = 20 # minimum number of years two grid cells must have in common to calculate a correlation coefficient between them (Hanski used 7; I chose 20, more conservative)
  #output = "grid.map" #or... CAR.data, sync.dist, grid.year, route.year
  
  # functions for identifying even and odd numbers
  is.even <- function(x) x %% 2 == 0
  is.odd <- function(x) x %% 2 != 0
  
  # combining the bird, route, and weather data
  # assigning routes to named grid cells 
  # creating a variable for abundance of target species
  brw = b %>%
    left_join(r,by = 'rid') %>%
    left_join(w,by = 'RouteDataID') %>%
    filter(RunType==1,RPID.x ==101) %>%   
    dplyr::select(RouteDataID, CountryNum = CountryNum.x, StateNum = StateNum.x, Route = Route.x,
                  Year = Year.x, AOU, StopTotal, rid = rid.x, RouteName, lat = Latitude, lon = Longitude, 
                  Stratum, BCR, ObsN) %>%
    mutate(sp.number = as.numeric(ifelse(str_detect(AOU,spcode),StopTotal,0))) %>%
    mutate(grid4_lat = ifelse(is.odd(trunc(lat)),trunc(lat)+1,trunc(lat))) %>%
    mutate(grid4_lon = ifelse(is.odd(trunc(lon)),trunc(lon)-1,trunc(lon))) %>%
    mutate(grid4 = paste(grid4_lat,grid4_lon,sep=""))
  
  # summarizing data for species of interest by route-year
  brw.route.year = brw %>%
    filter(Year>yearstart-1, Year<yearend+1) %>%
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
  
  # summarizing route-year data by route 
  # including # years species present & # years data on each route
  brw.route = brw.route.year %>%
    group_by(rid) %>%
    summarise(sp.yrs = sum(sp.number>0), yrs.data = length(sp.number), sp.mean = mean(sp.number),
              nobservers = max(ObsN.rid.count)) %>%
    left_join(distinct(dplyr::select(brw.route.year,-sp.number,-Year,-RouteDataID,-ObsN,-ObsN.rid,-ObsN.rid.year,-yr,-ObsN.rid.first,-ObsN.rid.count)),by="rid") %>%
    mutate(sp.present = as.numeric(ifelse(sp.yrs>0,1,0))) %>%
    arrange(lat) %>%
    ungroup() 
  
  # summarize data by 2x2 grid cell
  # including whether cells meet the "minyearsdata" and "minsproutes" criteria of Michel et al. (Ecography)
  grid4.summary = brw.route %>%
    group_by(grid4) %>%
    summarise(grid4_lat=mean(grid4_lat), grid4_lon=mean(grid4_lon), num.routes = length(sp.present), 
              sp.mean = mean(sp.mean), num_sp_rtes = sum(sp.present), 
              min_yr_rtes = sum(yrs.data>=minyearsdata), nobservers = max(nobservers)) %>%
    mutate(sp.min.rtes.yes = as.numeric(ifelse(num_sp_rtes >= minsproutes,1,0))) %>%
    mutate(min.years.yes = as.numeric(ifelse(min_yr_rtes > 0,1,0))) %>%
    mutate(include_strat = as.numeric(ifelse(sp.min.rtes.yes + min.years.yes == 2,1,0))) %>%
    mutate(nonzeroweight = num_sp_rtes/num.routes) %>%
    dplyr::select(-sp.min.rtes.yes,-min.years.yes)
  
  # add grid-level information to the route-year data
  # remove routes in cells excluded based on criteria in Michel et al. (Ecography)
  route.year.grid4 = brw.route.year %>%
    left_join(dplyr::select(grid4.summary, grid4, include_strat), by = "grid4") %>%
    filter(include_strat==1) %>%
    dplyr::select(count = sp.number, grid4, obser = ObsN.rid.count, firstyr = ObsN.rid.first, 
                  year = yr, -include_strat, grid4_lat, grid4_lon) %>%
    mutate(grid4 = gsub("-","_",grid4))
  
  # summarize route-year data by grid-year
  grid4.year = route.year.grid4 %>%
    group_by(grid4, year, grid4_lat, grid4_lon) %>%
    summarise(mean.count = mean(count)) %>%
    ungroup() %>%
    mutate(year = year+1965, grid4.year = paste(grid4,year,sep="."))
  
  ### create criteria for including only grid cells with < XX years of zeros
  grid4.include2 = grid4.year %>%
    group_by(grid4) %>%
    summarise(num.zeros = sum(mean.count==0)) %>%
    mutate(include2 = ifelse(num.zeros >maxzeros-1,0,1)) %>%
    ungroup()
  
  # subsetting route.year.grid4 based on the second inclusion criteria (number of years with zeros)
  # cells already removed based on the Michel et al. criteria
  # this is the data file to be used in the OpenBUGS conditional autoregressive (CAR) model via Smith et al. PlosONE 2015 code
  route.year.grid4b = route.year.grid4 %>%
    left_join(grid4.include2, by = "grid4") %>%
    filter(include2 == 1)
  
  # subsetting grid4.summary based on both criteria (Michel et al. and number of years with zeros) 
  grid4.summary2 = grid4.summary %>%
    mutate(grid4 = gsub("-","_",grid4)) %>%
    left_join(grid4.include2, by = "grid4") %>%
    filter(include_strat==1 & include2 == 1)
  
  ### subsetting the grid-year data based on 2nd criteria (number of years with zeros)
  # this had already been subsetted to exclude based on Michel et al. criteria 
  grid4.year = grid4.year %>%
    left_join(grid4.include2,by = "grid4") %>%
    filter(include2 == 1)
  
  ### calculating inter-grid distances for 2 X 2 degree lat/long grid
  grid4.spatial = SpatialPointsDataFrame(distinct(grid4.year[,c("grid4_lon","grid4_lat")]),
                                         distinct(grid4.year[,c("grid4_lon","grid4_lat","grid4")]))
  proj4string(grid4.spatial) = CRS("+init=epsg:4326")
  #plot(grid4.spatial)
  
  df.grid4 = data.frame(name = grid4.spatial$grid4, 
                        lat = coordinates(grid4.spatial)[,2],
                        lon = coordinates(grid4.spatial)[,1])
  
  #calculate Haversine distances between grid cells
  dmat.grid4 = round(distm(df.grid4[,c("lon","lat")])/1000)
  colnames(dmat.grid4) <- df.grid4$name; head(dmat.grid4)
  dmat.grid4B = data.frame(grid4B=df.grid4$name, dmat.grid4)
  
  dgrid4 = data.frame(reshape2::melt(dmat.grid4B)) %>%
    filter(value!=0) %>%
    dplyr::select(grid4B,grid4A = variable, dist = value) %>%
    mutate(grid4A = gsub("X","",grid4A)) %>%
    mutate(pairid = paste(grid4B,".",grid4A,sep=""))
  
  ### synchrony data formatting
  g4 = dgrid4
  
  # make a data frame of all possible grid/year combinations
  grid4.list = unique(grid4.year$grid4)
  year.list = (yearstart:yearend)
  grid4.year.list = expand.grid(grid4 = grid4.list,year = year.list) %>%
    arrange(grid4,year) %>%
    mutate(grid4.year = paste(grid4,year,sep="."))
  
  # attach the bird data to the complete list and format for correlation analysis
  grid4.year.full = grid4.year.list %>%
    left_join(grid4.year, by = "grid4.year") %>%
    dplyr::select(grid4=grid4.x, year=year.x, grid4.year, mean.count) %>%
    mutate(log.count = log(mean.count+1), logcountlag1 = lag(log.count)) %>% # log+1 transform count ala Hanski and create lag for AR and differencing analyses
    mutate(logcountlag1 = as.numeric(ifelse(year==yearstart,"NA",logcountlag1))) %>% # assign "NA" to first year which should have no lagged value
    filter(is.na(log.count) == F & is.na(logcountlag1) == F) %>% # temporarily remove NAs to get residuals and first diffs of logged counts
    mutate(change = log.count - logcountlag1) %>% # first diff of logged counts, ~ % change
    group_by(grid4) %>%
    mutate(AR1.R = resid(lm(log.count~logcountlag1))) %>% # calculate residuals from first order autoregression (Hanski-style)
    ungroup() %>%
    right_join(grid4.year.list, by = "grid4.year") %>% # merge back with the full grid/year list to allow matching in correlation analyses
    dplyr::select(grid4 = grid4.y, year = year.y, grid4.year, mean.count, log.count,logcountlag1, change, AR = AR1.R)
  
  ### put bird population data into 'wide' format for correlation analysis (rows = years, columns = grid cells)
  grid4.bird.mat = 
    grid4.year.full %>%
    dplyr::select(grid4,year,method) %>%
    spread(grid4,method) %>%
    filter(year > yearstart & year <= yearend) %>%
    # removes first year which have NAs for all calculations involving lagged counts
    dplyr::select(-year)
  
  # create matrix of Pearson correlations among grid cells
  grid4.bird.test = 
    grid4.bird.mat %>%
    corr.test(adjust="none", ci = F)
  
  # ntest checks if the n column of correlation test has length 1
  # allows us to deal with it differently in the next section
  # n has length 1 when all pairwise correlations have a sample size = yearend-yearstart (i.e., a quirk of the cortest function)
  # ...e.g., in Eastern Meadowlark, 1992-2017
  ntest = nrow(reshape2::melt(grid4.bird.test$n)) 
  
  # reshaping the giant correlation matrix into 'long' format for easier graphing and analysis
  grid4.bird.cor = 
    data.frame(reshape2::melt(grid4.bird.test$r), 
               reshape2::melt(grid4.bird.test$p)[3],
               n = ifelse(ntest > 1, reshape2::melt(grid4.bird.test$n)[3],
                          as.list(yearend-yearstart))[[1]]) %>%
    dplyr::select(grid4A=Var1, grid4B=Var2, cor=value, pval=value.1, n) %>%
    mutate(pairid = paste(grid4A,".",grid4B,sep="")) %>%
    filter(n>minyearscor-1) %>% # excludes grid pairs with less than XX years in common (Hanski's criteria was 7)
    filter(grid4A!=grid4B) %>% # excludes self-correlations
    filter(is.na(pval)==F) %>% # excludes NA correlations resulting from too many zeros if these weren't removed earlier - they were! So nothing changes.
    distinct(pval, .keep_all = TRUE) # should cut number in half, and yes it does! (eventually add check and warning message if it doesn't)
  
  # adding in the distance between grid cells again
  grid4.bird.cordist = grid4.bird.cor %>%
    left_join(g4,by='pairid') %>%
    dplyr::select(pairid, grid4A=grid4A.x, grid4B=grid4B.x, cor,pval,n,dist) %>%
    mutate(sig = ifelse(pval < 0.05,"Y","N"), sig1 = as.numeric(ifelse(pval < 0.05,"1","0"))) %>%
    distinct(cor, .keep_all = TRUE) 
  
  ### making spatial grid cell polygons to extract spatial data for use in the OpenBUGS CAR model via Smith et al. PlosONE 2015 code
  if(output == "grid.map" | output == "CAR.data"){
    
    # make a grid of spatial polygons to match up with data
    library(raster)
    rr <- raster::raster(ext = extent(-177, -52, 24, 69), res=c(2,2))
    values(rr) <- 1:ncell(rr)
    p <- rasterToPolygons(rr)
    attributes = tibble(grid4_lon = coordinates(p)[,1], grid4_lat = coordinates(p)[,2]) %>%
      mutate(join = 1, grid4 = paste(trunc(grid4_lat),trunc(grid4_lon),sep="")) %>%
      mutate(grid4 = gsub("-","_",grid4)) %>%
      left_join(grid4.summary2, by="grid4") 
    
    p$grid4 = attributes$grid4
    p$include_strat = attributes$include_strat
    p$include2 = attributes$include2
    p$nonzeroweight = attributes$nonzeroweight
    p$nobservers = attributes$nobservers
    p$strata = log(attributes$sp.mean)
    
    grid4.bird = p[p$include_strat==1 & is.na(p$include_strat)==F,]
    grid4.bird@data = grid4.bird@data %>%
      dplyr::select(layer, grid4, include_strat, include2, nonzeroweight, nobservers, strata) %>%
      mutate(strat = 1:nrow(grid4.bird@data))
    head(grid4.bird@data); nrow(grid4.bird@data)
    
    # project into Lambert conformal conic to allow plotting in OpenBUGS (with meters scale) if needed
    proj = "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=10000000 +y_0=10000000 +datum=NAD83 +units=m +no_defs"
    grid4.bird.lamb = spTransform(grid4.bird, proj)
    plot(grid4.bird.lamb, col="gray")
    head(grid4.bird.lamb); nrow(grid4.bird.lamb)
    
    # changing names of the ID variable of the polygons to grid labels so names are visible in OpenBUGS map
    # https://eco-data-science.github.io/spatial_analysis2_R/
    sapply(slot(grid4.bird.lamb, "polygons"), function(x) slot(x, "ID"))
    for (i in 1:length(slot(grid4.bird.lamb, "polygons"))){
      slot(slot(grid4.bird.lamb, "polygons")[[i]], "ID") = grid4.bird.lamb$grid4[i]
    }
    sapply(slot(grid4.bird.lamb, "polygons"), function(x) slot(x, "ID"))
  }
  
  # conditional statements telling the function what to provide as output (with code for making some)
  if(output == "grid.map"){
    return(grid4.bird.lamb)
  }
  
  if(output == "CAR.model"){
    # converting spatial data to geobugs format for CAR model
    sp2WB(grid4.bird.lamb, "grid4.bird.lamb.txt")
    
    # creating BUGS data for CAR model 
    count = route.year.grid4b$count
    strat = route.year.grid4b$strat
    obser = route.year.grid4b$obser
    firstyr = route.year.grid4b$firstyr
    year = route.year.grid4b$year
    
    ncounts = c(nrow(route.year.grid4b))
    nstrata = nrow(grid4.bird.lamb@data)
    ymin = 1
    ymax = 2017-1965
    
    nonzeroweight = grid4.bird.lamb@data$nonzeroweight
    strata = grid4.bird.lamb@data$strata
    nobservers = grid4.bird.lamb@data$nobservers
    
    require(R2WinBUGS)
    bugs.data(list("count"=count, "strat"=strat, "obser"=obser, "firstyr"=firstyr,
                   "year"=year,"ncounts"=ncounts,"nstrata"=nstrata,"ymin"=ymin,
                   "ymax"=ymax,"nonzeroweight"=nonzeroweight,"nobservers"=nobservers, 
                   "strata"=strata), 
              getwd(),12,data.file="bugs.data.txt")
    
    # note: additional "adj", "num", & "sumNumNeigh" are calculated in OpenBUGs from map file
  }
  
  if(output == "sync.dist"){ 
    return(grid4.bird.cordist)
  }
  
  if(output == "grid.year"){ 
    return(grid4.year.full %>% filter(year > yearstart))
  }
  
  if(output == "route.year"){ 
    return(brw.route.year)
  }
  
}

