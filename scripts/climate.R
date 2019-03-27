#######################################
### Climate data: unzip, combine, filter, and summarize temperature and precip data (0.5 to 2 degree grid)
#######################################
# Temp and precip data from:
#  Matsuura, Kenji & National Center for Atmospheric Research Staff (Eds). Last modified 20 Oct 2017. 
# "The Climate Data Guide: Global (land) precipitation and temperature: Willmott & Matsuura, 
# University of Delaware." 
# Retrieved from https://climatedataguide.ucar.edu/climate-data/global-land-precipitation-and-temperature-willmott-matsuura-university-delaware.

# functions for identifying even and odd numbers (for creating grid4 ID)
is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0

# untar("data//climate//precip_2017.tar.gz") 
# untar("data//climate//air_temp_2017.tar.gz") 
# do this once to extract all precip years, then only keep what you need

# read in all PRECIP data 1966-2017, filtered to include just lower 48 USA
# code from: https://stackoverflow.com/questions/29402528/append-data-frames-together-in-a-for-loop
precip.datalist = list()
for (i in 1966:2017) {
  dat <- read.table(paste("data//climate//precip.",i,sep="")) %>% 
    filter(V1 >= -154, V1 <= -53, V2 >= 24, V2 <= 70)
  dat$year <- i  # keep track of year
  precip.datalist[[i-1965]] <- dat # add it to list
}
precip = dplyr::bind_rows(precip.datalist) %>%
  rename(lon = V1, lat = V2, jan = V3, feb = V4, mar = V5, apr = V6, may = V7, jun = V8, jul = V9) %>%
  mutate(grid4_lat = ifelse(is.odd(trunc(lat)),trunc(lat)+1,trunc(lat))) %>%
  mutate(grid4_lon = ifelse(is.odd(trunc(lon)),trunc(lon)-1,trunc(lon))) %>%
  mutate(grid4 = gsub("-","_",paste(grid4_lat,grid4_lon,sep=""))) %>%
  dplyr::select(grid4, grid4_lat, grid4_lon, year, mar, apr, may, jun, jul) %>%
  mutate(marmay = mar + apr + may, junjul = jun + jul) %>%
  group_by(grid4, grid4_lat, grid4_lon, year) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(yearlag = year+1) %>%
  mutate(grid4.year = paste(grid4,".",yearlag,sep="")) %>%
  arrange(grid4, year)

# testing that grid coding makes sense in QGIS
#precip.spatial = SpatialPointsDataFrame(precip[,c("lon","lat")],precip)
#proj4string(precip.spatial) = CRS("+init=epsg:4326")
#plot(precip.spatial)
# writeOGR(precip.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "precip.spatial", driver="ESRI Shapefile")

# read in all TEMPERATURE data 1966-2017, filtered to include just lower 48 USA
air_temp.datalist = list()
for (i in 1966:2017) {
  dat <- read.table(paste("data//climate//air_temp.",i,sep="")) %>% 
    filter(V1 >= -154, V1 <= -53, V2 >= 24, V2 <= 70)
  dat$year <- i  # keep track of year
  air_temp.datalist[[i-1965]] <- dat # add it to list
}
air_temp = dplyr::bind_rows(air_temp.datalist)  %>%
  rename(lon = V1, lat = V2, jan = V3, feb = V4, mar = V5, apr = V6, may = V7, jun = V8, jul = V9) %>%
  mutate(grid4_lat = ifelse(is.odd(trunc(lat)),trunc(lat)+1,trunc(lat))) %>%
  mutate(grid4_lon = ifelse(is.odd(trunc(lon)),trunc(lon)-1,trunc(lon))) %>%
  mutate(grid4 = gsub("-","_",paste(grid4_lat,grid4_lon,sep=""))) %>%
  dplyr::select(grid4, grid4_lat, grid4_lon, year, mar, apr, may, jun, jul) %>%
  group_by(grid4, grid4_lat, grid4_lon, year) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  mutate(yearlag = year+1) %>%
  mutate(grid4.year = paste(grid4,".",yearlag,sep="")) %>%
  arrange(grid4, year)

# try drought? 0.5 degree data are here:
# https://climatedataguide.ucar.edu/climate-data/cru-sc-pdsi-self-calibrating-pdsi-over-europe-north-america


#######################################
### Match up bird and climate data, 1992-2017
#######################################
# should run code from bbs_grid4_synchrony.R first

grsp.climate.92.17 = grsp.grid.year.92.17 %>%
  left_join(precip, by = "grid4.year") %>%
  left_join(air_temp, by = "grid4.year") %>%
  select(grid4 = grid4.x, grid4_lat = grid4_lat.x, grid4_lon = grid4_lon.x, 
         marmayprec = marmay, junjulprec = junjul, yearbird = year.x, yearclim = year.y, 
         juntemp = jun.y, jultemp = jul.y, change, AR)

grsp.climgrid.92.17 = grsp.climate.92.17 %>%
  group_by(grid4) %>%
  summarize(juntemp.cor = cor(AR, juntemp), jultemp.cor = cor(AR, jultemp),
            marmayprec.cor = cor(AR,marmayprec), junjulprec.cor = cor(AR, junjulprec),
            grid4_lat = mean(grid4_lat), grid4_lon = mean(grid4_lon)) %>%
  mutate(juntemp.corcat = cut(juntemp.cor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1))) %>%
  mutate(jultemp.corcat = cut(jultemp.cor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1))) %>%
  mutate(marmayprec.corcat = cut(marmayprec.cor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1))) %>%
  mutate(junjulprec.corcat = cut(junjulprec.cor, breaks = c(-1,0,0.1,.2,.3,.4,.5,1)))

#######################################
### Map correlation between bird and climate data, 1992-2017
#######################################

# grsp AR v. marmayprec
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = marmayprec.corcat), size = 2, alpha = .7, data = grsp.climgrid.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. Lag1(Mar-May Precipitation), 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

# grsp AR v. junjulprec
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = junjulprec.corcat), size = 2, alpha = .7, data = grsp.climgrid.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. Lag1(Jun-Jul Precipitation), 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))
#ggsave("figures/grid4_maps/GRSPsync_v_junjulprec.cor.92.17.png")

x11(10,10)
mapsync %>%
  filter(sp=="GRSP",per=="1992-2017") %>%
  rename(grid4 = grid) %>%
  left_join(grsp.climgrid.92.17, by = "grid4") %>%
  ggplot() +
  geom_point(aes(x = junjulprec.cor, y = meancor)) +
  geom_smooth(aes(x = junjulprec.cor, y = meancor),method = "lm")
#ggsave("figures/grid4_maps/GRSPsync_v_junjulprec.cor.92.17.scatter.png")

lmdata = mapsync %>%
  filter(sp=="GRSP",per=="1992-2017") %>%
  rename(grid4 = grid) %>%
  left_join(grsp.climgrid.92.17, by = "grid4")

summary(lm(meancor~junjulprec.cor, data = lmdata))

# grsp AR v. juntemp
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = juntemp.corcat), size = 2, alpha = .7, data = grsp.climgrid.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. Lag1(Avg. June Temperature), 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

mapsync %>%
  filter(sp=="GRSP",per=="1992-2017") %>%
  rename(grid4 = grid) %>%
  left_join(grsp.climgrid.92.17, by = "grid4") %>%
  ggplot() +
  geom_point(aes(x = juntemp.cor, y = meancor)) +
  geom_smooth(aes(x = juntemp.cor, y = meancor),method = "lm")

# grsp AR v. jultemp
x11(10,10)
na + # na is defined in "mapping mean synchrony" above
  geom_point(aes(x = grid4_lon, y = grid4_lat, color = jultemp.corcat), size = 2, alpha = .7, data = grsp.climgrid.92.17) + 
  scale_colour_manual(values = c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")) +
  theme_classic() +
  theme(text = element_text(size=15)) +
  #  facet_wrap(~sp, ncol=2) +
  labs(colour="Cross\ncorrelation", x="", y="", title="GRSP v. Lag1(Avg. July Temperature), 1992-2017") + 
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.title.align=0.5) +
  theme(axis.line = element_line(color = "transparent")) +
  xlim(c(-142,-59))

(plotdata = mapsync %>%
    filter(sp=="GRSP",per=="1992-2017") %>%
    rename(grid4 = grid) %>%
    left_join(grsp.climgrid.92.17, by = "grid4") %>%
    ggplot() +
    geom_point(aes(x = jultemp.cor, y = meancor))) +
  geom_smooth(aes(x = jultemp.cor, y = meancor),method = "lm")





