### THESE WERE REMOVED FROM BBS_GRID4_SYNCHRONY SCRIPT FUNCTION BUT KEEPING HERE FOR FUTURE REFERENCE

# use this to make a shapefile of all BBS routes if desired
#brw.route.spatial = sp::SpatialPointsDataFrame(brw.route[,c("lon","lat")],brw.route)
#proj4string(brw.route.spatial) = CRS("+init=epsg:4326")
#plot(brw.route.spatial)
# writeOGR(brw.route.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_routes1", driver="ESRI Shapefile")

# use this to make a point shapefile of centroids for all grid cells included based on criteria in Michel et al. (Ecography)
# this is to compare with the map with exclusions based on # of zeros below
#grid4.summary.spatial = sp::SpatialPointsDataFrame(grid4.summary[grid4.summary$include_strat==1,c("grid4_lon","grid4_lat")],grid4.summary[grid4.summary$include_strat==1,])
#proj4string(grid4.summary.spatial) = CRS("+init=epsg:4326")
#plot(grid4.summary.spatial); nrow(grid4.summary.spatial@data)
# writeOGR(grid4.summary.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "bbs_grid4", driver="ESRI Shapefile")

# plot or make shapefile of grid cell centroids included based on all criteria (if needed)
#grid4.grsp.spatial = SpatialPointsDataFrame(grid4.summary2[,c("grid4_lon","grid4_lat")],grid4.summary2)
#proj4string(grid4.grsp.spatial) = CRS("+init=epsg:4326")
#plot(grid4.grsp.spatial); nrow(grid4.grsp.spatial@data)
# writeOGR(grid.grsp.spatial, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid", driver="ESRI Shapefile")

# plot or save shapefile of included grid cells if you want
plot(grid4.bird, col="gray")
#writeOGR(grid4.grsp, "C:\\Users\\Mike\\Documents\\Research\\Chapter 3 - Synchrony\\GLBsynchrony\\data", "grsp_grid4_poly", driver="ESRI Shapefile")
