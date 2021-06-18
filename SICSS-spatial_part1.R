## ---- message = FALSE, warning = FALSE, results = 'hide'-----------------------------------------------------------------
pkgs <- c("sf", "gstat", "mapview", "nngeo", "rnaturalearth", "dplyr") 
lapply(pkgs, require, character.only = TRUE)



## ------------------------------------------------------------------------------------------------------------------------
sessionInfo()



## ------------------------------------------------------------------------------------------------------------------------
# Coordinate pairs of two locations
coords1 <- c(51.752595, -1.262801)
coords2 <- c(51.753237, -1.253904)
coords <- rbind(coords1, coords2)

# Conventional data frame
nuffield.df <- data.frame(name = c("Nuffield College", "Radcliffe Camera"),
                          address = c("New Road", "Radcliffe Sq"),
                          lat = coords[,1], lon = coords[,2])

head(nuffield.df)

# Combine to spatial data frame
nuffield.spdf <- st_as_sf(nuffield.df, 
                          coords = c("lon", "lat"), # Order is important
                          crs = 4326) # EPSG number of CRS

# Map
mapview(nuffield.spdf, zcol = "name")



## ------------------------------------------------------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
st_crs(world)

# Extract a country and plot in current CRS (WGS84)
ger.spdf <- world[world$name == "Germany", ]
plot(st_geometry(ger.spdf))

# Now, lets transform Germany into a CRS optimized for Iceland
ger_rep.spdf <- st_transform(ger.spdf, crs = 5325)
plot(st_geometry(ger_rep.spdf))



## ------------------------------------------------------------------------------------------------------------------------
# Create subdir 
dn <- "_data"
ifelse(dir.exists(dn), "Exists", dir.create(dn))

# Download zip file and unzip
tmpf <- tempfile()
boundary.link <- "https://data.london.gov.uk/download/statistical-gis-boundary-files-london/9ba8c833-6370-4b11-abdc-314aa020d5e0/statistical-gis-boundaries-london.zip"
download.file(boundary.link, tmpf)
unzip(zipfile = tmpf, exdir = paste0(dn))
unlink(tmpf)

# We only need the MSOA layer for now
msoa.spdf <- st_read(dsn = paste0(dn, "/statistical-gis-boundaries-london/ESRI"),
                     layer = "MSOA_2011_London_gen_MHW" # Note: no file ending
                     )
head(msoa.spdf)



## ------------------------------------------------------------------------------------------------------------------------
mapview(msoa.spdf[, "POPDEN"])



## ------------------------------------------------------------------------------------------------------------------------

# Download
hp.link <- "https://data.london.gov.uk/download/average-house-prices/bdf8eee7-41e1-4d24-90ce-93fe5cf040ae/land-registry-house-prices-MSOA.csv"
hp.df <- read.csv(hp.link)
hp.df <- hp.df[which(hp.df$Measure == "Median" &
                       hp.df$Year == "Year ending Sep 2017"), ]
hp.df$Value <- as.numeric(hp.df$Value)

# Merge (as with conventional df)
msoa.spdf <- merge(msoa.spdf, hp.df,
                   by.x = "MSOA11CD", by.y = "Code",
                   all.x = TRUE, all.y = FALSE)

mapview(msoa.spdf[, "Value"])



## ------------------------------------------------------------------------------------------------------------------------
# Download zip shapefile
tmpf <- tempfile()
trees.link <- "https://data.london.gov.uk/download/curio-canopy/4fd54ef7-195f-43dc-a0d1-24e96e876f6c/shp-hexagon-files.zip"
download.file(trees.link, tmpf)
unzip(zipfile = tmpf, exdir = paste0(dn))
unlink(tmpf)
trees.spdf <- st_read(dsn = paste0(dn, "/shp-hexagon-files"),
                      layer = "gla-canopy-hex")

mapview(trees.spdf[, "canopy_per"])



## ------------------------------------------------------------------------------------------------------------------------
culture.link <- "https://data.london.gov.uk/download/cultural-infrastructure-map/bf7822ed-e90a-4773-abef-dda6f6b40654/CulturalInfrastructureMap.gpkg"

# This time, we have Geopackage format (gpkg)
tmpf <- tempfile(fileext = ".gpkg")
download.file(culture.link, destfile = tmpf, mode = "wb")

# And we only load the layer containing pubs
st_layers(tmpf)
pubs.spdf <- st_read(dsn = tmpf, layer = "Pubs")
unlink(tmpf)

mapview(st_geometry(pubs.spdf))



## ------------------------------------------------------------------------------------------------------------------------
st_crs(msoa.spdf) == st_crs(trees.spdf)
st_crs(trees.spdf) == st_crs(pubs.spdf)

# MSOA in different crs --> transform
msoa.spdf <- st_transform(msoa.spdf, crs = st_crs(trees.spdf))



## ------------------------------------------------------------------------------------------------------------------------
# Subset msoa and combine to one unit
westminster.spdf <- msoa.spdf[msoa.spdf$LAD11NM == "Westminster",]
westminster.spdf <- st_union(westminster.spdf)

# Subset to points in this area
sub1.spdf <- pubs.spdf[westminster.spdf, ] # or:
sub1.spdf <- st_filter(pubs.spdf, westminster.spdf)
mapview(sub1.spdf)



## ------------------------------------------------------------------------------------------------------------------------
# Subset to points not in this area
sub2.spdf <- pubs.spdf[westminster.spdf, , op = st_disjoint] # or:
sub2.spdf <- st_filter(pubs.spdf, westminster.spdf, .predicate = st_disjoint)
mapview(list(sub1.spdf, sub2.spdf), col.regions = list("red", "blue"))



## ------------------------------------------------------------------------------------------------------------------------
# Assign MSOA to each point
pubs_msoa.join <- st_join(pubs.spdf, msoa.spdf, join = st_within)

# Count N by MSOA code (drop geometry to speed up)
pubs_msoa.join <- dplyr::count(st_drop_geometry(pubs_msoa.join),
                               MSOA11CD = pubs_msoa.join$MSOA11CD,
                               name = "pubs_count")
sum(pubs_msoa.join$pubs_count)

# Merge and replace NAs with zero (no matches, no pubs)
msoa.spdf <- merge(msoa.spdf, pubs_msoa.join,
                   by = "MSOA11CD", all.x = TRUE)
msoa.spdf$pubs_count[is.na(msoa.spdf$pubs_count)] <- 0



## ------------------------------------------------------------------------------------------------------------------------
hist(trees.spdf$canopy_per)

# Select areas with at least 50% tree coverage
trees_high.spdf <- trees.spdf[trees.spdf$canopy_per >= 50, ]
mapview(trees_high.spdf[, "canopy_per"])

# Use geometric centroid of each MSOA
cent.sp <- st_centroid(msoa.spdf[, "MSOA11CD"])

# Get K nearest neighbour with distance
knb.dist <- st_nn(cent.sp, trees_high.spdf,
                  k = 1, returnDist = TRUE,
                  progress = FALSE)
msoa.spdf$dist_trees50 <- unlist(knb.dist$dist)
summary(msoa.spdf$dist_trees50)



## ------------------------------------------------------------------------------------------------------------------------
# Create buffer (1km radius)
cent.buf <- st_buffer(cent.sp, dist = 1000)
mapview(cent.buf)

# Calculate intersection between buffers and tree-cover hexagons
buf_hex.int <- st_intersection(cent.buf, trees.spdf)
dim(buf_hex.int)

# We could also calculate the area of overlap for each pair (to calculate weighted averages)
buf_hex.int$area <- as.numeric(st_area(buf_hex.int))

# Or we just use the simple average per each MSOA
buf_hex.int <- aggregate(list(canopy_per = buf_hex.int$canopy_per),
                         by = list(MSOA11CD = buf_hex.int$MSOA11CD),
                         mean)

# Merge back to spatial data.frame
msoa.spdf <- merge(msoa.spdf, buf_hex.int, by = "MSOA11CD", all.x = TRUE)



## ------------------------------------------------------------------------------------------------------------------------
# Download
traffic.link <- "https://dft-statistics.s3.amazonaws.com/road-traffic/downloads/aadf/region_id/dft_aadf_region_id_6.csv"
traffic.df <- read.csv(traffic.link)
traffic.df <- traffic.df[which(traffic.df$year == 2017 &
                       traffic.df$road_type == "Major"), ]
names(traffic.df)

# Transform to spatial data
traffic.spdf <- st_as_sf(traffic.df,
                         coords = c("longitude", "latitude"),
                         crs = 4326)

# Transform into common crs
traffic.spdf <- st_transform(traffic.spdf,
                             crs = st_crs(msoa.spdf))

# Map
mapview(traffic.spdf[, "all_motor_vehicles"])



## ------------------------------------------------------------------------------------------------------------------------
# Set up regular grid over London with 1km cell size
london.sp <- st_union(st_geometry(msoa.spdf))
grid <- st_make_grid(london.sp, cellsize = 1000)
mapview(grid)

# Reduce to London bounds
grid <- grid[london.sp]
mapview(grid)




## ------------------------------------------------------------------------------------------------------------------------
# IDW interpolation
all_motor.idw <- idw(all_motor_vehicles ~ 1,
                     locations = traffic.spdf,
                     newdata = grid,
                     idp = 2) # power of distance decay
mapview(all_motor.idw[, "var1.pred"])


## ------------------------------------------------------------------------------------------------------------------------
# Variogram
all_motor.var <- variogram(all_motor_vehicles ~ 1,
                          traffic.spdf)
plot(all_motor.var)

# Fit variogram
all_motor.varfit <- fit.variogram(all_motor.var,
                                  vgm(c("Nug", "Exp")), # use vgm() to see options
                                  fit.kappa = TRUE, fit.sills = TRUE, fit.ranges = TRUE)

plot(all_motor.var, all_motor.varfit)

# Parameters
all_motor.varfit


### Ordinary Kriging (ATTENTION: This can take some time!)

all_motor.kg <- krige(all_motor_vehicles ~ 1,
                      locations = traffic.spdf,
                      newdata = grid,
                      model = all_motor.varfit)

# Look at results
mapview(all_motor.kg[, "var1.pred"])



## ------------------------------------------------------------------------------------------------------------------------
# Calculate intersection
smoa_grid.int <- st_intersection(msoa.spdf, all_motor.kg)

# average per MSOA
smoa_grid.int <- aggregate(list(traffic = smoa_grid.int$var1.pred),
                         by = list(MSOA11CD = smoa_grid.int$MSOA11CD),
                         mean)

# Merge back
msoa.spdf <- merge(msoa.spdf, smoa_grid.int, by = "MSOA11CD", all.x = TRUE)



## ------------------------------------------------------------------------------------------------------------------------
# Save
save(msoa.spdf, file = "_data/msoa_spatial.RData")


