#' ---
#' title: "Practical exercises: solutions"
#' author: "Tobias Rüttenauer"
#' date: "June 19, 2021"
#' output_dir: docs
#' output: 
#'   html_document:
#'     theme: flatly
#'     highlight: haddock
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'       smooth_scroll: false
#'     toc_depth: 2
#' theme: united
#' bibliography: sicss-spatial.bib
#' link-citations: yes
#' ---
#' 
#' ### Required packages
#' 
## ---- message = FALSE, warning = FALSE, results = 'hide'------------------------------------------------------------------------------------------------------------------
pkgs <- c("sf", "mapview", "spdep", "spatialreg", "tmap", "tmaptools", 
          "gstat", "randomForest", "nomisr") # note: load spdep first, then spatialreg
lapply(pkgs, require, character.only = TRUE)


#' 
#' ### Session info
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()


#'   
#' ### Load spatial data
#' 
#' See previous file.
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("_data/msoa_spatial.RData")


#' 
#' # Exercise 1
#' 
#' Add a control for the distance to Boris' home (10 Downing St, London SW1A 2AB, UK). 
#' 
#'     * You could either find the coordinates manually or you use `tmaptools` function `geocode_OSM()` using OpenStreetMaps. 
#'     
#'     * There are also Google Maps APIs like `ggmap` or `mapsapi` but they require registration.
#'     
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Geocode an address using tmaptools
adr <- "10 Downing St, London SW1A 2AB, UK"
boris.spdf <- geocode_OSM(adr, as.sf = TRUE)
mapview(boris.spdf)

# Transform into same projection
boris.spdf <- st_transform(boris.spdf, crs = st_crs(msoa.spdf))

# Copute distances betweens msoas and Boris 
msoa.spdf$dist_boris <- st_distance(msoa.spdf, boris.spdf)

# Contiguity (Queens) neighbours weights
queens.nb <- poly2nb(msoa.spdf, 
                     queen = TRUE, 
                     snap = 1)
queens.lw <- nb2listw(queens.nb,
                      style = "W")

# Add to an SLX model (using queens.lw from previous script)
hv_2.slx <- lmSLX(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic + dist_boris,  
                  data = msoa.spdf, 
                  listw = queens.lw, 
                  Durbin = as.formula( ~ log(POPDEN) + canopy_per + 
                                         pubs_count + traffic)) # we omit dist from the lagged vars
summary(hv_2.slx)



#' 
#' Compared to previous results, adding the distance to Downing street, changed the effect of the lagged population density. So, our previous explanation based on centrality seems to make sense.
#'     
#' # Exercise 2
#' 
#' How could we improve the interpolation of traffic counts based on `idw()`? 
#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load traffic data 
load("_data/traffic.RData")

# Set up regular grid over London with 1km cell size
london.sp <- st_union(st_geometry(msoa.spdf))
grid <- st_make_grid(london.sp, cellsize = 1000)
mapview(grid)

# Reduce to London bounds
grid <- grid[london.sp]

all_motor.idw <- idw(all_motor_vehicles ~ 1,
                     locations = traffic.spdf,
                     newdata = grid,
                     idp = 0.9, # Set decay lower
                     maxdist = 2000, # Define sharp cutoff 
                     nmin = 5, # but use at least 5 points
                     force = TRUE) # even if there are less within maxdist
mapview(all_motor.idw[, "var1.pred"])


#' 
#' 
#' # Exercise 3    
#' 
#' Can you think of a way of testing if higher order neighbours add additional information to the SLX regression model?
#' 
#'     * The function `nblag()` is helpful.
#'     
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create neighbours of orders 1 to 3
queens.lag <- nblag(queens.nb, maxlag = 3)

# Create listwise of 1st, 2nd and 3rd order neighbours
queens.lw1 <- nb2listw(queens.lag[[1]], style = "W")
queens.lw2 <- nb2listw(queens.lag[[2]], style = "W")
queens.lw3 <- nb2listw(queens.lag[[3]], style = "W")

### Create lagged variables for different orders of neighbours
msoa.spdf$log_POPDEN <- log(msoa.spdf$POPDEN)
vars <- c("Value", "log_POPDEN", "canopy_per", "pubs_count", "traffic")
# lag1
w_vars <- create_WX(st_drop_geometry(msoa.spdf[, vars]),
                    listw = queens.lw1,
                    prefix = "w1")
msoa.spdf <- cbind(msoa.spdf, w_vars)
# lag2
w_vars <- create_WX(st_drop_geometry(msoa.spdf[, vars]),
                    listw = queens.lw2,
                    prefix = "w2")
msoa.spdf <- cbind(msoa.spdf, w_vars)
# lag3
w_vars <- create_WX(st_drop_geometry(msoa.spdf[, vars]),
                    listw = queens.lw3,
                    prefix = "w3")
msoa.spdf <- cbind(msoa.spdf, w_vars)


### LM to test imporatance

lm.mod <- lm(log(Value) ~ log_POPDEN + canopy_per + pubs_count + traffic +
               w1.log_POPDEN + w1.canopy_per + w1.pubs_count + w1.traffic +
               w2.log_POPDEN + w2.canopy_per + w2.pubs_count + w2.traffic +
               w3.log_POPDEN + w3.canopy_per + w3.pubs_count + w3.traffic, 
             data = st_drop_geometry(msoa.spdf))
summary(lm.mod)


### Random forest to test importance

# Train
rf.mod <- randomForest(log(Value) ~ log_POPDEN + canopy_per + pubs_count + traffic +
                         w1.log_POPDEN + w1.canopy_per + w1.pubs_count + w1.traffic +
                         w2.log_POPDEN + w2.canopy_per + w2.pubs_count + w2.traffic +
                         w3.log_POPDEN + w3.canopy_per + w3.pubs_count + w3.traffic, 
                       data = st_drop_geometry(msoa.spdf), 
                       mtry = 2, 
                       ntree = 1000,
                       importance = TRUE)

# Inspect the mechanics of the model
importance(rf.mod)
varImpPlot(rf.mod)


#'     
#'     
#' # Exercise 4    
#' 
#' Add the amount of particulate matter (PM10, e.g. 2017) from [Defra](https://uk-air.defra.gov.uk/data/pcm-data) and check if pollution influences the house values 
#' 
#'     * Note the following important sentence: "The coordinate system is OSGB and the coordinates represent the centre of each 1x1km cell".
#' 
#'     * `st_buffer()` with the options `nQuadSegs  = 1, endCapStyle = 'SQUARE'` provides an easy way to get points into grids
#'     
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Download
pm10.link <- "https://uk-air.defra.gov.uk/datastore/pcm/mappm102017g.csv"
pm10.df <- read.csv(pm10.link, skip = 5, na.strings = "MISSING")

# Convert to sf point data
pm10.spdf <- st_as_sf(pm10.df, 
                      coords = c("x", "y"), # Order is important
                      crs = 27700) # EPSG number of CRS

plot(pm10.spdf[1:100, "pm102017g"])

# Add qudratic buffer with 500m distance (1km/2) to get grid
pm10.spdf <- st_buffer(pm10.spdf, dist = 500, nQuadSegs  = 1, 
                       endCapStyle = 'SQUARE')

plot(pm10.spdf[1:100, "pm102017g"])

# Subset to London
pm10.spdf <- pm10.spdf[msoa.spdf, ]

plot(pm10.spdf[, "pm102017g"])


### Add to msoa data

# Calculate intersection
smoa_pm10.int <- st_intersection(msoa.spdf, pm10.spdf)

# average per MSOA
smoa_pm10.int <- aggregate(list(pm10 = smoa_pm10.int$pm102017g),
                         by = list(MSOA11CD = smoa_pm10.int$MSOA11CD),
                         mean)

# Merge back
msoa.spdf <- merge(msoa.spdf, smoa_pm10.int, by = "MSOA11CD", all.x = TRUE)


### Add to our SLX regression

hv_3.slx <- lmSLX(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic + 
                    dist_boris + pm10,  
                  data = msoa.spdf, 
                  listw = queens.lw, 
                  Durbin = as.formula( ~ log(POPDEN) + canopy_per + 
                                         pubs_count + traffic + pm10)) # we omit dist from the lagged vars
summary(hv_3.slx)

#'     
#' 
#' # Exercise 5
#' 
#' Add some more demographic variables from the Census 2011.
#' 
#'     * The `nomisr` package provides an API to nomis. See the [Vignette](https://cran.r-project.org/web/packages/nomisr/vignettes/introduction.html). Make sure to restrict your request to London only (Guest users are limited to 25,000 rows per query).
#'     
#'     * You can browse the available data online, such as the Census 2011 [key statistics](https://www.nomisweb.co.uk/census/2011/key_statistics_uk).
#'     
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### For larger request, register and set key
# Sys.setenv(NOMIS_API_KEY = "XXX")
# nomis_api_key(check_env = TRUE)

x <- nomis_data_info()

# Get London ids
london_ids <- msoa.spdf$MSOA11CD

### Get key statistics ids
# select requires tables (https://www.nomisweb.co.uk/census/2011/key_statistics_uk)

# Get internal ids
stats <- c("KS401UK", "KS402UK")
oo <- which(grepl(paste(stats, collapse = "|"), x$name.value))
ksids <- x$id[oo]
ksids # This are the inernal ids


### Lets look at meta information
q <- nomis_overview(ksids[1])
head(q)
a <- nomis_get_metadata(id = ksids[1], concept = "GEOGRAPHY", type = "type")
a # TYPE297 is MSOA level

b <- nomis_get_metadata(id = ksids[1], concept = "MEASURES", type = "TYPE297")
b # 20100 is the measure of absolute numbers


### Query data in loop
for(i in ksids){
  
  # Check if table is urban-rural devided
  nd <- nomis_get_metadata(id = i)
  if("RURAL_URBAN" %in% nd$conceptref){
    UR <- TRUE
  }else{
    UR <- FALSE
  }
  
  # Query data
  if(UR == TRUE){
    ks_tmp <- nomis_get_data(id = i, time = "2011", 
                             geography = london_ids, # replace with "TYPE297" for all
                             measures = 20100, RURAL_URBAN = 0)
  }else{
    ks_tmp <- nomis_get_data(id = i, time = "2011", 
                             geography = london_ids, # replace with "TYPE297" for all
                             measures = 20100)
  }
  

  
  # Make lower case names
  names(ks_tmp) <- tolower(names(ks_tmp))
  names(ks_tmp)[names(ks_tmp) == "geography_code"] <- "msoa"
  names(ks_tmp)[names(ks_tmp) == "geography_name"] <- "name"
  
  # replace weird cell codes
  onlynum <- which(grepl("^[[:digit:]]+$", ks_tmp$cell_code))
  if(length(onlynum) != 0){
    code <- substr(ks_tmp$cell_code[-onlynum][1], 1, 7)
    ks_tmp$cell_code[onlynum] <- paste0(code, "_", ks_tmp$cell_code[onlynum])
  }
  
  # save codebook
  ks_cb <- unique(ks_tmp[, c("date", "cell_type", "cell", "cell_code", "cell_name")])
  
  ### Reshape
  ks_res <- tidyr::pivot_wider(ks_tmp, id_cols = c("msoa", "name"),
                                names_from = "cell_code",
                                values_from = "obs_value")
  
  ### Merge
  if(i == ksids[1]){
    census_keystat.df <- ks_res
    census_keystat_cb.df <- ks_cb
  }else{
    census_keystat.df <- merge(census_keystat.df, ks_res, by = c("msoa", "name"), all = TRUE)
    census_keystat_cb.df <- rbind(census_keystat_cb.df, ks_cb)
  }
  
}

# Descirption are saved in the codebook
census_keystat_cb.df

### Merge with MSOA
msoa.spdf <- merge(msoa.spdf, census_keystat.df, 
                   by.x = "MSOA11CD", by.y = "msoa", all.x = TRUE)


#' 
#'     
#' 
