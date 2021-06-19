## ---- message = FALSE, warning = FALSE, results = 'hide'-----------------------------------------------------------------
pkgs <- c("sf", "mapview", "spdep", "spatialreg", "tmap", "GWmodel", "viridisLite") # note: load spdep first, then spatialreg
lapply(pkgs, require, character.only = TRUE)



## ------------------------------------------------------------------------------------------------------------------------
sessionInfo()



## ------------------------------------------------------------------------------------------------------------------------
load("_data/msoa_spatial.RData")



## ------------------------------------------------------------------------------------------------------------------------
# Contiguity (Queens) neighbours weights
queens.nb <- poly2nb(msoa.spdf, 
                     queen = TRUE, 
                     snap = 1) # we consider points in 1m distance as 'touching'
summary(queens.nb)

# Lets plot that
plot(st_geometry(msoa.spdf), border = "grey60")
plot(queens.nb, st_centroid(st_geometry(msoa.spdf)), 
     add = TRUE, pch = 19, cex = 0.6)

# We can also transform this into a matrix W
W <- nb2mat(queens.nb, style = "B")
print(W[1:10, 1:10])



## ------------------------------------------------------------------------------------------------------------------------
# Crease centroids
coords <- st_geometry(st_centroid(msoa.spdf))

# Neighbours within 3km distance
dist_3.nb <- dnearneigh(coords, d1 = 0, d2 = 3000)
summary(dist_3.nb)

# Lets plot that
plot(st_geometry(msoa.spdf), border = "grey60")
plot(dist_3.nb, coords, 
     add = TRUE, pch = 19, cex = 0.6)



## ------------------------------------------------------------------------------------------------------------------------
queens.lw <- nb2listw(queens.nb,
                      style = "W") # W ist row-normalization
summary(queens.lw)


## ------------------------------------------------------------------------------------------------------------------------
idw.lw <- nb2listwdist(dist_3.nb,
                       x = coords, # needed for idw
                       type = "idw", # inverse distance weighting
                       alpha = 1, # the decay parameter for distance weighting
                       style = "minmax") # for eigenvalue normalization
summary(idw.lw)



## ------------------------------------------------------------------------------------------------------------------------
mp1 <- tm_shape(msoa.spdf) +
  tm_fill(col = "Value", 
          #style = "cont",
          style = "fisher", n = 8,
          title = "Median", 
          palette = viridis(n = 8, direction = -1, option = "C"),
          legend.hist = TRUE) +
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            #legend.position = c("right", "bottom"),
            legend.outside = TRUE,
            main.title = "House values 2017", 
            main.title.position = "center",
            title.snap.to.legend = TRUE) 

mp1 



## ------------------------------------------------------------------------------------------------------------------------
# Global Morans I test of housing values based on contiguity weights
moran.test(msoa.spdf$Value, listw = queens.lw, alternative = "two.sided")

# Global Morans I test of housing values based on idw
moran.test(msoa.spdf$Value, listw = idw.lw, alternative = "two.sided")



## ------------------------------------------------------------------------------------------------------------------------
hv_1.sar <- lagsarlm(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + log(traffic),  
                     data = msoa.spdf, 
                     listw = queens.lw,
                     Durbin = FALSE) # we could here extend to SDM
summary(hv_1.sar)



## ------------------------------------------------------------------------------------------------------------------------
hv_1.sem <- errorsarlm(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic,  
                     data = msoa.spdf, 
                     listw = queens.lw,
                     Durbin = FALSE) # we could here extend to SDEM
summary(hv_1.sem)



## ------------------------------------------------------------------------------------------------------------------------
hv_1.slx <- lmSLX(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic,  
                  data = msoa.spdf, 
                  listw = queens.lw, 
                  Durbin = TRUE) # use a formula to lag only specific covariates
summary(hv_1.slx)



## ------------------------------------------------------------------------------------------------------------------------
# Loop through vars and create lagged variables
msoa.spdf$log_POPDEN <- log(msoa.spdf$POPDEN)
vars <- c("Value", "log_POPDEN", "canopy_per", "pubs_count", "traffic")
for(v in vars){
  msoa.spdf[, paste0("w.", v)] <- lag.listw(queens.lw,
                                            var = st_drop_geometry(msoa.spdf)[, v])
}

# Alternatively:
w_vars <- create_WX(st_drop_geometry(msoa.spdf[, vars]),
                    listw = queens.lw,
                    prefix = "w")

head(w_vars)



## ------------------------------------------------------------------------------------------------------------------------
hv_1.lm <- lm (log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic +
                 w.log_POPDEN + w.canopy_per + w.pubs_count + w.traffic,
               data = msoa.spdf)
summary(hv_1.lm)



## ------------------------------------------------------------------------------------------------------------------------
hv_1.sar.imp <- impacts(hv_1.sar, listw = queens.lw, R = 300)
summary(hv_1.sar.imp, zstats = TRUE, short = TRUE)

# Alternative with traces (better for large W)
W <- as(queens.lw, "CsparseMatrix")
trMatc <- trW(W, type="mult")
hv_1.sar.imp2 <- impacts(hv_1.sar, tr = trMatc, R = 300, Q = 10)
summary(hv_1.sar.imp2, zstats = TRUE, short = TRUE)



## ------------------------------------------------------------------------------------------------------------------------
print(impacts(hv_1.slx, listw = queens.lw))



## ------------------------------------------------------------------------------------------------------------------------
# Search for the optimal bandwidth 
set.seed(123)
hv_1.bw <- bw.gwr(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic,
                  data = as_Spatial(msoa.spdf),
                  kernel = "boxcar",
                  adaptive = TRUE) 
hv_1.bw


### GWR 
hv_1.gwr <- gwr.robust(log(Value) ~ log(POPDEN) + canopy_per + pubs_count + traffic,
                      data = as_Spatial(msoa.spdf), 
                      kernel = "boxcar", 
                      adaptive = TRUE, 
                      bw = hv_1.bw, 
                      longlat = FALSE)
print(hv_1.gwr)


## ------------------------------------------------------------------------------------------------------------------------
# Spatial object
gwr.spdf <- st_as_sf(hv_1.gwr$SDF)
gwr.spdf <- st_make_valid(gwr.spdf)

# Map
tmap_mode("view")

mp2 <- tm_shape(gwr.spdf) +
  tm_fill(col = "canopy_per", 
          style = "hclust", n = 8,
          title = "Coefficient", 
          palette = inferno(n = 8, direction = 1),
          legend.hist = TRUE) +
  tm_borders(col = "black", lwd = 1) +
  tm_layout(legend.frame = TRUE, legend.bg.color = TRUE,
            #legend.position = c("right", "bottom"),
            legend.outside = TRUE,
            main.title = "Coefficient of tree cover", 
            main.title.position = "center",
            title.snap.to.legend = TRUE) 

mp2 


