## ----setup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, 
                      warning = FALSE, cache = TRUE)


## ----------------------------------------------------------------------
library(ggplot2)
library(sf)


## ----loaddata----------------------------------------------------------
load("data/elk_processed.rda")
is(elk_gps)
is(elk_sf)


## ----subsetdata--------------------------------------------------------
str(elk_gps)
elk_gps <- elk_gps |> subset(id %in% levels(id)[1:9])
elk_sf <- elk_sf |> subset(id %in% levels(id)[1:9])


## ----------------------------------------------------------------------
E4049 <- subset(elk_gps, id == "4049")
with(E4049, plot(lon, lat))


## ----------------------------------------------------------------------
with(E4049, plot(lon, lat, type = "o"))


## ----------------------------------------------------------------------
long2UTM <- function(long) (floor((long + 180)/6) %% 60) + 1


## ----------------------------------------------------------------------
long2UTM(-115)


## ----E4049_xyplot------------------------------------------------------
E4049.xy <-  elk_sf |> subset(id == 4049) |> st_transform(32611) |> st_coordinates()
plot(E4049.xy, asp = 1, type = "o")


## ----addXY-------------------------------------------------------------
require(sf)
elk_gps <-  elk_gps |> data.frame(elk_sf |> st_transform(32611) |> st_coordinates())
str(elk_gps)


## ----scantrack1, echo = -1---------------------------------------------
par(mar = c(0,4,0,0), oma = c(4,0,5,2), xpd=NA)
layout(rbind(c(1,2), c(1,3)))
with(E4049, {
  plot(lon, lat, asp = 1, type="o", ylab="Latitude", xlab="Longitude")
  plot(datetime, lon, type="o", xaxt="n", ylab="Longitude", xlab="")
  plot(datetime, lat, type="o", ylab="Latitude", xlab="Datetime")
  title(paste("ID", id[1]), outer = TRUE)
})



## ----scan_track function-----------------------------------------------
scan_track <- function(dataframe, x = "lon", y = "lat",
                             time = "datetime", id = "id", ...){
  par(mar = c(0,4,0,0), oma = c(4,0,5,2), xpd=NA)
  layout(rbind(c(1,2), c(1,3)))
  with(dataframe, {
    plot(get(x), get(y), asp = 1, type="o", ylab=y, xlab=x, ...)
    plot(get(time), get(x), type="o", xaxt="n", ylab=y, xlab="", ...)
    plot(get(time), get(y), type="o", ylab=y, xlab=time, ...)
    title(paste("ID", id[1]), outer = TRUE)
  })
}


## ----------------------------------------------------------------------
scan_track(elk_gps |> subset(id == id[1]))


## ----------------------------------------------------------------------
myelk <- elk_gps |> subset(id == id[1])
scan_track(myelk, x = "X", y = "Y", 
                 col = topo.colors(nrow(myelk)))


## ----ScanTrack5Elk, fig.height = 3-------------------------------------
require(plyr)
elk_gps |> subset(id %in% levels(id)[1:5]) |> d_ply("id", scan_track)


## ----------------------------------------------------------------------
library(ggplot2)


## ----ggplot------------------------------------------------------------
ggplot(data = E4049, aes(x = lon, y = lat)) +
   geom_path(size = 0.5) +
   geom_point(aes(color = datetime)) + theme_classic()


## ----facetted_elk, fig.height=9, fig.width=9---------------------------
ggplot(data = elk_gps, aes(x = lon, y = lat)) +
  geom_path(size = 0.5, color = "darkgrey") +
  geom_point() +
  theme_classic() +
  facet_wrap(~id, scale="free", ncol=3)


## ----all_the_elk_on_one_ggplot, fig.height=9, fig.width=9--------------
ggplot(data = elk_gps, aes(x = lon, y = lat, col = id, group =id)) +
  geom_path(size = 0.5, color = "darkgrey") +
  geom_point() +
  theme_classic() 


## ----fig.height=9, fig.width=9-----------------------------------------
ggplot(data = elk_gps, aes(x = datetime, y = lat)) +
  geom_path(size = 0.5) +
  xlab("DateTime") + ylab("Latitude") +
  theme_classic() +
  facet_wrap(~id, scale="free", ncol = 2)


## ----------------------------------------------------------------------
plot(elk_sf[,"id"], pch = 19)


## ----------------------------------------------------------------------
require(basemaps)


## ----------------------------------------------------------------------
elk_bbox <- st_bbox(elk_sf)
elk_bbox


## ----------------------------------------------------------------------
elk_bbox <- st_bbox(elk_sf) + c(-.25,-.25,.25,.25)


## ----------------------------------------------------------------------
elk_basemap <- basemap_raster(ext = elk_bbox, 
                              map_service = "esri", 
                              map_type = "world_imagery",
                              zoom = 10)


## ----basemap, fig.width= 4, fig.height = 8-----------------------------
require(raster)
plotRGB(elk_basemap)


## ----fig.width= 4, fig.height = 8--------------------------------------
plotRGB(elk_basemap)
plot(elk_sf[,"id"] |> st_transform(crs(elk_basemap)), 
     add = TRUE, pch = 19, cex = 0.8)


## ----------------------------------------------------------------------
require(ggspatial)
ggplot() +
  layer_spatial(elk_basemap) +  # Plot the raster basemap
  geom_sf(data = elk_sf, aes(col = id)) +  # Add elk data
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()


## ----elk basemap, fig.width = 4, fig.height = 6------------------------
elk_basemap <- ggplot() +
  annotation_map_tile(type = 'osm', zoom = 10) +
  annotation_scale() +
  annotation_north_arrow(height=unit(0.5,"cm"), width=unit(0.5,"cm"), 
                         pad_y = unit(1,"cm")) +
  shadow_spatial(elk_bbox) +
  ylab("Latitude") + xlab("Longitude")  


## ----------------------------------------------------------------------
require(dplyr)
elk_tracks <- elk_sf |> 
  group_by(id) |> 
  summarize(do_union=FALSE) |> 
  st_cast("LINESTRING")


## ----------------------------------------------------------------------
elk_tracks


## ----addElkLines, fig.width = 5, fig.height = 6------------------------
elk_basemap +  
  geom_sf(data=elk_tracks, aes(col = id)) 


## ----eval = FALSE------------------------------------------------------
# jpeg(file="ElkMap.jpg", units="in", width=4, height=7,res=300)
#   elk_basemap +
#     geom_sf(data=elk_sf, aes(col = id))
# dev.off()


## ----mapview_tracks----------------------------------------------------
library(mapview)
mapview(elk_tracks, zcol="id")

