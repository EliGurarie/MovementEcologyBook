load("C:/Users/egurarie/teaching/movementbook/SpatialAndMovementEcologyCourse/EFB796/data/elk_gps.rda")

c("YL80", "YL91", "YL94")

myelk <- elk_gps |> subset(id == "YL94")

myelk_sf <- myelk |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
  st_transform(32611)

myelk_sf <- cbind(myelk_sf, 
                       st_coordinates(myelk_sf))

require(mapview)
mapview(myelk_sf)


# MCP

XY <- st_coordinates(myelk_sf)
chull(XY)

plot(XY, asp = 1, cex = 0.3)
XY[1,]
XY[chull(XY),]
points(XY[chull(XY),], col = 2, pch = 19)
polygon(XY[chull(XY),])

Z <- XY[,1] + 1i * XY[,2]
plot(Z, asp = 1, cex = 0.5)

# centroid

Z.centroid <- mean(Z)
points(Z.centroid, col = 2, pch = 4, lwd =3, cex = 4)

Distances <- Mod(Z - Z.centroid)
cutoff <- quantile(Distances, .95)
plot(Distances)
abline(h = cutoff, col = 3)

Z_trimmed <- Z[Distances < cutoff]
X_trimmed <- Re(Z_trimmed)
Y_trimmed <- Im(Z_trimmed)
XY_trimmed <- cbind(X_trimmed, Y_trimmed)  
plot(Z, asp = 1, cex = 0.5)
polygon(XY_trimmed[chull(XY_trimmed),])

# adehabitatHR

  install.packages("adehabitatHR")
  require(adehabitatHR)
  
  myelk_sp <- myelk_sf |> as_Spatial()
  
  mcp100 <- mcp(myelk_sp, percent = 100, unout = "km2")
  mcp95 <- mcp(myelk_sp, percent = 95, unout = "km2")
  
  
  st_as_sf(mcp100) # convert back to simple feature
  
 mapview(myelk_sf) + 
   mapview(st_as_sf(mcp100)) + 
   mapview(st_as_sf(mcp95))


# Kernel density estimate
 
 ?kernelUD
myelk_kernel <- kernelUD(myelk_sp, grid = 200)

require(raster)
mapview(raster(myelk_kernel))



plot(raster(myelk_kernel), xlim = c(590e3, 610e3), 
     ylim = c(573e4, 574.5e4))
contour(raster(myelk_kernel), add = TRUE)


kernel95_polygon <- 
  getverticeshr(myelk_kernel, percent = 95) |> st_as_sf()

mapview(raster(myelk_kernel)) + 
        mapview(kernel95_polygon)

st_area(kernel95_polygon)/1e6





# LoCoH

require("amt")

myelk_track <- make_track(myelk_sf, 
                          .x = X, .y = Y, .t = datetime,
                          crs = st_crs(myelk_sf))

myelk_locoh <- hr_locoh(myelk_track, type = "k", 
                        n = 30, levels = c(.5,.75,.95))

plot(myelk_locoh$locoh[3:1,"level"])

