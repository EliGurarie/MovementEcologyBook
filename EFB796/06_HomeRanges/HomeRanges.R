## ----setup, include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache = TRUE)


## -------------------------------------------------------------------------------
library(ggplot2)
library(plyr)
library(sf)
library(mapview)


## ----loadData, echo = -1--------------------------------------------------------
load("../data/elk_processed.rda")
head(elk_gps)


## ----SelectResidentElk----------------------------------------------------------
elk_res <- elk_gps |>
  subset(id %in% c("YL80", "YL91", "YL94")) |>
  mutate(id = droplevels(id))

elk_res_sf <- elk_res |> 
  st_as_sf(coords = c("lon","lat"), crs = 4326) |> 
  st_transform(32611) 

elk_res <- cbind(elk_res, st_coordinates(elk_res_sf))
str(elk_res)


## ----ThreeElkMapped-------------------------------------------------------------
ggplot(data=elk_res, aes(x=X, y=Y, col = id)) +
  geom_path() + coord_fixed() + ggtitle("Some (mostly) resident elk")


## ----cache = FALSE--------------------------------------------------------------
require(dplyr)
elk_tracks <- elk_res_sf |> 
  dplyr::group_by(id) |> 
  summarize(do_union=FALSE) |> 
  st_cast("LINESTRING")

mapview(elk_tracks,zcol="id")


## -------------------------------------------------------------------------------
xy <- (elk_res |> subset(id == id[1]))[,c("X","Y")]
head(xy)


## -------------------------------------------------------------------------------
chull(xy)


## -------------------------------------------------------------------------------
mcp <- xy[chull(xy),]


## ----echo = -1------------------------------------------------------------------
par(bty = "l", mgp = c(1,.1,0), tck = 0.01, cex.axis = 0.75)
plot(xy, cex = 0.5, asp = 1)
polygon(mcp, bor = 2, lwd = 2)


## -------------------------------------------------------------------------------
library(pracma)
polyarea(mcp[,1], mcp[,2])


## -------------------------------------------------------------------------------
-polyarea(mcp[,1], mcp[,2])/(1e6)


## -------------------------------------------------------------------------------
library(adehabitatHR)


## -------------------------------------------------------------------------------
myelk_sf <- elk_res_sf |> subset(id == "YL94")
myelk_sp <- myelk_sf |> as_Spatial()


## ----mcp95----------------------------------------------------------------------
(mcp95 <- mcp(myelk_sp, percent = 95, unout = "km2"))


## ----mcp100---------------------------------------------------------------------
(mcp100 <- mcp(myelk_sp, percent = 100, unout = "km2"))


## -------------------------------------------------------------------------------
mcp95 <- st_as_sf(mcp95)
mcp100 <- st_as_sf(mcp100)

ggplot() + 
  geom_sf(data = mcp100, color = "blue", fill = alpha("blue",.2)) +
    geom_sf(data = mcp95, color = "red", fill = alpha("red",.2)) + 
     geom_sf(data = myelk_sf, alpha = 0.2)


## ----findMCPFunction------------------------------------------------------------
findMCP <- function(id_sf, percent){
  id_sp <- id_sf |> as_Spatial()
  id_mcp <- mcp(id_sp, percent, unout="km2")
  return(st_as_sf(id_mcp))
}
findMCP(elk_res_sf |> subset(id == "YL80"), 95)
findMCP(elk_res_sf |> subset(id == "YL91"), 95)
findMCP(elk_res_sf |> subset(id == "YL94"), 95)


## ----cache = FALSE--------------------------------------------------------------
elk_sp <- elk_res_sf |> mutate(datetime = NULL) |> as_Spatial(IDs = "id") 
elk_mcps <- mcp(elk_sp, percent = 95, unout = "km2") |> st_as_sf()
mapview(elk_mcps, zcol = "id")


## ----PickAnElk------------------------------------------------------------------
myelk_sf <- elk_res_sf |> subset(id == "YL80")
myelk_sp <- myelk_sf |>  as_Spatial() 
myelk_kernelud <- kernelUD(myelk_sp, grid = 200)


## -------------------------------------------------------------------------------
plot(myelk_kernelud)


## -------------------------------------------------------------------------------
require(raster)
my_kernel <- raster(myelk_kernelud)
my_kernel
plot(my_kernel)


## ----mapviewKernel, cache = FALSE-----------------------------------------------
mapview(my_kernel)


## ----contourPlotKernel----------------------------------------------------------
plot(my_kernel)
contour(my_kernel, add = TRUE)


## ----PerspPlotKernel------------------------------------------------------------
persp(my_kernel, border = NA, shade = TRUE)


## ----getKDEpolygon--------------------------------------------------------------
myelk_kde_poly <- getverticeshr(myelk_kernelud, percent = 95) |>
  st_as_sf()
myelk_kde_poly


## ----getArea--------------------------------------------------------------------
st_area(myelk_kde_poly)


## ----mapview2, cache =FALSE-----------------------------------------------------
mapview(myelk_kde_poly) + mapview(myelk_sf)


## ----getKernelPolyFunction------------------------------------------------------
getKernelPoly <- function(sf, percent = 95, idcol = "id", ...){
  sp <- sf |> mutate(id = droplevels(get(idcol))) 
  as_Spatial(sp[,"id"], cast = TRUE, IDs = "id") |> kernelUD(...) |>
  getverticeshr(percent = 95) |> 
  st_as_sf()
}


## ----ComparingKernels-----------------------------------------------------------
kde_poly_norm <- getKernelPoly(elk_sf |> subset(id == "YL91"), kern = "bivnorm") 
kde_poly_epa <- getKernelPoly(elk_sf |> subset(id == "YL91"), kern = "epa") 

kde_compare_kernels <- rbind(
  kde_poly_norm |> mutate(type = "Bivariate Normal"),
  kde_poly_epa |> mutate(type = "Epanechnikov"))

ggplot(kde_compare_kernels) + geom_sf(aes(fill = type), alpha = .5) + 
    geom_sf(data = elk_sf |> subset(id == "YL91"), alpha = .2, size = 1)


## ----computeAllMCPs-------------------------------------------------------------
MCP_allElks <- getKernelPoly(elk_sf |> st_transform(32611), kern = "epa") 


## ----mapAllMCPs-----------------------------------------------------------------
ggplot(MCP_allElks) + 
  geom_sf(aes(fill = id, color = id), alpha = .2)


## ----AndSummarize---------------------------------------------------------------
MCP_allElks
summary(MCP_allElks$area)


## ----computeTwoKernels----------------------------------------------------------
elk1 <- elk_res_sf |> subset(id == "YL80")
elk2 <- elk_res_sf |> subset(id == "YL94")

kernel1 <- elk1 |>  as_Spatial() |>  kernelUD(grid = 200)
kernel2 <- elk2 |>  as_Spatial() |>  kernelUD(grid = 200)


## ----Countour2Kernels-----------------------------------------------------------
contour(kernel1, 
        xlim = c(590e3, 610e3), ylim = c(5720e3, 5745e3), 
        col = "blue")
 contour(kernel2, add = TRUE, col = "red")
axis(1); axis(2)


## -------------------------------------------------------------------------------
threeElks_sp <- as_Spatial(elk_res_sf[,"id"], cast = TRUE, IDs = "id")
kerneloverlap(threeElks_sp, method = "VI", kern = "epa")


## -------------------------------------------------------------------------------
require(amt)


## ----PrepareTLoCoHData----------------------------------------------------------
myelk <- elk_res |> subset(id == "YL80") |> 
  make_track(.x = X, .y = Y, .t = datetime, crs = st_crs(elk_res_sf))


## ----computeLoCoH---------------------------------------------------------------
myelk_locoh <- hr_locoh(myelk, n = 30, type = "k", levels = c(0.5,.75,.95))


## -------------------------------------------------------------------------------
plot(myelk_locoh$locoh[3:1,"level"])


## ----cache = FALSE--------------------------------------------------------------
myelk_sf <- elk_res_sf |> subset(id == "YL80")
mapview(myelk_locoh$locoh, zcol = "level") + mapview(myelk_sf, cex = 0.3)


## -------------------------------------------------------------------------------
st_area(myelk_locoh$locoh)


## ----eval = FALSE, echo = FALSE-------------------------------------------------
# knitr::purl("HomeRanges_MCP_KDE.Rmd")

