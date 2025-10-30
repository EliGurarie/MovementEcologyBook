
a <- load("data/Cranberry_Data.rda")
# b <- load("data/ProcessedNEC.rda")

require(terra)
require(mapview)
require(raster)

  mapview(CranNEC_sf)
  plot(Cranberry_covariates)
  Cranberry_covariates <- dropLayer(Cranberry_covariates, 3)
  Cranberry_covariates[Cranberry_covariates == 0] <- NA
  

# transform movement data to raster covariates

  CranNEC_sf <- CranNEC_sf |> st_transform(crs(Cranberry_covariates))
  plot(Cranberry_covariates[[1]])
  plot(CranNEC_sf, add = TRUE, col = "darkred", pch = 19, cex = 0.5)
  
  
# create used set
  
  covars <- extract(Cranberry_covariates, CranNEC_sf)
  used_df <- data.frame(CranNEC_sf |> st_drop_geometry(), 
                        st_coordinates(CranNEC_sf),
                        covars, 
                        Used = TRUE)
  
# create null set

## from MCP
  
  xy <- st_coordinates(CranNEC_sf)
  mcp <- xy[chull(xy),] 
  mcp <- rbind(mcp, mcp[1,])
  mcp_sf <- st_polygon(list(mcp))

  plot(Cranberry_covariates[[1]])
  plot(CranNEC_sf, add = TRUE, col = "darkred", pch = 19, cex = 0.5)
  plot(mcp_sf, add = TRUE)
  
# can also buffer MCP if desired
  
  mcp_buffered <- st_buffer(mcp_sf, 20)
  plot(mcp_buffered, add = TRUE)
  
  
# sample points from within buffer
  
  # How many points to choose?   
  # The more - kind of - the better.  But not necessary
  # Let's do a 5 to 1 ratio
  
  n_used <- nrow(used_df)
  n_available <- n_used * 5
  
  available_pts <- st_sample(mcp_buffered, n_available)
  plot(null_pts, add = TRUE, col = "grey")
  
# obtain available covariates
  
  covars_available <- extract(Cranberry_covariates, st_coordinates(null_pts))
  available_df <- data.frame(Used = FALSE, 
                             st_coordinates(null_pts), 
                             covars_available)
  
# combine into a single data frame

  require(gtools)
  NEC_RSF_df <- smartbind(used_df, available_df)
  
  # remove all the NA's
  NEC_RSF_df <- NEC_RSF_df |> subset(!is.na(Canopy_cover))
  
  dim(NEC_RSF_df)
  head(NEC_RSF_df)


# pairs plot of covariates
  
  NEC_RSF_df[,c("Canopy_cover", "invasive_stem")] |> 
    pairs()
  
  
# compare empirical distributions of used and avialable
  
  with(NEC_RSF_df, {
          plot(density(Canopy_cover[Used]), col = "darkred", 
               main = "Canopy Cover")
          lines(density(Canopy_cover[!Used]), col = "darkgrey")
       })
      
  
  with(NEC_RSF_df, {
          plot(density(invasive_stem[Used]), col = "darkred", 
               main = "Invasive Stem")
          lines(density(invasive_stem[!Used]), col = "darkgrey")
       })
  

# Fit RSFs
  
  M1 <- glm(Used ~ Canopy_cover + invasive_stem, data = NEC_RSF_df, 
      family = "binomial")
  
  summary(M1)

# working with lists
  
  formulae <- list(
    M0 = Used~1,
    Canopy = Used~Canopy_cover,
    Invasive = Used~invasive_stem,
    CanInv = Used~Canopy_cover + invasive_stem,
    CanByInv = Used~Canopy_cover * invasive_stem
  )

  
  fits <- lapply(formulae, function(f)
                glm(f, data = NEC_RSF_df, family = "binomial"))  
  
  require(plyr)
  ldply(fits, AIC)

  
# make more complex model
  moreformulae <- list(
    C2 = Used ~ poly(Canopy_cover, 2),
    I2 = Used ~ poly(invasive_stem, 2),
    C2I2 = Used ~ poly(Canopy_cover, 2) + poly(invasive_stem, 2)
  )
  
  morefits <- lapply(moreformulae, function(f)
    glm(f, data = NEC_RSF_df, family = "binomial"))  
  
  ldply(morefits, AIC)
  
# combinations (even with two variables) get complex. 
# let's pull out a "big gun"  ... the dredge
  
  require(MuMIn)
  global.model <-  glm(Used ~ 
                         (Canopy_cover + I(Canopy_cover^2)) *
                         (invasive_stem + I(invasive_stem^2)), 
                       data = NEC_RSF_df, family = "binomial", 
                       na.action = "na.fail")

  NEC_dredge <- dredge(global.model)
  plot(NEC_dredge)
  
# good model
  
  fit <- glm(Used ~ (Canopy_cover + I(Canopy_cover^2)) +
               (invasive_stem + I(invasive_stem^2)) + 
               Canopy_cover * I(invasive_stem^2), 
             family = "binomial",
             data = NEC_RSF_df)
  
    
  summary(fit)

# Make prediction map
  
  covars_predict <- values(Cranberry_covariates)
  RSF_predict <- predict(fit, newdata = data.frame(covars_predict))
  
  RSF_raster <- Cranberry_covariates[[1]]*NA
  values(RSF_raster) <- RSF_predict

  plot(RSF_raster)
  
# break into 10 quantiles - worst to best
  
  breaks <- quantile(values(RSF_raster), seq(0,1, length = 11), na.rm = TRUE)
  
  plot(RSF_raster, breaks = breaks, col = gplots::rich.colors(10))
  plot(CranNEC_sf, add = TRUE)  
  
  RSF_raster_discrete <- cut(RSF_raster, breaks)
  plot(RSF_raster_discrete, col = gplots::rich.colors(10))
  plot(CranNEC_sf, add = TRUE, col = "white", pch = 4)  
  
  
  
# A plot of predictions showing interactions
  
  # 1. Effect of Canopy_cover (holding invasive_stem constant)
  canopy_seq <- seq(0, 60, length.out = 200)
  invasive_levels <- seq(0,180,30)
  colors <- gplots::rich.colors(length(invasive_levels))
  
  plot(canopy_seq, canopy_seq, type = "n", ylim = c(0, 1),
       xlab = "Canopy cover", ylab = "Predicted response",
       main = "Effect of Canopy cover")
  
  for(i in seq_along(invasive_levels)) {
    newdata <- data.frame(Canopy_cover = canopy_seq, 
                          invasive_stem = invasive_levels[i])
    pred <- predict(fit, newdata = newdata, type = "response")
    lines(canopy_seq, pred, col = colors[i], lwd = 2)
  }
  legend("right", legend = invasive_levels, title = "invasive stem",
         col = colors, lwd = 2, bty = "n")
  