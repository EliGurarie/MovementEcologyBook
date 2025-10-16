whale <- read.csv("whale.csv")

require(sf)
require(lubridate)

whale_sf <- whale |> 
  mutate(time = ymd_hms(datetime)) |> 
  st_as_sf(coords = c("Longitude", "Latitude"), 
           crs = 4326) |> st_transform(32619) 

whale_sf <- cbind(whale_sf, st_coordinates(whale_sf)) |> 
  mutate(Z = X +1i*Y)

with(whale_sf, scan_track(x = X, y = Y, time = time))

whale_sweep <- with(whale_sf, 
                    sweepRACVM(Z = Z, 
                               T = time, 
                               windowsize = 360, 
                               time.unit = "mins", 
                               windowstep = 10, 
                               model = "RACVM", 
                               progress = TRUE))
plotWindowSweep(whale_sweep)

whale_cp <- findCandidateChangePoints(windowsweep = whale_sweep, 
                                      clusterwidth = 35) 
abline(v = whale_cp)
whale_table <- getCPtable(whale_cp, modelset = "all", iterate = FALSE, 
                          criterion = "AIC") 
whale_phases <- estimatePhases(whale_table)
whale_summary <- summarizePhases(whale_phases)

plotPhaseList(whale_phases, cex = 4)

whale_sf$phase <- 
  cut(whale_sf$time |> as.numeric(), 
    c(0, whale_cp |> as.numeric(), Inf), 
    labels = whale_summary$phase) 

whale_with_phase <- merge(whale_sf, whale_summary, by = "phase")



library(ggplot2)
library(dplyr)

# Prepare data with speed and combined label
whale_plot_data <- whale_with_phase |>
  mutate(
    # Calculate speed (handling NAs in mu.x and mu.y)
    speed = sqrt(ifelse(is.na(mu.x), 0, mu.x)^2 + 
                   ifelse(is.na(mu.x), 0, mu.x)^2 + 
                   eta^2),
    
    # Create combined label for legend
    phase_model = paste(phase, model, sep = " - ")
  ) |>
  arrange(time)  # Ensure chronological order

# Create the plot
ggplot(whale_plot_data) +
  geom_path(aes(x = X, y = Y, 
                color = phase_model,
                linewidth = speed,
                group = 1),  # Single trajectory
            lineend = "round",
            linejoin = "round") +
  scale_color_viridis_d(name = "Phase - Model") +
  scale_linewidth_continuous(name = "Speed (m/s)",
                             range = c(1.5, 3)) +  # Adjust width range as needed
  coord_sf() +
  theme_minimal()




