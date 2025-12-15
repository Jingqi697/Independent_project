#### Libaries
library(data.table)
library(ggplot2)
library(maps)       
library(SeqArray)   
library(foreach)    
library(doParallel)

#### Data loading
meta_file <- "/scratch/cqh6wn/Class/full_sample_metadata.90Sept2025_ExpEvo.csv"
gds_file  <- "/scratch/cqh6wn/Class/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds"
### Filter metadata
samps <- fread(meta_file, fill=T)
valid_types <- c("wild", "Wild flies", "F1")
samps <- samps[Recommendation == "Pass" & fly_type %in% valid_types & set != "ExpEvo" & city != "Charlottesville"]
samps <- samps[order(jday)] # Order by date collected
candidates <- samps[, .SD[c(1, .N)], by = .(locality, year)] # group data by "locality, year", subset by selecting the first and last row
## Define spring and fall by the gap between jday
candidates[, gap := max(jday) - min(jday), by = .(locality, year)]
final_samps <- candidates[gap > 60]
final_samps[, season_defined := ifelse(jday == min(jday), "Spring", "Fall"), by = .(locality, year)]

#### Contrast
final_samps <- final_samps[order(sampleId)]
final_samps[, contrast_code := ifelse(season_defined == "Spring", 1, -1)] # If the season is Spring, assign value 1

#### Map coordinates
map_samps <- final_samps[!is.na(lat) & !is.na(long)] # Filter to non-NA latitude and longitude 
map_data <- map_samps[, .N, by = .(locality, lat, long)] # Aggregate by locality, lat, long, show the number of samples fall into each combination
map_data[, lat := as.numeric(lat)]
map_data[, long := as.numeric(long)]

world_map <- map_data("world")
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "gray60", linewidth = 0.3) +
  geom_point(data = map_data, aes(x = long, y = lat, size = N), 
             color = "darkred", alpha = 0.8) +
  scale_size_continuous(range = c(1, 3)) + # Adjusting point size
  labs(x = "Longitude", y = "Latitude", size = "Seasonal\nSamples") +
  theme_bw() 

#### Unique continents
unique(final_samps$continent)
