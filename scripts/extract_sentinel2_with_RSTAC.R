library(sf)
library(tidyverse)
library(lubridate)
library(foreach)
library(doParallel)
library(R.utils)
library(terra)
library(rstac)

source("scripts/FUNCTIONS.R")

workers <- length(future::availableWorkers())
print(workers)

# Sequence of dates to query the imagery
starting_date <- '2023-01-01'
ending_date <- '2023-12-31'
minclouds <- 80
dl_months <- 5:9

polygon_buffer <- 10000 # if AOI's are points give the buffer radius in metres. If your AOI's are polygons, ignore.

# Direction where imagery will be stores (in area-specific sub-directories)
base_dir <- "/scratch/project_2003061/sentinel"

utmall <- st_read("aoi/utm_zones.gpkg") %>% 
  filter(ZONE != 0)

# Read in the study area points

aoi_all <- st_read("aoi/CA_plots_summer2023_1km_buffer.gpkg") %>%
  st_transform(crs = 4326) %>% 
  rename(name = Farm)

for(iAREA in aoi_all$name){
  # iAREA <- "CA-HV010"
  print(iAREA)
  print(Sys.time())
  
  aoi <- aoi_all %>% 
    filter(name == iAREA)
  
  # WGS84 UTM zones to set the correct projection
  utm <- utmall[aoi,] # Which zone the study points falls in
  
  if(nrow(utm) == 0){
    utm <- utmall[st_nearest_feature(aoi, utmall),]
  }
  
  lat <- st_coordinates(st_centroid(aoi))[,"Y"]
  utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0",utm$ZONE), utm$ZONE)
  epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
  
  # From point to polygon (20km x 20km)
  if(st_geometry_type(aoi) %in% c("POINT","MULTIPOINT")){
    aoi <- aoi %>% st_transform(crs = epsg) %>% 
      st_buffer(polygon_buffer) %>% st_bbox() %>% 
      st_as_sfc() %>% st_as_sf() %>% 
      mutate(name = aoi$name)
  } else {
    aoi <- aoi %>% st_transform(crs = epsg)
  }
  
  # Create sub-directory where imagery will be saved if does not exists
  area_dir <- paste0(base_dir,"/",aoi$name)
  if(!dir.exists(area_dir)){
    dir.create(area_dir)
  }
  
  # First round of downloads
  
  # Excisting files
  tifs <- list.files(area_dir, pattern = "GMT.tif$")
  tif_dates <- ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][6])))
  downloaded_dates <- tibble(date = unique(c(tif_dates)))
  
  e <- try({
    suppressMessages({
      extract_sentinel2_stac(aoi = aoi,
                           epsg = epsg,
                           excl_dates = downloaded_dates,
                           site_name = iAREA,
                           start_date = starting_date, 
                           end_date = ending_date, 
                           months = dl_months,
                           minclouds = minclouds,
                           area_dir = area_dir,
                           workers = workers)
    })
  })
  
  if(class(e)[[1]] == "try-error"){
    tifs <- list.files(area_dir, pattern = "GMT.tif$")
    tif_dates <- ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][6])))
    downloaded_dates <- tibble(date = unique(c(tif_dates)))
    e <- try({
      suppressMessages({
        extract_sentinel2_stac(aoi = aoi,
                               epsg = epsg,
                               excl_dates = downloaded_dates,
                               site_name = iAREA,
                               start_date = starting_date, 
                               end_date = ending_date, 
                               months = dl_months,
                               minclouds = minclouds,
                               area_dir = area_dir,
                               workers = workers)
      })
    })
  }
  
  # List downloaded images
  tifs <- list.files(area_dir, pattern = "GMT.tif$")
  
  if(length(tifs) > 0){
    
    tif_dates <- ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][6])))
    image_df <- tibble(file = tifs,
                       date = tif_dates)
    
    # Check if all rasters are working. Remove corrupted files.
    img_remove <- unlist(mclapply(image_df$file, check_raster, image_dir = area_dir, mc.cores = workers))
    unlink(paste0(area_dir,"/",img_remove))
    
    e <- try({
      suppressMessages({
        extract_sentinel2_stac(aoi = aoi,
                               epsg = epsg,
                               excl_dates = image_df,
                               site_name = iAREA,
                               start_date = starting_date, 
                               end_date = ending_date, 
                               months = dl_months,
                               minclouds = minclouds,
                               area_dir = area_dir,
                               workers = workers)
      })
    })
    
    tifs <- list.files(area_dir, pattern = "GMT.tif$")
    tif_dates <- ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][6])))
    image_df <- tibble(file = tifs,
                       date = tif_dates)
    
    img_remove <- unlist(mclapply(image_df$file, check_raster, image_dir = area_dir, mc.cores = workers))
    
    if(length(img_remove) > 0){
      print(paste0(length(img_remove), " raster(s) not functional. REMOVED!!"))
    }
    unlink(paste0(area_dir,"/",img_remove))
    
    image_df <- image_df %>% 
      filter(!file %in% img_remove)
    
  }
  unlink(tmpFiles())
}

# Calculate NDVIs, stack them by area and write to the base folder

for(iAREA in aoi_all$name){
  # iAREA <- "CA-HV068"
  print(iAREA)
  print(Sys.time())
  
  aoi <- aoi_all %>% 
    filter(name == iAREA)
  area_dir <- paste0(base_dir,"/",aoi$name)
  
  # From point to polygon (20km x 20km)
  if(st_geometry_type(aoi) %in% c("POINT","MULTIPOINT")){
    aoi <- aoi %>% st_transform(crs = epsg) %>% 
      st_buffer(polygon_buffer) %>% st_bbox() %>% 
      st_as_sfc() %>% st_as_sf() %>% 
      mutate(name = aoi$name)
  } else {
    aoi <- aoi %>% st_transform(crs = epsg)
  }
  
  # List downloaded images
  tifs <- list.files(area_dir, pattern = "GMT.tif$")
  
  if(length(tifs) > 0){
    
    tif_dates <- ymd(unlist(lapply(tifs, function(x) str_split(x, "_")[[1]][6])))
    image_df <- tibble(file = tifs,
                       date = tif_dates) %>% 
      arrange(date)
    
    # Calculate NDVIs
    rs <- lapply(image_df$file, function(image){
      # image <- image_df$file[1]
      r <- rast(paste0(area_dir, "/", image))
      ndvi <- (r$B08 - r$B04) / (r$B08 + r$B04)
      ndvi <- mask(ndvi, vect(aoi %>% st_transform(crs(r))))
      ndvi[!r$SCL %in% c(4,5)] <- NA
      ndvi[ndvi > 1] <- NA
      ndvi[ndvi < (-1)] <- NA
      names(ndvi) <- str_split(image, "_")[[1]][6]
      varnames(ndvi) <- "NDVI"
      
      if(sum(!is.na(values(ndvi, mat = FALSE))) >= 1000){
        return(ndvi)
      } else {
        return(NULL)
      }
    }) %>% rast
    
    # If same date has multiple images, replace those by their average
    if(any(duplicated(names(rs)))){
      ddates <- names(table(names(rs))[which(table(names(rs)) > 1)])
      for(dd in ddates){
        rt <- mean(rs[[names(rs) == dd]], na.rm = TRUE)
        names(rt) <- dd
        varnames(ndvi) <- "NDVI"
        rs <- c(rs[[names(rs) != dd]], rt)
      }
    }
    
    # Order by names (date)
    rs <- rs[[names(rs) %>% sort]]
    
    writeRaster(rs, paste0(base_dir,"/", iAREA, "_NDVI.tif"), overwrite = TRUE)
  }
  unlink(tmpFiles())
}

