
# image <- image_df$file[3]
# rast(paste0(area_dir,"/",image_df$file[1]))
#
check_raster <- function(image, image_dir){
  
  xx <- try(terra::rast(paste0(image_dir,"/",image)))
  
  if(class(xx)[1] == "try-error"){
    return(image)
  } else {
    
    xx <- try(terra::values(terra::rast(paste0(image_dir,"/",image))[[1]]))
    
    if(class(xx)[1] == "try-error"){
      return(image)
    } else {
      return(NULL)
    }
  }
}


make_vsicurl_url <- function(base_url) {
  paste0(
    "/vsicurl", 
    "?pc_url_signing=yes",
    "&pc_collection=sentinel-2-l2a",
    "&url=",
    base_url
  )
}

extract_sentinel2_stac <- function(aoi,
                                 epsg,
                                 excl_dates,
                                 site_name,
                                 start_date, 
                                 end_date, 
                                 months = 1:12,
                                 minclouds,
                                 area_dir,
                                 workers){
  sf_use_s2(TRUE)
  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  
  it_obj <- s_obj %>% 
    stac_search(collections = "sentinel-2-l2a",
                bbox = st_bbox(aoi %>% st_transform(4326)),
                datetime = paste0(start_date,"/",end_date),
                limit = 1000) %>%
    get_request()
  
  if(length(it_obj$features) == 1000){
    dr <- tibble(date = seq.Date(as.Date(start_date), as.Date(end_date), by = "day")) %>% 
      mutate(gr = cut_number(row_number(.), 20)) %>% 
      group_by(gr) %>% 
      group_split()
    
    it_obj <- s_obj %>% 
      stac_search(collections = "sentinel-2-l2a",
                  bbox = st_bbox(aoi %>% st_transform(4326)),
                  datetime = paste0(min(dr[[1]]$date),"/",max(dr[[1]]$date)),
                  limit = 1000) %>%
      get_request()
    
    for(i in 2:length(dr)){
      # i <- 2
      it_obj2 <- s_obj %>% 
        stac_search(collections = "sentinel-2-l2a",
                    bbox = st_bbox(aoi %>% st_transform(4326)),
                    datetime = paste0(min(dr[[i]]$date),"/",max(dr[[i]]$date)),
                    limit = 1000) %>%
        get_request()
      
      it_obj$features <- c(it_obj$features, it_obj2$features)
      
    }
    
  }
  
  it_obj <- it_obj %>%
    items_filter(filter_fn = function(x) {x$properties$`eo:cloud_cover` < minclouds}) %>% 
    assets_select(asset_names = c("B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B11","B12","SCL")) %>%
    items_filter(filter_fn = function(x) {month(ymd_hms(x$properties$datetime)) %in% months}) %>%
    items_filter(filter_fn = function(x) {!as_date(ymd_hms(x$properties$datetime)) %in% excl_dates$date})
  
  if(length(it_obj$features) > 0){
    
    itst <- items_as_sf(it_obj) %>% st_make_valid()
    e <- try({
      itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
    })
    if(class(e)[[1]] == "try-error"){
      sf_use_s2(FALSE)
      itst <- st_crop(itst, aoi %>% st_transform(st_crs(itst)))
      sf_use_s2(TRUE)
    }
    
    
    itst$area <- as.numeric(st_area(itst))/(1000*1000)
    
    if(NA %in% names(itst)){
      itst <- itst[,!is.na(names(itst))]
    }
    
    tt <- itst %>% 
      mutate(date = as_date(ymd_hms(datetime))) %>% 
      group_by(date, instruments, `s2:datatake_id`, .add = TRUE) %>% 
      mutate(n = n()) %>% 
      arrange(desc(n)) %>% mutate(gid = cur_group_id()) %>% group_split()
    
    suppressMessages({
      tt <- lapply(tt, function(x){
        # x <- tt[[1]]
        if(nrow(x) > 1){
          if(nrow(x) < 3){
            if(diff(x %>% pull(area)) == 0){
              if(epsg %in% (x %>% pull(`proj:epsg`))){
                x <- x %>% filter(`proj:epsg` == epsg) %>% slice(1)
              } else {
                x <- x %>% slice(1)
              }
            } else {
              e <- try({
                stint <- st_intersection(x) %>% st_make_valid() %>% slice(2) %>% st_make_valid() %>% st_area %>% as.numeric()
              }, silent = TRUE)
              if(class(e)[[1]] == "try-error"){
                sf_use_s2(FALSE)
                stint <- st_intersection(x) %>% st_make_valid() %>% slice(2) %>% st_make_valid() %>% st_area %>% as.numeric()
                sf_use_s2(TRUE)
              }
              x <- x %>% filter(area > stint/(1000*1000)+0.1)
            }
          } else {
            x <- x %>% arrange(desc(area), (`s2:mean_solar_zenith`)) %>% slice(1)
          }
        }
        return(x)
      }) %>% bind_rows()
    })
    
    if(nrow(tt) > 0){
      itst <- itst %>% 
        filter(`s2:product_uri` %in% (tt %>% pull(`s2:product_uri`)))
      
      it_obj <- it_obj %>%
        items_filter(filter_fn = function(x) {x$properties$`s2:product_uri` %in% (itst %>% pull(`s2:product_uri`))})
      
      juuh <- mclapply(it_obj$features, function(ft){
        # ft <- it_obj$features[[1]]
        
        nm <- paste0(site_name, "_",
                     paste(str_split(ft$id, "_")[[1]][2], collapse = "_"), "_",
                     paste(str_split(ft$id, "_")[[1]][1], collapse = "_"), "_",
                     paste(str_split(ft$id, "_")[[1]][4:5], collapse = "_"), "_",
                     paste(substr(str_split(ft$id, "_")[[1]][3], 1, 8), collapse = "_"), "_",
                     gsub(":","",substr(ft$properties$datetime, 12, 19)), "GMT.tif")
        
        if(!file.exists(paste0(area_dir,"/", nm))){
          # print(ft$id)
          full_url <- make_vsicurl_url(assets_url(ft) %>% sort)
          file_names <- gsub("TIF$","tif",basename(full_url))
          
          juuh <- lapply(seq_len(length(full_url)), function(nr){
            e <- try({
              gdal_utils(
                "warp",
                source = full_url[[nr]],
                destination = paste0(tempdir(),"/",file_names[[nr]]),
                options = c(
                  "-t_srs", sf::st_crs(aoi)$wkt,
                  "-te", sf::st_bbox(aoi),
                  "-tr", c(10, 10)
                )
              )
            }, silent = TRUE)
            if(class(e)[[1]] == "try-error"){
              return(FALSE)
            } else {
              return(TRUE)
            }
          })
          
          err <- file_names[!unlist(juuh)]
          if(length(err) > 0){
            ll <- lapply(err, function(xx){
              #xx <- err[1]
              r <- rast(paste0(tempdir(),"/",file_names[unlist(juuh)][1]))
              r[] <- NA
              names(r) <- gsub(".tif","",xx)
              writeRaster(r, paste0(tempdir(),"/",xx), datatype = "INT2U", overwrite = TRUE)
            })
          }
          
          # if(length(file_names[grepl("_ST_", file_names)]) == 0){
          #   r <- rast(paste0(tempdir(),"/",file_names[grepl("_SR_B5", file_names)]))
          #   r[] <- NA
          #   names(r) <- gsub("_SR_B5.tif","_ST_B10",file_names[grepl("_SR_B5", file_names)])
          #   writeRaster(r, paste0(tempdir(),"/",gsub("_SR_B5.tif","_ST_B10.tif",file_names[grepl("_SR_B5", file_names)])), datatype = "INT2U", overwrite = TRUE)
          #   file_names <- c(file_names, gsub("_SR_B5.tif","_ST_B10.tif",file_names[grepl("_SR_B5", file_names)]))
          # }
          
          r <- rast(c(paste0(tempdir(),"/",sort(file_names[grepl("_B", file_names)])),
                      paste0(tempdir(),"/",file_names[grepl("_SCL", file_names)])))
          names(r) <- lapply(names(r), function(x) paste(str_split(x, "_")[[1]][3], collapse = "_")) %>% unlist
          
          writeRaster(r, paste0(area_dir,"/", nm), datatype = "INT2U", overwrite = TRUE)
          
          unlink(c(paste0(tempdir(),"/",file_names)))
        }
      }, mc.cores = workers)
      
      print(paste0(length(juuh %>% unlist()), " images downloaded!"))
    } else {
      print("No new Landsat scenes downloaded!")
    }
  } else {
    print("No new Landsat scenes downloaded!")
  }
}