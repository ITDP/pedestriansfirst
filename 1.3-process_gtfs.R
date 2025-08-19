# transform frequencies gtfs files into stop_times

library(gtfsio) # remotes::install_github('kauebraga/gtfsio')
library(gtfstools)
library(dplyr)
library(sf)


# hdc = "05472"; year = 2024 # jakarta
# hdc = "07910"; year = 2024 # dc
# hdc = "01156"; year = 2023 # trujillo

trim_gtfs <- function(hdc, year) {
  
  
  files <- dir(sprintf("input_data/gtfs/%s", year), pattern = sprintf("gtfs_%s", hdc), full.names = TRUE)
  # open gtfs
  if (length(files) == 1) {
    
    # open gtfs
    gtfs <- read_gtfs(files)
    
  } else {
    
    # open gtfs
    gtfs <- lapply(files, read_gtfs)
    gtfs <- purrr::reduce(gtfs, gtfstools::merge_gtfs)
    
    
  }
  
  # trim to boundaries
  boundaries <- st_read(sprintf("/media/kauebraga/data/pedestriansfirst/cities_out/ghsl_region_%s/debug/area_for_osm_extract.gpkg", hdc))
  # mapview::mapview(boundaries)
  
  gtfs_new <- gtfstools::filter_by_spatial_extent(gtfs, boundaries)
  
  # save
  gtfstools::write_gtfs(gtfs_new, sprintf("data/%s/gtfs_%s_%s.zip", year, hdc, year))
  
}

# run only for the hdc with gtfs
gtfs_run <- dir(sprintf("input_data/gtfs/%s", year))
gtfs_run <- stringr::str_extract(gtfs_run, pattern = "\\d{5}")

library(furrr)
plan(multicore)
furrr::future_walk(gtfs_run, trim_gtfs)

  
  # file <- files[1]
  
  process_gtfs_hdc <- function(file) {
    
    # open gtfs
    gtfs <- read_gtfs(file)
    # convert to st
    gtfs1 <- frequencies_to_stop_times(gtfs)
    # save
    save <- stringr::str_remove(file, ".zip")
    save <- paste0(save, "_fixed.zip")
    write_gtfs(gtfs1, save)
    
    # move the old GTFS to the raw folder
    file.copy(from = file, to = paste0(dirname(file), "/raw/", basename(file)))
    file.remove(file)
    
  }
  
  purrr::walk(files, process_gtfs_hdc)
  
}





process_gtfs <- function(hdc, year) {
  
  
  files <- dir(sprintf("input_data/gtfs/%s", year), pattern = sprintf("gtfs_%s", hdc), full.names = TRUE)
  # check if we the gtfs are frequency based
  files_files <- purrr::map_lgl(files, function(x) isTRUE(nrow(zip::zip_list(x) %>% filter(filename == "frequencies.txt")) == 1))
  # then filter
  files <- files[files_files]
  
  # file <- files[1]
  
  process_gtfs_hdc <- function(file) {
    
    # open gtfs
    gtfs <- read_gtfs(file)
    # convert to st
    gtfs1 <- frequencies_to_stop_times(gtfs)
    # save
    save <- stringr::str_remove(file, ".zip")
    save <- paste0(save, "_fixed.zip")
    write_gtfs(gtfs1, save)
    
    # move the old GTFS to the raw folder
    file.copy(from = file, to = paste0(dirname(file), "/raw/", basename(file)))
    file.remove(file)
    
  }
  
  purrr::walk(files, process_gtfs_hdc)
  
}