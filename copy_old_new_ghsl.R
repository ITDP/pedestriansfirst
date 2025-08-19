library(sf)
library(data.table)
library(dplyr)
library(mapview)
sf::sf_use_s2(FALSE)



# -------------------------------------------------------------------------


# filter only the cities we had preivously?
# extract the folder  from the old hdcs
folders <- dir("../atlas/data-raw/data_final/cities", full.names = TRUE)
# extract hdc
hdcs <- stringr::str_extract(folders, "ghsl_region_\\d{1,}")
hdcs <- stringr::str_remove(hdcs, "ghsl_region_")
# hdcs <- stringr::str_pad(hdcs, width = 5, side = "left", pad = 0)
  



# open the ghsl boundaries file
ghsl_old <- st_read("input_data/ghsl/SMOD_V1s6_opr_P2023_v1_2020_labelUC_DB_release.gpkg") %>%
  st_transform(crs = 4326) %>%
  # filter(POP_2020 > 50000) %>%
  select(hdc = ID_UC_G0, NAME_MAIN) %>%
  # filter(hdc %in% hdcs) %>%
  st_make_valid()

ghsl_new <- st_read("input_data/ghsl/ghsl_2024.gpkg") %>%
  st_transform(crs = 4326) %>%
  # filter(POP_2020 > 50000) %>%
  select(hdc = ID_UC_G0, NAME_MAIN) %>%
  st_make_valid()

# join
ghsl_join <- st_join(ghsl_new, ghsl_old, suffix = c("_new", "_old"), largest = TRUE)
  # st_set_geometry(NULL)

# test 1: new cities that don't intersect with the old cities
ghsl_join_prob1 <- filter(ghsl_join, is.na(hdc_old))
mapview(ghsl_join_prob1)

# test 2: cities that don't match the name
ghsl_join_prob2 <- filter(ghsl_join, NAME_MAIN_new != NAME_MAIN_old) %>%
  # filter some that I know are alright
  filter(!(NAME_MAIN_new %in% c("La Paz", "Panama City", "Taipei",
                              "Valparaíso", "San Francisco",
                              "Manila", "Osaka", "Brasilia",
                              "Islamabad", "Niterói", "New Delhi",
                              "New York City", "Vitória", "Guangzhou", "Guwahati")))


ghsl_join_prob2_old <- ghsl_old %>%
  filter(hdc %in% ghsl_join_prob2$hdc_old)
mapview(ghsl_join_prob2) + ghsl_join_prob2_old

filter(ghsl_join, is.na(hdc_old))
filter(ghsl_join, NAME_MAIN_new != NAME_MAIN_old)


ghsl_ok <- ghsl_join %>% 
  mutate(run_again = ifelse(hdc_new %in% ghsl_join_prob1$hdc_new | NAME_MAIN_new %in% ghsl_join_prob2$NAME_MAIN_new, TRUE, FALSE)) %>%
  st_set_geometry(NULL)

# save
data.table::fwrite(ghsl_ok, "input_data/hdc_to_run.csv")

# test 3: cities in x/y that match multiple cities
# multiple1 <- ghsl_ok$hdc_old[duplicated(ghsl_ok$hdc_old)]


# now, copy the city files ------------
ghsl_ok <- fread("input_data/hdc_to_run.csv",
                 colClasses = c("character", "character", "integer", "character", "logical")
                 )

# x <- folders[1]
create_files <- function(x) {
  
  a <- dir(x, recursive = TRUE, full.names = TRUE)
  # remove files
  library(stringr)
  a <- a[!(str_detect(a, "debug"))]
  a <- a[!(str_detect(a, "geodata/population"))]
  a <- a[!(str_detect(a, "indicator_values"))]
  a <- a[!(str_detect(a, "temp/boundaries.geojson"))]
  
  a <- data.frame(filename_old = a)
  
  
  
}

# create df
folders_files <- purrr::map_dfr(folders, create_files)

folders_files <- folders_files %>% 
  mutate(hdc_old = stringr::str_extract(filename_old, "ghsl_region_\\d{1,}")) %>%
  mutate(hdc_old = stringr::str_remove(hdc_old, "ghsl_region_"))
# bring the new hdc
folders_files <- folders_files %>%
  mutate(hdc_old = as.integer(hdc_old)) %>%
  left_join(ghsl_ok, by = c("hdc_old" = "hdc_old")) %>%
  filter(!run_again)


# extract file name
folders_files <- folders_files %>%
  mutate(filename_new = stringr::str_replace(filename_old, "(.*)/(ghsl_region_)(\\d{1,})/(.*)", 
                                             sprintf("\\2%s/\\4", hdc_new))) %>%
  # create new dir
  mutate(filename_new = sprintf("cities_out/%s", filename_new))
  

# copy and paste data
purrr::walk(unique(dirname(folders_files$filename_new)), dir.create, recursive = TRUE)
purrr::walk2(folders_files$filename_old, folders_files$filename_new, file.copy, overwrite = TRUE)
# remove the indicators files / pop

# # remove log files
# purrr::walk(dir("input_data/gtfs/2024", full.names = TRUE, pattern = ".txt"), file.remove)
# 
# 
# 
# # now, copy the gtfs files ----------------------
# # extract the gtfs files from the old hdcs
# gtfs_files <- dir("../atlas/data-raw/data_final/cities", full.names = TRUE)
# # extract hdc
# hdcs <- stringr::str_extract(gtfs_files, "ghsl_region_\\d{1,}")
# hdcs <- stringr::str_remove(hdcs, "ghsl_region_")
# gtfs_files <- paste0(gtfs_files, "/temp/gtfs")
# gtfs_files <- purrr::map(gtfs_files, dir, full.names = TRUE)
# names(gtfs_files) <- hdcs
# gtfs_files <- gtfs_files[lengths(gtfs_files)!=0]
# # create df
# gtfs_files <- lapply(gtfs_files, data.frame) %>% rbindlist(idcol = "hdc")
# gtfs_files <- gtfs_files %>% as.data.frame() %>% 
#   rename(filename = 2)
# # bring the new hdc
# gtfs_files <- gtfs_files %>%
#   mutate(hdc = as.integer(hdc)) %>%
#   left_join(ghsl_join, by = c("hdc" = "hdc_old"))
# 
# filter(gtfs_files, is.na(hdc_new))
# 
# # extract file name
# gtfs_files <- gtfs_files %>%
#   mutate(filename_new = stringr::str_replace(filename, "gtfs/(.*)", 
#                                              sprintf("gtfs/gtfs_%s_202404_\\1", hdc_new))) %>%
#   # create new dir
#   mutate(filename_new = sprintf("../pedestriansfirst/input_data/gtfs/2024/%s", basename(filename_new)))
# 
# # copy and paste data
# purrr::walk2(gtfs_files$filename, gtfs_files$filename_new, file.copy)
# # remove log files
# purrr::walk(dir("input_data/gtfs/2024", full.names = TRUE, pattern = ".txt"), file.remove)
# 
# # # open gtfs boundaries
# # gtfs_boundaries <- fread("input_data/gtfs/gtfs_sources.csv") %>%
# #   filter(!is.na(gtfs_boundaries$location.bounding_box.maximum_latitude))
# # bbox <- purrr::pmap(list("xmin" = gtfs_boundaries$location.bounding_box.minimum_latitude,
# #           "xmax" = gtfs_boundaries$location.bounding_box.maximum_latitude,
# #         "ymin" = gtfs_boundaries$location.bounding_box.minimum_longitude,
# #         "ymax" = gtfs_boundaries$location.bounding_box.maximum_longitude),
# #         c)
# # # create sf from boundaries
# # gtfs_boundaries <- lapply(bbox, function(x) st_as_sfc(st_bbox(x), crs=4326)) %>% rbindlist()
# library(sf)
# library(mapview)
# library(dplyr)
# 
# oi <- st_read("cities_out/ghsl_region_08154/geodata/allbike/allbike_latlon_2024.geojson")
# oi <- st_read("cities_out/ghsl_region_08154/geodata/protectedbike/protectedbike_latlon_2024.geojson")
# 
# mapview(oi)
# 
# a <- st_read("cities_out/ghsl_region_08154/indicator_values.gpkg") %>%
#   # filter(admin_level == 10) %>%
#   filter(is.na(admin_level))
# 
# 
# 
# mapview(a, zcol = "pnpb_2024") + oi
# oi1 <- st_read("cities_out/ghsl_region_08154/geodata/pnft/pnft_latlon_2024.geojson")
# mapview(oi1)
# 
# mapview(a, zcol = "h.s_2024")
