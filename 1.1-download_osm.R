library(sf)
# library(mapview)
library(data.table)
library(dplyr)
library(purrr)

# open ghsl table
ghsl_table <- fread("input_data/pbf/pbf_hdc_country.csv")


# create the file to open
# extract date
# year <- format(Sys.Date(), "%Y")
year <- "2024"
# month <- format(Sys.Date(), "%m")
month <- "12"

# create the url to download countries from ghsl
ghsl_countries <- distinct(ghsl_table, region, country1) %>%
  mutate(file_dir = sprintf("input_data/pbf/%s/%s_%s%s.osm.pbf", year, country1, year, month)) %>%
  # create url
  mutate(url = sprintf("https://download.geofabrik.de/%s/%s-latest.osm.pbf", region, country1)) %>%
  # for russia only, it doesn't have the region
  mutate(url = ifelse(country1 == "russia", sprintf("https://download.geofabrik.de/%s-latest.osm.pbf", country1), url)) %>%
  # mutate(url = sprintf("https://osm-internal.download.geofabrik.de/%s/%s-internal.osh.pbf", region, country1)) %>%
  arrange(region)


# # check if the urls work
# library(purrr)
# response <- lapply(ghsl_countries$url,
#   RCurl::url.exists)
# names(response) <- ghsl_countries$country
# response_error <- Filter(isFALSE, response)



# to download - south america to test
# url_download <- ghsl_countries %>% filter(country1 == "russia")
url_download <- ghsl_countries %>% filter(region == "central-america")
response_download <- purrr::map2(url_download$url,
                                 url_download$file_dir,
                                 possibly(curl::curl_download, quiet = FALSE, otherwise = "erro")
)


# # check if everything was downloaded
# countries <- dir("input_data/pbf/2024")
# countries <- sub(pattern = "(\\w)_(.*$)", replacement = "\\1", x = countries) 
# left <- setdiff(ghsl_countries$country1, countries)

