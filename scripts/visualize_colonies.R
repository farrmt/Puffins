#---
#"Visualizing TUPU Colony Data"
#"Lisanne Petracca"
#"Nov 2022"
#---

#This is script summarizing the available datasets for where TUPU colonies occur in space (and for what temporal window)
#A note that many of these colonies no longer have TUPU and are historical

library(here)
library(ggplot2)
library(sf)
library(tidyverse)
library(sf)
library(ggmap)
library(MetBrewer)
library(RColorBrewer)
library(lubridate)

getwd()

#reading in data from WA and AK
data_AK <- read.csv("data/AK/AMNWR_tufted_puffin_colonies.csv", header=T)
data_WA1 <- read.csv("data/WA/S_Krock/QRY_export_for_UW_20220601.csv", header=T)
data_WA2 <- read.csv("data/WA/J_Brusa/TUPUPelagic.csv", header=T)
data_WA3 <- read.csv("data/WA/J_Brusa/TUPUSalishSea.csv", header=T) #no coordinates
data_WA4 <- read.csv("data/WA/J_Brusa/TUPUZone2No2013.csv", header=T)
data_WA5 <- read.csv("data/WA/S_Pearson/TUPU_Master_2016_15March2017.csv", header=T)
data_OR1 <- read.csv("data/OR/S_Stephensen/OR_SpeciesSurvey.csv", header=T)
data_OR2 <- read.csv("data/OR/S_Stephensen/OR_ColonyLocations.csv", header=T)
data_BC <- read.csv("data/BC/BC_SeabirdColonyDatabase.csv", header=T)
data_CA1 <- read.csv("data/CA/AllTUPUdata_McChesneyRevisions_2019-09-19.csv", header=T)
data_CA2 <- read.csv("data/CA/ColonyLocations_withCCN.csv", header=T)

#reading in regional datasets
NPSDP <- read.csv("data/N_Pacific_Seabird_Data_Portal/seabird_data_download/status_record_download.csv", header=T)
head(NPSDP)
unique(NPSDP$country)
NPSDP <- NPSDP %>% filter(country=="United States" | country=="Canada")

#getting coordinate column names
head(data_AK)
tail(data_AK)
head(data_WA1)
head(data_WA2)
head(data_WA3) #no coordinates in Salish Sea dataset
head(data_WA4) #no coordinates in TUPU Zone 2 dataset for "all years"; coords wo 2013 only
head(data_WA5)
head(data_OR1)
head(data_OR2)
head(NPSDP)
head(data_BC)
head(data_CA1)
head(data_CA2)

#renaming column names for the merge 
#also getting new column for data type
newColumnName <- "source"

#need to combine Oregon data into 1
data_OR1 <- data_OR1 %>% rename(survey_year = intSurveyYear) %>% 
  mutate(!!newColumnName := "All OR Colony Register")
data_OR2 <- data_OR2 %>% rename(long = dblLongitude) %>% rename(lat = dblLatitude) 
data_OR1 <- left_join(data_OR1, data_OR2, by="strColonyCode")
data_OR1 <- data_OR1 %>% rename(site_name = strColonyCode)

#need to combine CA data into 1
head(data_CA1)
data_CA1 <- data_CA1 %>% mutate(RevisedSurveyDate_format = as.Date(RevisedSurveyDate, format = "%d-%b-%y")) %>% 
  mutate(survey_year = year(RevisedSurveyDate_format)) 
data_CA2 <- data_CA2 %>% rename(long = Longitude) %>% rename(lat = Latitude) 
data_CA1 <- left_join(data_CA2, data_CA1, by="USFWSColonyCode")
data_CA1 <- data_CA1 %>% rename(site_name = USFWSColonyCode) %>% mutate(!!newColumnName := "All CA Colony Register")
tail(data_CA1)

#now we rename columns from the others
data_AK <- data_AK %>% rename(long = location_midpoint_lng) %>% rename(lat = location_midpoint_lat) %>%  
  rename(site_name = location_index) %>% rename(survey_year = survey_end_year) %>% 
  mutate(!!newColumnName := "TUPU AK Colony Register")
data_WA1 <- data_WA1 %>% rename(long = LongDecDeg) %>% rename(lat = LatDecDeg ) %>%  
  rename(site_name = SITE_NAME) %>% rename(survey_year = YEAR_) %>% 
  mutate(!!newColumnName := "All WA Colony Register")
data_WA2 <- data_WA2 %>% rename(long = Longitude) %>% rename(lat = Latitude ) %>%  
  rename(survey_year = YYYY) %>% 
  mutate(!!newColumnName := "All Pelagic At-Sea")
data_WA4 <- data_WA4 %>% rename(long = Longitude) %>% rename(lat = Latitude ) %>%  
  rename(survey_year = YYYY) %>% 
  mutate(!!newColumnName := "All Zone 2 Nearshore")
data_WA5 <- data_WA5 %>% rename(long = Long) %>% rename(lat = Lat ) %>%  
  rename(site_name = Island) %>% rename(survey_year = Year) %>% 
  mutate(!!newColumnName := "All WDFW At-Sea")
NPSDP <- NPSDP %>% rename(long = location_midpoint_lng) %>% rename(lat = location_midpoint_lat ) %>%  
  rename(site_name = location_id) %>% rename(survey_year = survey_end_year) %>% 
  mutate(!!newColumnName := "All N Pacific Seabird Data Portal")
data_BC <- data_BC %>% rename(long = Longitude) %>% rename(lat = Latitude) %>%  
  rename(site_name = SiteName1) %>% rename(survey_year = SurveyYear) %>% 
  mutate(!!newColumnName := "All BC Colony Register")

#get rid of rows without locations in two datasets
data_WA4 <- data_WA4 %>% filter(!is.na(long)) %>% filter(!is.na(lat))
data_WA5 <- data_WA5 %>% filter(!is.na(long)) %>% filter(!is.na(lat))

#let's make versions where there are TUPU records

head(data_WA4) #no coordinates in TUPU Zone 2 dataset for "all years"; coords wo 2013 only

data_WA1_TUPU <- data_WA1 %>% filter(SPECIES_CO=="TUPU") %>%  mutate(!!newColumnName := "TUPU WA Colony Register")
data_WA2_TUPU <- data_WA2 %>% filter(Species=="TUPU") %>%  mutate(!!newColumnName := "TUPU Pelagic At-Sea")
data_WA4_TUPU <- data_WA4 %>% filter(Species=="TUPU") %>%  mutate(!!newColumnName := "TUPU Zone 2 Nearshore")
data_WA5_TUPU <- data_WA1 %>% filter(SPECIES_CO=="TUPU") %>%  mutate(!!newColumnName := "TUPU WDFW At-Sea")
NPSDP_TUPU <- NPSDP %>% filter(species_common_name=="Tufted Puffin") %>%  mutate(!!newColumnName := "TUPU N Pacific Seabird Data Portal")
NPSDP_AK <- NPSDP %>% filter(state_province == "Alaska") %>% mutate(!!newColumnName := "All AK Colony Register")
data_OR1_TUPU <- data_OR1 %>% filter(strAOUCode =="TUPU") %>% mutate(!!newColumnName := "TUPU OR Colony Register")
data_BC_TUPU <- data_BC %>% filter(SpeciesID =="TUPU") %>% mutate(!!newColumnName := "TUPU BC Colony Register")
data_CA1_TUPU <- data_CA1 %>% filter(AOUCode =="TUPU") %>% mutate(!!newColumnName := "TUPU CA Colony Register")
head(data_CA1_TUPU)

#let's merge these datasets together
df_list <- list(NPSDP_AK, data_AK, 
                data_WA1, data_WA1_TUPU,
                data_WA2, data_WA2_TUPU,
                data_WA4, data_WA4_TUPU,
                data_WA5, data_WA5_TUPU,
                data_OR1, data_OR1_TUPU,
                data_BC, data_BC_TUPU,
                data_CA1, data_CA1_TUPU)      

#merge all data frames together
data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  
head(data)

data$lat <- as.numeric(data$lat)
data <- data %>% filter(!is.na(data$lat))

####---- LET'S DO A MAP WITH ALL DATA SOURCES ----####

#let's transfer to an sf object and assign a coordinate system
TUPU <- st_as_sf(data, coords = 
                         c("long", "lat"), crs = 4326)

#these need to be transformed to display properly w google basemap (weird)
TUPU_proj <- st_transform(TUPU, crs =3857)

#get bounding box in WGS84
area <- c(left = -180, bottom = 28, right = -118, top = 73)

#get stamen map for bounding box
basemap <- get_stamenmap(area, zoom = 3, maptype = "terrain") 

#here is a function that will allow sf objects to plot properly on google basemaps
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
  # Convert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

#using the function
myMap <- ggmap_bbox(basemap) 

head(TUPU_proj)
#plotting the colonies

#select only columns of interest for mapping
TUPU_proj <- TUPU_proj %>% dplyr::select(survey_year, source, site_name)

colors <- c('#e6550d','#1f78b4','#fc8d59','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#a1d99b','#6a3d9a','#ffff99','#b15928',
            '#000000', '#808080', '#fa9fb5', '#c51b8a')
dev.off()
ggmap(myMap) +
  geom_sf(data = TUPU_proj, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources in AK, BC, WA, OR, and CA")
ggsave("data/a_figures/AK_and_BC_and_WA_and_OR_and_CA.jpg")

####---- LETS DO A MAP WITH NPSDP ONLY ####

#NPSDP corresponds to the North Pacific Seabird Data Portal

#let's merge these datasets together
df_list2 <- list(NPSDP, NPSDP_TUPU)      

#merge all data frames together
data2 <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list2)  

data2$lat <- as.numeric(data2$lat)
data2 <- data2 %>% filter(!is.na(data2$lat))

#let's transfer to an sf object and assign a coordinate system
NPSDP_TUPU <- st_as_sf(data2, coords = 
                         c("long", "lat"), crs = 4326)

#these need to be transformed to display properly w google basemap (weird)
NPSDP_TUPU_proj <- st_transform(NPSDP_TUPU, crs =3857)

#lets map
colors <- c("darkmagenta", "green3")
ggmap(myMap) +
  geom_sf(data = NPSDP_TUPU_proj, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "",
                     labels=c("all data", "TUPU only")) + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("North Pacific Seabird Data Portal")#+
#theme(legend.position="none")
ggsave("data/figures/NPSDP.jpg")

####---- LETS DO A AK MAP WITH SOURCES ####

#get rid of AK for WA map
TUPU_proj_AK <- TUPU_proj %>% filter(source == "All AK Colony Register" |
                                       source == "TUPU AK Colony Register")

#get bounding box for AK only in WGS84
AK_area <- c(left = -180, bottom = 50, right = -125, top = 73)
plot(AK_area)
#get stamen map for bounding box
basemap <- get_stamenmap(AK_area, zoom = 6, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- c('#d95f02', '#fecc5c')
ggmap(myMap) +
  geom_sf(data = TUPU_proj_AK, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources for AK")
ggsave("data/a_figures/AK.jpg")

####---- LETS DO A WA MAP WITH SOURCES ####

#get rid of AK for WA map
TUPU_proj_noAK <- TUPU_proj %>% filter(source != "All AK Colony Register" &
                                        source != "All BC Colony Register" &
                                        source != "All OR Colony Register" &
                                         source != "TUPU AK Colony Register" &
                                         source != "TUPU BC Colony Register" &
                                         source != "TUPU OR Colony Register")
unique(TUPU_proj_noAK$source)

#get bounding box **for WA only** in WGS84
WA_area <- c(left = -126.5, bottom = 46, right = -122, top = 49.5)

#get stamen map for bounding box
basemap <- get_stamenmap(WA_area, zoom = 7, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- c('#1f78b4','#33a02c','#000000','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
ggmap(myMap) +
  geom_sf(data = TUPU_proj_noAK, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources for WA")
ggsave("data/a_figures/WA.jpg")

####---- LETS DO AN OR MAP WITH SOURCES ####

#get rid of AK for WA map
TUPU_proj_OR <- TUPU_proj %>% filter(source == "All OR Colony Register" |
                                     source == "TUPU OR Colony Register")

#get bounding box **for WA AND OR only** in WGS84
OR_area <- c(left = -126.5, bottom = 41.75, right = -122, top = 46.5)
plot(OR_area)
#get stamen map for bounding box
basemap <- get_stamenmap(OR_area, zoom = 7, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- c('#d95f02', '#fecc5c')
ggmap(myMap) +
  geom_sf(data = TUPU_proj_OR, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources for OR")
ggsave("data/a_figures/OR.jpg")

####---- LETS DO A BC MAP WITH SOURCES ####

#get rid of AK for WA map
TUPU_proj_BC <- TUPU_proj %>% filter(source == "All BC Colony Register" |
                                       source == "TUPU BC Colony Register")

#get bounding box **for WA AND OR only** in WGS84
BC_area <- c(left = -135, bottom = 48, right = -120, top = 55)
plot(BC_area)
#get stamen map for bounding box
basemap <- get_stamenmap(BC_area, zoom = 6, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- c('#d95f02', '#fecc5c')
ggmap(myMap) +
  geom_sf(data = TUPU_proj_BC, aes(color=as.factor(source)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources for BC")
ggsave("data/a_figures/BC.jpg")

####---- LETS DO A CA MAP WITH SOURCES ####

#get rid of AK for WA map
TUPU_proj_CA <- TUPU_proj %>% filter(source == "All CA Colony Register" |
                                       source == "TUPU CA Colony Register")

#get bounding box **for WA AND OR only** in WGS84
CA_area <- c(left = -125, bottom = 32, right = -113.5, top = 42.5)
plot(CA_area)
#get stamen map for bounding box
basemap <- get_stamenmap(CA_area, zoom = 6, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- c('#d95f02', '#fecc5c')
ggmap(myMap) +
  geom_sf(data = TUPU_proj_CA, aes(color=as.factor(source)), size=2, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Data source") + 
  coord_sf(crs = st_crs(3857))+
  ggtitle("Puffin Data Sources for CA")
ggsave("data/a_figures/CA.jpg")

####---- LETS GET FINAL YEAR THERE WAS PUFFIN DATA ####

#using "data" created around L. 113

head(data)

#this will select cols of interest
data_simple <- data %>% dplyr::select(long, lat, survey_year, site_name, source) 

#this will subset to those w TUPU sightings
data_simple <- data_simple[!grepl("All", data_simple$source),]
unique(data_simple$source)

#get rid of those without survey year
data_simple <- data_simple %>% filter(!is.na(survey_year))

#need to group by last year for the colony sites; at-sea can remain as is
data_atsea <- data_simple[grepl("At-Sea|Nearshore", data_simple$source),]
data_colony <-  data_simple[!grepl("At-Sea|Nearshore", data_simple$source),]
data_colony <- data_colony %>% group_by(site_name) %>% slice_max(survey_year)

data_simple2 <- rbind(data_colony, data_atsea)

data_simple2$year_cat <- cut(data_simple2$survey_year, 
                   breaks=c(-Inf, 1899, 1949, 1999, 2009, Inf), 
                   labels=c("pre-1900","pre-1950","1950-2000", "2001-2010", "2011-present"))

#let's transfer to an sf object and assign a coordinate system
TUPU_year <- st_as_sf(data_simple2, coords = 
                   c("long", "lat"), crs = 4326)

#these need to be transformed to display properly w google basemap (weird)
TUPU_proj <- st_transform(TUPU_year, crs =3857)

#get bounding box in WGS84
area <- c(left = -180, bottom = 28, right = -118, top = 73)

#get stamen map for bounding box
basemap <- get_stamenmap(area, zoom = 3, maptype = "terrain") 

#here is a function that will allow sf objects to plot properly on google basemaps
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
  # Convert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

#using the function
myMap <- ggmap_bbox(basemap) 

head(TUPU_proj)
unique(TUPU_proj$year_cat)
#plotting the colonies

colors <- met.brewer("Austria", 5)
ggmap(myMap) +
  geom_sf(data = TUPU_proj, aes(color=as.factor(year_cat)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Last year of TUPU record") + 
  coord_sf(crs = st_crs(3857))#+
  #ggtitle("Puffin Data Sources in AK, WA, and OR")
ggsave("data/a_figures/LastYr_AK_and_BC_and_WA_and_OR_and_CA.jpg")

#subset to WA, BC, and OR
#get bounding box **for WA AND OR only** in WGS84
new_area <- c(left = -134, bottom = 41.5, right = -122, top = 54.5)
plot(new_area)
#get stamen map for bounding box
basemap <- get_stamenmap(new_area, zoom = 6, maptype = "terrain") 

#using the function
myMap <- ggmap_bbox(basemap) 

colors <- met.brewer("Austria", 5)
ggmap(myMap) +
  geom_sf(data = TUPU_proj, aes(color=as.factor(year_cat)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Last year of TUPU record") + 
  coord_sf(crs = st_crs(3857))#+
  #ggtitle("Puffin Data Sources in AK, WA, and OR")
ggsave("data/a_figures/LastYr_WA_and_BC_and_OR.jpg")

####---- LET'S MAKE 5 x 5 GRID FOR TUPU COLONIES ONLY ----####

#using data_colony on L. 258

#let's transfer to an sf object and assign a coordinate system
TUPU_colony <- st_as_sf(data_colony, coords = 
                        c("long", "lat"), crs = 4326)

TUPU_colony$year_cat <- cut(TUPU_colony$survey_year, 
                             breaks=c(-Inf, 1899, 1949, 1999, 2009, Inf), 
                             labels=c("pre-1900","pre-1950","1950-2000", "2001-2010", "2011-present"))

#these need to be transformed to display properly w google basemap (weird)
TUPU_proj <- st_transform(TUPU_colony, crs =3857)

#get bounding box in WGS84
area <- c(left = -180, bottom = 41, right = -118, top = 73)

#get stamen map for bounding box
basemap <- get_stamenmap(area, zoom = 3, maptype = "terrain") 

#here is a function that will allow sf objects to plot properly on google basemaps
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
  # Convert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

#using the function
myMap <- ggmap_bbox(basemap) 

head(TUPU_proj)

#plotting the colonies

colors <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')
ggmap(myMap) +
  geom_sf(data = TUPU_proj, aes(color=as.factor(year_cat)), size=1.5, inherit.aes=F)+
  scale_color_manual(values = colors, name = "Last year of TUPU record") + 
  coord_sf(crs = st_crs(3857))#+
#ggtitle("Puffin Data Sources in AK, WA, and OR")
ggsave("data/a_figures/LastYr_ColoniesOnly_AK_and_WA_and_OR.jpg")

#theyre already in meters, which is awesome (3857)
#3857 is Pseudo-Mercator (Web Mercator; Spherical Mercator) [Google Maps, etc.]

#ok doesn't work though bc there are some islands that technically go over the 180 deg longitude line
#will alaska albers work?
TUPU_proj <- st_transform(TUPU_proj, crs =3338)

#let's make grid over those colonies
#what if we wanted to do hexagons instead?
TUPU_grid_1600km2 <- st_make_grid(
  TUPU_proj,
  cellsize = 40000, #this makes grids of 1600 km2
  crs = 3338,
  what = "polygons",
  square = TRUE)

extent <- st_bbox(TUPU_grid_1600km2)

library(rnaturalearth)
test <- ne_download(scale=50, type="states", category="cultural",returnclass="sf")
states <- test %>% filter(name_en=="Alaska" | name_en=="Oregon" |
                             name_en=="Washington" | name_en=="British Columbia")
# states <- st_transform(states, crs =3338)

ggplot() +
  geom_sf(data=states, fill=NA, color="black", size=1)+
  geom_sf(data = TUPU_proj, color="blue", size=1.5, inherit.aes=F)+
  geom_sf(data=TUPU_grid_1600km2, fill=NA, color = "darkgrey", size=0.25)+
  coord_sf(crs=3338, xlim=c(extent[[1]], extent[[3]]), ylim=c(extent[[2]], extent[[4]]))+
  ggtitle("1600 km2 Grid Over Puffin Colonies in N America")
ggsave("data/a_figures/ColoniesOnly_w1600km2Grid.jpg")