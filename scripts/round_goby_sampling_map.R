### August 4, 2019
### Make a map of sampling sites in the Great Lakes/ Finger lakes region

setwd("/Users/kbja10/Github/eDNA_goby_popgen/")

# install and load all of the packages we will be using
# install.packages(c("raster", "maptools", "rgdal", "maps", "usmap"))
library(raster) 
library(maptools)
library(rgdal)
library(maps)
library(usmap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(inlmisc)

# Go to the website to download lakes data (shapefile):
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/
# scroll down to Lakes + Reservoirs
# click "Download lakes" and save zip file to your Downloads
# click on the zip file to unzip it (should be a file named ne_10m_lakes)

# Make lake objects 
lakes <- readOGR("/Users/kbja10/Downloads/ne_10m_lakes") # file path to lake shapefile
oneida <- subset(lakes, name=="Oneida Lake") # make object for Oneida Lake
finger_lakes <- subset(lakes, name_alt=="Finger Lakes") # make object for Finger Lakes
great_lakes <- subset(lakes, name_alt=="Great Lakes") # make object for Great Lakes

# Make canal objects
url.river_data <- url("http://sharpsightlabs.com/wp-content/datasets/usa_rivers.RData")
load(url.river_data)
lines.streams <- subset(lines.rivers, (FEATURE %in% c("Right Bank", "Left Bank", "Stream")))
lines.canals <- subset(lines.rivers, (FEATURE %in% "Canal"))
lines.streams <- subset(lines.streams, (STATE %in% c("WI","IL","MI","IN","OH","PA","NY")))
lines.canals <- subset(lines.canals, (STATE %in% c("WI","IL","MI","IN","OH","PA","NY")))
df.usa_streams <- fortify(lines.streams)
df.usa_canals <- fortify(lines.canals)
table(lines.rivers$FEATURE)

# Make a data frame of sampling points longitude/latitude
sites <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/Field_sites_and_data/Round_goby_field_collections_data_summer2019.csv", header=TRUE)
sites <- sites[sites$Sample_type=="eDNA",]
sampling_points <- data.frame(Site_ID=unique(sites$Site_ID),
                              lat=unique(sites$Latitude),
                              long=unique(sites$Longitude))
sampling_points <- sampling_points[-grep("BUF", sampling_points$Site_ID),]
sampling_points <- sampling_points[-grep("SEN", sampling_points$Site_ID),]
sampling_points$Site_ID <- factor(sampling_points$Site_ID, 
                                  levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE",
                                           "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
# Make state object 
all_states <- map_data("state")  

# PLOT
supp.labs <- c("Cayuga Lake","Cross Lake","Erie Canal","Lake Erie (east)","Lake Erie (west)",
               "Lake Huron","Lake Michigan (west)","Lake Michigan (east)","Lake St. Clair",
               "Oneida Lake","Onondaga Lake","Lake Ontario (east)","Lake Ontario (west)")
names(supp.labs) <- c("CAY","CRO","ECAN","ERIE","ERIW","LHS","LMK","LMM",
                      "LSC","ONE","ONO","OSW","ROC")
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
# cols <- c(tol(n = 12), "#5E4FA2")

map <- ggplot() +
    geom_polygon(data=all_states, aes(x=long, y=lat, group=group), color="lightgray", fill="white") +
    geom_polygon(data=great_lakes, aes(x=long, y=lat, group=group), fill = "slategray") +
    geom_polygon(data=finger_lakes, aes(x=long, y=lat, group=group), fill = "slategray") +
    geom_polygon(data=oneida, aes(x=long, y=lat, group=group), fill = "slategray") +
    geom_path(data=df.usa_streams, aes(x=long, y=lat, group=group), color = "lightgray", size = 0.1) +
    geom_path(data=df.usa_canals, aes(x=long, y=lat, group=group), color = "slategray", size = 0.8) +
    geom_point(shape=21, size=5, color="grey23", data=sampling_points, aes(x=long, y=lat, fill=(Site_ID))) +
    scale_fill_manual(values=cols, labels=supp.labs,  name = "") +
    coord_map(projection="mercator", xlim = c(-88.5,-75), ylim = c(41,47)) +
    ggsn::scalebar(location = "topright", x.min = -75.5, x.max = -85, 
                   y.min = 43, y.max = 46.8, 
                   dist = 200, dist_unit="km", transform = TRUE,
                   st.dist = 0.1, height = 0.03, border.size=0.5) +
    theme_bw() +
    theme(panel.grid = element_blank() ,axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

map
ggsave(plot=map, filename="figures/goby_sampling_sites_summer_2019.pdf", width=7, height=5)
