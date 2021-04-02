### August 4, 2019
### Make a map of sampling sites in the Great Lakes/ Finger lakes region

setwd("Downloads/") # set the working directory to Downloads

# install and load all of the packages we will be using
install.packages(c("raster", "maptools", "rgdal", "maps", "usmap"))
library(raster) 
library(maptools)
library(rgdal)
library(maps)
library(usmap)

# Go to the website to download lakes data (shapefile):
# https://www.naturalearthdata.com/downloads/10m-physical-vectors/
# scroll down to Lakes + Reservoirs
# click "Download lakes" and save zip file to your Downloads
# click on the zip file to unzip it (should be a file named ne_10m_lakes)

# Make objects for Canada and US 
canada <- getData("GADM", country = "CAN", level = 1) # level = 1 provides province data
usa <- getData("GADM", country = "USA", level = 1)

# make lake objects 
lakes <- readOGR("ne_10m_lakes") # file path to lake shapefile
oneida <- subset(lakes, name=="Oneida Lake") # make object for Oneida Lake
finger_lakes <- subset(lakes, name_alt=="Finger Lakes") # make object for Finger Lakes
great_lakes <- subset(lakes, name_alt=="Great Lakes") # make object for Great Lakes

# make the map and save it in your Downloads
png(filename = "goby_sampling_sites_summer_2019.png", width = 10, height = 7.5, res = 300, units = "in")
plot(canada, border = "gray55", xlim = c(-89, -75), ylim = c(40.5, 47)) # plot Canada
plot(usa, border = "gray55", add = TRUE) # plot USA
plot(oneida, col = "slategray1", add = TRUE) # add Oneida Lake
plot(finger_lakes, col = "slategray1", add = TRUE) # add Finger Lakes
plot(great_lakes, col = "slategray1", add = TRUE) # add Great Lakes 

# make a data frame of sampling points longitude/latitude
sampling_points <- data.frame(x = c(-82.7839, -82.5390, -86.3392, -76.951418, -76.750299, -76.484041, -76.241005, -83.181959, -79.595408, -87.8212, -76.005611, -78.8784, -76.5105, -77.6453, -77.6069),
                      y = c(42.5642, 43.4309, 43.2235, 42.873401, 42.897155, 43.149889, 43.115009, 42.065604, 42.341587, 42.5847, 43.179738, 42.8864, 43.4553, 43.1167, 43.2595))
# make a data frame of future sampling sites longitude/latitude
# future_sites <- data.frame(x = c(-76.005611, -76.939996, -78.878974, -75.2327),
#                            y = c(43.179738, 43.238647, 43.305242, 43.1009))

# add sampling points to map
points(x = sampling_points$x, y = sampling_points$y, pch = 19, cex = 1.7, col = "black")
# points(x = future_sites$x, y = future_sites$y, pch = 19, cex = 1.7, col = "black")
points(x = sampling_points$x, y = sampling_points$y, pch = 19, cex = 1.5, col = "darkorange")
# points(x = future_sites$x, y = future_sites$y, pch = 19, cex = 1.5, col = "purple")

map.scale(y=41, x=-89, ratio=FALSE, relwidth = 0.3) # add a scale bar
box() # add a box around the map
dev.off() # save the map
