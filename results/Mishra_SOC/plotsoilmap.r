# Load packages
library(maps)
library(maptools)
library(mapproj)
library(RColorBrewer)
library(classInt)
library(gpclib)
library(mapdata)

# Define vector with the values that you would like to see plotted at desired lat/long. Your csv input file loaded as dataframe (Var) must feature the following columns (Site is optional, but useful for labeling) ID,Para,Lat,Lon
#Var = read.csv('Site_info.csv')
Var = read.csv('allsite.csv')
site = read.csv('Ameriflux_info.csv')
#plotvar <- Var$Landcover
plotvar <- Var$Landcover[Var$Soil.Order == 'Gelisol']
plotvar2 <- site$Landcover[site$Soil.Order == 'Gelisol']

# Define number of colours to be used in plot
nclr <- 7

# Define colour palette to be used
plotclr <- brewer.pal(nclr,"RdPu")

# Define colour intervals and colour code variable for plotting
class <- classIntervals(plotvar, nclr, style = "pretty")
colcode <- findColours(class, plotclr)
#colcode <- findColours(Var$Suborder[Var$Soil.Order == 'Gelisol'], plotclr)

class2 <- classIntervals(plotvar2, nclr, style = "pretty")
colcode2 <- findColours(class2, plotclr)

# Plot the map with desired lat/long coordinates and data points with colour coding and legend
#map('world2','USA:alaska')
#map("worldHires", xlim=c(-170,-120), ylim=c(50,72), col='gray90', fill=TRUE)
m <- map("worldHires", xlim=c(-180,-135), ylim=c(55,72), col='gray90', fill=TRUE)
map.grid(m)#, lim = c(-185,-130, 55, 75), nx=3, ny=3)#, c(-180, -135, 55, 72))
#points(Var$Lon, Var$Lat, pch = 16, col= colcode, cex = 2)
points(Var$Lon, Var$Lat, pch = 16, col= Var$Suborder[Var$Soil.Order == 'Gelisol'], cex = 2)
# text(Var$Lon, Var$Lat, labels = Var$Profile.ID, pos = 4)
# pointLabel(Var$Lon, Var$Lat, labels = Var$Profile.ID, method = c("SANN", "GA"), allowSmallOverlap = TRUE, cex = 1)
#points(site$Lon, site$Lat, pch = 8, col= colcode2, cex = 2)
#text(site$Lon, site$Lat, labels = site$Profile.ID, pos = 4)
legend("bottomright", legend = c('Orthel','Turbel'), fill=1:length(Var$Suborder[Var$Soil.Order == 'Gelisol']), cex = 1.5, bty = "y")
#legend("bottomright", legend = names(attr(colcode, "table")), fill = attr(colcode, "palette"), cex = 0.7, bty = "n")
