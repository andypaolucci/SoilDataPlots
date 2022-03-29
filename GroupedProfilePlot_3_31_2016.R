# load required packages
library(lattice)
library(soilDB)
library(aqp)
library(rgdal)
library(plyr)
library(sp)
library(cluster)
library(ape)
library(RColorBrewer)
library(latticeExtra)
library(plotrix)
library(classInt)
library(sharpshootR)

# setup plotting style for later
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen'), lwd=2))

# Get map unit data
mu <- readOGR(dsn = 'C:/Paolucci/Projects/CA628/soilmu_a_ca628.shp', stringsAsFactors=FALSE)

# get pedon data
x <- fetchNASIS()

str(x)

plotSPC(x, plot.order=order(x$tsectstopnum), print.id=TRUE, label='tsectstopnum', color='dry_soil_color')

str(x)

par(mar=c(1,1,5,1))
groupedProfilePlot(x, groups = 'hillslopeprof', alt.label='taxonname', alt.label.col='black', label='tsectstopnum', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.6, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, color='clay')
title('Gravels', cex.main=1)
abline(h=c(50, 100, 150), lty=2, col='red')

par(mar=c(1,1,5,1))
groupedProfilePlot(x, groups = 'Age', alt.label='taxonname', alt.label.col='black', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.6, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, color='clay')
title('7207', cex.main=1)

# keep only those pedons with real coordinates
good.idx <- which(!is.na(x$x_std))
x <- x[good.idx, ]

#initalize spatial coordinates
coordinates(x) <- ~ x_std + y_std

# get/set spatial reference system
proj4string(x) <- '+proj=longlat +datum=WGS84'

# extract spatial data + site level attributes
# transform to CRS of mu data
x.sub <- as(x, 'SpatialPointsDataFrame')
x.sub <- spTransform(x.sub, CRS(proj4string(mu)))

# subset seleted map unit
idx <- grep('7212', mu$MUSYM, ignore.case = TRUE)
s.sub <- mu[idx , ]

# copy MUSYM to pedons using spatial overlay
x.sub$musym <- over(x.sub, s.sub)$MUSYM

# Graph
par(mar=c(0,0,0,0))
plot(s.sub)
points(x.sub, col=2, cex=0.5)

#remove NAs
idx <- which(!is.na(x.sub$musym))
y <- x.sub[idx, ]

# get profiles in x whose peiid are in subset y
idx <- which(profile_id(x) %in% y$peiid)
x.mu <- x[idx, ]


# SDC
sdc <- getSoilDepthClass(x.mu)
site(x.mu) <- sdc
str(x.mu)

# Plot ordering by depth
par(mar=c(1,1,5,1))
plotSPC(x.mu, plot.order=order(x.mu$slope_field), label='pedon_id', cex.id=0.7, cex.names=0.7, col.label=' ', width=0.15, y.offset=0, id.style='side', x.idx.offset=0.2, shrink=FALSE, shrink.cutoff=1)
title('7212 Grouped Profile Plot')

# set up classes using field aspect
plotvar <- x.mu$aspect_field
class <- classIntervals(plotvar, 5, style = "fixed", fixedBreaks = c(0, 45, 135, 225, 315, 360))
#extract classes
onion <-findCols(class)
#add to site data
x.mu$aspect_class <-onion

# make a new column to store aspect class and fill with NA
x.mu$generalized_aspect <- rep(NA, times=length(x.mu))

# Sort classes
x.mu$generalized_aspect[x.mu$aspect_class %in% c('1','5')] <- 'North'
x.mu$generalized_aspect[x.mu$aspect_class %in% c('2')] <- 'East'
x.mu$generalized_aspect[x.mu$aspect_class %in% c('3')] <- 'South'
x.mu$generalized_aspect[x.mu$aspect_class %in% c('4')] <- 'West'

# subset only those pedons with generalized bedrock
x.mu <- subsetProfiles(x.mu, s="!is.na(generalized_aspect)")

# make a new column to store aspect class and fill with NA
x.mu$MU_aspect <- rep(NA, times=length(x.mu))

# Sort classes
x.mu$MU_aspect[x.mu$generalized_aspect %in% c('North','East')] <- 'North/East'
x.mu$MU_aspect[x.mu$generalized_aspect %in% c('South', 'West')] <- 'South/West'

# subset only those pedons with generalized bedrock
x.mu <- subsetProfiles(x.mu, s="!is.na(MU_aspect)")

#Check grouping with 
groupedProfilePlot(x.mu, groups = 'taxonname', label='pedon_id', color='clay', alt.label='ecositeid', alt.label.col='black', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.5, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, col.label='')
title('7212 Grouped Profile Plot: Color=Clay')


# check sorted taxonnames
(taxname.tab <- table(x.mu$generalized_aspect))

# extract spatial data + site level attributes
x.mu.sp <- as(x.mu, 'SpatialPointsDataFrame')

#Alternatively Load Shapefile
MU3046<-readOGR(dsn='S:/NRCS/Archive_Andy_Paolucci/CA630/ComponentReports/3046/3046_SuesTransectNotInNasis_7_5_2016.shp', layer='3046_SuesTransectNotInNasis_7_5_2016')


# Load curvature class raster
curvature <- raster('l:/Geodata/DEM_derived/curvature_classes_15.tif')
library(raster)
elevation <- raster('L:/Geodata/elevation/10_meter/ca630_elev')


# sample raster using series data points
x.mu.elev <- extract(elevation, MU3046)

# Add values back to pedon data??
MU3046$elevation <- x.mu.elev

head(MU3046)
table(MU3046$elevation, MU3046$IDENT)

MU3046.df<- as.data.frame(MU3046, row.names = NULL)

str(plotclr)
  # set up colors for Caly comparison
  plotclr <- brewer.pal(3, "PuBuGn")
  plotvar <- x.mu$clay
  class <- classIntervals(plotvar, 3, style = "fixed", fixedBreaks = c(1, 18, 35, 55))
  colcode <- findColours(class, plotclr)
  x.mu$colcode <- colcode
  class.names <- names(attr(colcode, "table"))


# Plot 2  Grouped Profile plot
par(mar=c(1,1,3,1))
str(x.mu)
  #remove NAs
  groupedProfilePlot(x.mu, groups = 'geompos_hill', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.4, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, col.label='', color='clay')
  groupedProfilePlot(x.mu, groups = 'generalized_aspect', label='pedon_id', alt.label='aspect_field', alt.label.col='green', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.6, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, col.label='')
  title('7210 Grouped Profile Plot')

par(mar=c(1,1,5,1))
groupedProfilePlot(x.mu, groups = 'ecositeid', alt.label='taxonname', alt.label.col='green', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.5, cex.names=0.5, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, col.label='')
title('8178:Ecosite Group Plot', cex.main=1)

# Do not run
axis(3, line=-30, at=1:length(x.mu), labels=x.mu$taxon_kind, tick=FALSE, las=2, cex.axis=0.8)
axis(3, line=-40, at=1:length(x.mu), labels=x.mu$pedon_id, tick=TRUE, las=3, cex.axis=0.8)




legend("bottom", horiz=TRUE, legend = c("0 - 18", "18 - 35", "35 - 55"), fill=plotclr, cex=1)


(sort(table(x.mu$taxonname), decreasing = TRUE))

# make a new column to store generalized taxonname, and fill with NA
x.mu$generalized_taxname <- rep(NA, times=length(x.mu))

# Sort pedons 
x.mu$generalized_taxname[x.mu$taxonname %in% c('Nedsgulch')] <- 'Nedsgulch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Sheepranch')] <- 'Sheepranch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Wallyhill')] <- 'Wallyhill'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Fricot')] <- 'Fricot'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Ultic Haploxeralfs', 'Loamy-skeletal, mixed, super, acid, thermic Lithic Xerorthen')] <- 'Other'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Arpatutu')] <- 'Arpatutu'

# subset only those pedons with generalized bedrock
x.mu <- subsetProfiles(x.mu, s="!is.na(generalized_taxname)")

# check sorted taxonnames
(taxname.tab <- table(x.mu$generalized_taxname))

groupedProfilePlot(x.mu, groups = 'generalized_taxname', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.4, cex.names=0.4, y.offset=0, group.name.offset=-8, group.name.cex=0.5, width=0.15, group.line.col='red', group.line.lwd=2, shrink=FALSE, shrink.cutoff=1)
title('6038 Taxonomic Group Profile Plot', cex.main=1)

plotSPC(x.mu, plot.order=order(x.mu$bedrckdepth), label='pedon_id', print.id = 
        TRUE, id.style = 'side', cex.id=0.4, cex.names=0.4, y.offset=0, 
        group.name.offset=-8, group.name.cex=0.5, width=0.15, group.line.col='red', 
        group.line.lwd=2, shrink=FALSE, shrink.cutoff=1)



# plot 4 hillslope position grouped profile plot
par(mar=c(2,1,3,1))
groupedProfilePlot(x.mu, groups = 'hillslope_pos', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.3, y.offset=0.4, group.name.offset=-14, group.name.cex=0.4, width=0.2, group.line.col='red')
title('8034 Taxonomic Group Profile Plot', cex.main=1)

# if NA make table to show hillslope/taxonomic relationship
table(x.mu$generalized_taxname, x.mu$hillslope_pos)

write.csv(pedons, file='S:/NRCS/Archive_Andy_Paolucci/CA630pedonlist_6_7_2016.csv', row.names=FALSE)


