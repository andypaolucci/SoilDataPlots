install.packages("sharpshootR", repos="http://R-Forge.R-project.org")
install.packages("soilDB", repos="http://R-Forge.R-project.org")

# setup selected set in NASIS

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
library(raster)
# setup plotting style for later
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen'), lwd=2))

# Get map unit data
mu <- readOGR(dsn = 'L:/NRCS/MLRAShared/CA630/FG_CA630_OFFICIAL.gdb', layer = 'ca630_a', stringsAsFactors=FALSE)

# get pedon data
x <- fetchNASIS()

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
idx <- grep('3058', mu$MUSYM, ignore.case = TRUE)
s.sub <- mu[idx , ]

# copy MUSYM to pedons using spatial overlay
x.sub$musym <- over(x.sub, s.sub)$MUSYM

#remove NAs
idx <- which(!is.na(x.sub$musym))
y <- x.sub[idx, ]

# get profiles in x whose peiid are in subset y
idx <- which(profile_id(x) %in% y$peiid)
x.mu <- x[idx, ]

# add MU Id to site info
site(x.mu) <- x.sub$musym


# Plot 2  Grouped Profile plot
par(mar=c(1,1,3,1))
groupedProfilePlot(x.mu, groups = 'ecositeid', label='pedon_id', print.id = TRUE, id.style = 'side', cex.id=0.6, cex.names=0.6, y.offset=7, axis.line.offset=-2.5, group.name.cex=0.5, group.name.offset = -6, col.label='', color='clay')
plotSPC(x.mu, plot.order=order(x.mu$slope_field), label='pedon_id', cex.id=0.8, cex.names=0.8, color='clay', col.label=' ', width=0.15, y.offset=0, id.style='side', x.idx.offset=0.2, shrink=FALSE, shrink.cutoff=1)

# regression graph
site.x.mu <- site(x.mu)
xyplot(x.mu$depth ~ slope_field, data=site.x.mu, grid=TRUE, type= c("p", "r"), ylab='Depth (cm)', xlab='Slope %', xlim=c(0, 90), ylim=c(25, 150))
head(site.x.mu)
histogram(x.mu$depth, breaks=3)
histogram(x.mu$slope_field, breaks=3)


# tabulate number of pedons in each taxon
(sort(table(x.mu$taxonname), decreasing = TRUE))

# SDC
sdc <- getSoilDepthClass(x.mu)
site(x.mu) <- sdc

#PSC
# convert PSC to factor
x.mu$part_size_class <- factor(x.mu$part_size_class)
# work-around: new function in sharpshootR
hp <- multinominal2logical(x.mu, 'part_size_class')
site(x.mu) <- hp

#Hillslope Position
# work-around: new function in sharpshootR
hp2 <- multinominal2logical(x.mu, 'hillslope_pos')
site(x.mu) <- hp2

#Slope Position
x.mu$slope_position <- factor(x.mu$slope_position)
hp3 <- multinominal2logical(x.mu, 'slope_position')
site(x.mu) <- hp3

#Shape Across
hp4 <- multinominal2logical(x.mu, 'shapeacross')
site(x.mu) <- hp4

#Geomorphic Position (3D)
hp5 <- multinominal2logical(x.mu, 'geompos_hill')
site(x.mu) <- hp5

# properties to consider, no need to convert to factors
v <- c('Concave', 'Linear', 'Convex')
# Diagnostic Plot
x <- diagnosticPropertyPlot(x.mu, v, k=4, grid.label='pedon_id', dend.label = 'taxonname')
x.mu$musym

#Shapedown
hp6 <- multinominal2logical(x.mu, 'shapedown')

str(x.mu)

## taxonname 
this.data <- xtabs(~ geomorphons + taxonname, data=site(x.mu))
this.data[this.data == 0] <- NA
# col.palette <- colorRampPalette(cols, bias=1/skewness(as.vector(this.data), na.rm = TRUE))

levelplot(this.data, col.regions=col.palette, colorkey=list(tick.number=15), xlab='10m Geomorphons', ylab='Field Described Hillslope Position', main='CA630:8172', scales=list(alternating=3), panel=function(x, y, z, ...) {
  panel.levelplot(x, y, z, ...)
  idx <- which(!is.na(z))
  panel.text(x[idx], y[idx], z[idx], font=2)
  panel.abline(h=seq(from=0.5, to=length(y), by=1), col=grey(0.45))
  panel.abline(v=seq(from=0.5, to=length(x), by=1), col=grey(0.45))
})




