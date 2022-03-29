library(soilDB)
library(aqp)
library(plyr)
library(lattice)
library(maps)
library(raster)
library(rgdal)
library(reshape)
library(reshape2)

#define plotting style 
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen'), lwd=2))

#Load KSSL data using large bounding box
# Had to do two seperate bounding boxes 
kssl.1 <- fetchKSSL(bbox=c(-116, 31, -125, 39))
kssl.2 <- fetchKSSL(bbox=c(-116, 39, -125, 42))
kssl.3 <- fetchKSSL(bbox=c(-116, 42, -125, 45))
kssl.4 <- fetchKSSL(bbox=c(-116, 45, -125, 49))

#combine GAVE ERROR
kssl.california <- rbind(kssl.1, kssl.2)

#define coordinates/system for group 1
coordinates(kssl.1)<- ~ x + y
proj4string(kssl.1) <- '+proj=longlat +datum=WGS84'
kssl.1.sp <- as(kssl.1, "SpatialPointsDataFrame")

#define coordinates/system for group 2
coordinates(kssl.2)<- ~ x + y
proj4string(kssl.2) <- '+proj=longlat +datum=WGS84'
kssl.2.sp <- as(kssl.2, "SpatialPointsDataFrame")

#define coordinates/system for group 3
coordinates(kssl.3)<- ~ x + y
proj4string(kssl.3) <- '+proj=longlat +datum=WGS84'
kssl.3.sp <- as(kssl.3, "SpatialPointsDataFrame")

#define coordinates/system for group 4
coordinates(kssl.4)<- ~ x + y
proj4string(kssl.4) <- '+proj=longlat +datum=WGS84'
kssl.4.sp <- as(kssl.4, "SpatialPointsDataFrame")

#combine
kssl.california.sp <-rbind(kssl.1.sp, kssl.2.sp)
kssl.northcoast.sp <-rbind(kssl.3.sp, kssl.4.sp)

#Map
map <- map('state', xlim=c(-125, -116), ylim=c(31, 50))
points(kssl.california.sp)
points(kssl.northcoast.sp)


# save to shapefile
writeOGR(kssl.northcoast.sp, dsn='C:/Users/Andrew.Paolucci/Desktop/CZSS', layer='North_Coast_kssl_3_19_2018AJP2', driver='ESRI Shapefile', overwrite_layer=TRUE)




# tabulate number of pedons / generalized taxonname
(pedon.tab <- table(kssl$taxonname))



# aggregate properties by generalized taxonname 
a <- slab(x.22, taxonname ~ clay + bs82 + ph_h2o)

# compute mid-point between slice top and bottom depth for plotting
a$mid <- with(a, (top+bottom)/2)

# append the number of pedons to the taxonname label
a$taxonname <- factor(a$taxonname, levels=names(taxname.tab), labels=paste(names(taxname.tab), ' (', taxname.tab, ')', sep=''))

levels(a$variable)
levels(a$taxonname)

#re-name soil property labels
levels(a$variable) <- c('Clay (%)', 'Base Saturation at pH 8.2 (%)', 'pH 1:1 H2O')

#plot aggregate data

xyplot(top ~ p.q50 | variable, groups=taxonname, data=a, ylab='Depth',
       xlab='median bounded by 5th and 95th percentiles',
       lower=a$p.q25, upper=a$p.q75, ylim=c(155,-5),
       panel=panel.depth_function, alpha=0.25, sync.colors=TRUE,
       prepanel=prepanel.depth_function,
       cf=a$contributing_fraction,
       strip=strip.custom(bg=grey(0.85)),
       layout=c(3,1), scales=list(x=list(alternating=1, relation='free'), y=list(alternating=3)), par.settings=tps, auto.key=list(columns=2, lines=TRUE, points=FALSE))


# extract BS8.2 and BS7 at select depths
slices.x.22 <- slice(x.22, c(50, 100, 150) ~ bs7 + bs82 + ph_h2o + cec82 + clay + oc, strict = FALSE, just.the.data = TRUE)

# reshape slices
slices.long <- melt(slices.x.22, id.vars=c('pedon_key', 'hzn_top'), measure.vars=c('bs7', 'bs82', 'ph_h2o','cec82','clay','oc'))
slices.7 <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'bs7')
slices.82 <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'bs82')
slices.pH <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'ph_h2o')
slices.cec <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'cec82')
slices.clay <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'clay')
slices.oc <- dcast(slices.x.22, pedon_key ~ hzn_top, value.var = 'oc')

names(slices.7)[2:4] <- paste0('bs7.', c(50, 100, 150))
names(slices.82)[2:4] <- paste0('bs82.', c(50, 100, 150))
names(slices.pH)[2:4] <- paste0('pH.', c(50, 100, 150))
names(slices.cec)[2:4] <- paste0('cec82.', c(50, 100, 150))
names(slices.clay)[2:4] <- paste0('clay.', c(50, 100, 150))
names(slices.oc)[2:4] <- paste0('oc.', c(50, 100, 150))


# merge slices into lab site data
site(x.22) <- slices.7
site(x.22) <- slices.82
site(x.22) <- slices.pH
site(x.22) <- slices.cec
site(x.22) <- slices.clay
site(x.22) <- slices.oc

str(x.22)
x.22.sp<- as(x.22, 'SpatialPointsDataFrame')
str(x.22.sp)
# subset columns of lab data
x.22.s1 <- x.22.sp[, c('pedon_id', 'taxonname', 'waterbalance', 'Elevation', 'Precipitation', 'Air_temp', 'bs82.50', 'bs82.100', 'bs82.150', 'bs7.50', 'bs7.100', 'bs7.150', 'pH.50', 'pH.100', 'pH.150', 'cec82.50', 'cec82.100', 'cec82.150', 'clay.50', 'clay.100', 'clay.150', 'oc.50', 'oc.100', 'oc.150')]
# keep only the attribute data
x.22.site.data <- x.22.s1@data
str(x.22.site.data)



# plot relationships
# pH
xyplot(bs82.100 ~ pH.100, data=x.22.site.data, grid=TRUE, type= c("p", "smooth"), ylab='8.2 Base Saturation (%)', xlab='pH', xlim=c(6.75, 4.75), ylim=c(15, 85))
# water balance
xyplot(bs82.100 ~ waterbalance, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='8.2 Base Saturation (%)', xlab='waterbalance')
# precip
xyplot(bs82.100 ~ Precipitation, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='8.2 Base Saturation (%)', xlab='Precip (mm)')
# water blance related to precip
xyplot(waterbalance ~ Precipitation, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='Water Balance', xlab='Precip (mm)')

# elevation
xyplot(bs82.100 ~ Elevation, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='8.2 Base Saturation (%)', xlab='Elevation (m)')

cloud(bs82.100 ~ waterbalance + clay.100, data=x.22.site.data, auto.key=list(columns=3, points=TRUE, lines=FALSE))


