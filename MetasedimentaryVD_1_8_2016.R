library(soilDB)
library(aqp)
library(plyr)
library(lattice)
library(maps)
library(raster)
library(rgdal)
library(reshape)
library(reshape2)
library(sharpshootR)
library(RColorBrewer)
library(latticeExtra)
library(plotrix)
library(classInt)
library(vegan)

#define plotting style 
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen', 'yellow'), lwd=2))

#fetch KSSL data by fuzzy-matching of series name
jocal <- fetchKSSL('Jocal')
josephine <- fetchKSSL('Josephine')
schaad <- fetchKSSL('Schaad')
nedsgulch <- fetchKSSL('Nedsgulch')

# standardize taxonname 
jocal$taxonname <- 'jocal'
josephine$taxonname <- 'josephine'
schaad$taxonname <- 'schaad'
nedsgulch$taxonname <- 'nedsgulch'



#combine data
x <- rbind(jocal, josephine, schaad, nedsgulch)

#subset points in MLRA 22A
x.22 <- subsetProfiles(x, s="mlra == '22A'")

groupedProfilePlot(x.22, groups = 'taxonname', label='pedon_id', color='bs82', name='hzn_desgn')


# Read CSV file with climate attributes
u <- read.csv(file = 'C:/Users/Andrew.Paolucci/Desktop/Jocal_Joseph_Neds_Scha_Report_1_12_2016/Jocal22_sitedata6.csv')
# Add climate data to site data
site(x.22) <- u

## Diagnostic Property Plot ##
# load nasis data
nasis <- fetchNASIS()
# subset nasis sites us labsampnum
nasis.lab <- which(nasis$pedlabsampnum %in% x.22$pedlabsampnum)
nasis.sub <- nasis[nasis.lab, ]

#iniitalize spatial coordinates
coordinates(nasis.sub) <- ~ x_std + y_std
# get/set spatial reference system
proj4string(nasis.sub) <- '+proj=longlat +datum=WGS84'

# tabulate number of pedons / generalized taxonname
(taxname.tab2 <- table(nasis.sub$taxonname))

#merge in climate and current taxon info
# Read CSV file with climate attributes
nsd <- read.csv(file = 'C:/Users/Andrew.Paolucci/Desktop/Jocal_Joseph_Neds_Scha_Report_1_12_2016/Jocal22_sitedata5.csv')

# Add climate data to site data
site(nasis.sub) <- nsd

# convert PSC to factor
nasis.sub$part_size_class <- factor(nasis.sub$part_size_class)
# work-around: new function in sharpshootR
hp <- multinominal2logical(nasis.sub, 'part_size_class')
site(nasis.sub) <- hp

# convert tax_order to factor
nasis.sub$Current_Order <- factor(nasis.sub$Current_Order)
# work-around: new function in sharpshootR
kp <- multinominal2logical(nasis.sub, 'Current_Order')
site(nasis.sub) <- kp
# convert current taxonname to factor
nasis.sub$Current_Taxoname <- factor(nasis.sub$Current_Taxoname)

# tabulate number of pedons / Current taxonname
(taxname.tab3 <- table(nasis.sub$Current_Taxoname))

#Get soil depth class
sdc <- getSoilDepthClass(nasis.sub, depth.classes = c(shallow = 50, mod.deep = 100, deep = 150, very.deep = 1000))
site(nasis.sub) <- sdc

# Diagnostic Property Plot
v <- c('Ultisol', 'Alfisols', 'argillic.horizon','deep', 'very.deep', 'fine', 'fine.loamy', 'coarse.loamy','Inceptisols')
diagnosticPropertyPlot(nasis.sub, v, k=5, dend.label = 'Current_Taxoname', grid.label='pedlabsampnum')






## Plot Locations of Sites ##
#iniitalize spatial coordinates of x.22
coordinates(x.22) <- ~ x + y
# get/set spatial reference system
proj4string(x.22) <- '+proj=longlat +datum=WGS84'

# extract spatial data + site level attributes
x.22.sp <- as(x.22, 'SpatialPointsDataFrame')

#graphic settings
par(mfrow=c(1,1), pin=c(7,7))
#generate a basemap of northern California with county boundaries
par(mar=c(0,0,0,0))
map('county', 'California', xlim=c(-122, -119), ylim=c(37, 40))
#add long/lat axes
map.axes()
#add locations of lab sites
points(y ~ x, data=site(jocal), pch=21, bg='Red')
points(y ~ x, data=site(josephine), pch=21, bg='Blue')
points(y ~ x, data=site(nedsgulch), pch=21, bg='Yellow')
points(y ~ x, data=site(schaad), pch=21, bg='Green')
# add title
title('Location of lab Sites in MLRA 22A', cex.main=1.5, line=5)
# add a simple legend
legend('topright', pch=21, pt.bg=c('Red','Blue', 'Yellow','Green'), legend=c('Jocal','Josephine', 'Nedsgulch','Schaad'))
# save to shapefile
writeOGR(x.22.sp, dsn='C:/Paolucci/GIS/CA630', layer='VD_MSD_BaseSatFigure', driver='ESRI Shapefile', overwrite_layer=TRUE)

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
xyplot(bs82.100 ~ pH.100, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='8.2 Base Saturation (%)', xlab='pH', xlim=c(6.0, 4.75), ylim=c(0, 70))
# water balance
xyplot(bs82.100 ~ waterbalance, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='8.2 Base Saturation (%)', xlab='waterbalance')
# precip
xyplot(bs82.100 ~ Precipitation, data=x.22.site.data, grid=TRUE, type= c("p", "smooth"), ylab='8.2 Base Saturation (%)', xlab='Precip (mm)')
# water blance related to precip
xyplot(waterbalance ~ Precipitation, data=x.22.site.data, grid=TRUE, type= c("p", "r"), ylab='Water Balance', xlab='Precip (mm)')



