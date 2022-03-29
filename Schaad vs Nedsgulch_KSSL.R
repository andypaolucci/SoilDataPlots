library(soilDB)
library(aqp)
library(plyr)
library(lattice)
library(latticeExtra)
library(maps)
library(MASS)
library(rgdal)
library(raster)
library(reshape2)
library(sharpshootR)



#define plotting style 
tps <- list(add.line=list(col=c('DarkRed'), lwd=1))

#fetch KSSL data by fuzzy-matching of series name
s <- fetchKSSL(series='Antigo', returnMorphologicData = TRUE, simplifyColors = TRUE)

#extract pedons into SPC
pedons <- s$SPC

#Check out list
str(t,2)

# explicit string matching
idx <- which(pedons$state == 'Wisconsin')
pedons <- pedons[idx, ]

pedons$taxonname <- 'Antigo'

#Plot
par(mar=c(1,1,10,2))

plotSPC(pedons, color="moist_soil_color", name='hzn_desgn', groups = 'taxonname', label='pedon_id', alt.label='slope_field', alt.label.col='green', 
                   print.id = TRUE, id.style = 'side', cex.id=0.9, cex.names=0.8, axis.line.offset=-2.5)
title('Antigo Series', cex.main=2, line= 3)
title('Moist Soil Colors', cex.main=1, line= 1)

plotSPC(pedons, color="ph_h2o", name='hzn_desgn', groups = 'taxonname', label='pedon_id', alt.label='slope_field', alt.label.col='green', 
        print.id = TRUE, id.style = 'side', cex.id=0.9, cex.names=0.9, axis.line.offset=-2.5, col.label='Soil pH',  col.legend.cex = 1.1)
title('Antigo Series', cex.main=1.25, line= 4)



groupedProfilePlot(x, label = 'suborder', groups='taxonname', print.id = TRUE, id.style = 'side', name='hzn_desgn', color='clay', cex.id=0.6, y.offset=0.3, col.label='')


#generate a basemap of northern California with county boundaries
map <- map('county', 'california', xlim=c(-123.25, -118.75), ylim=c(36.5, 41))

x.sp <- SpatialPoints(x)

#Subset Lab Pedons in MLRA 22A 
#Read MLRA Shapefile
MLRA22A <- readOGR('L:/NRCS/MLRAShared/Geodata/mlra/MLRA22A.shp', layer ='MLRA22A') 

#add long/lat axes
map.axes()

#add locations of lab sites
points(y ~ x, data=site(x), pch=21, bg='DarkRed')

# add a simple legend
legend('topright', pch=21, pt.bg=c('DarkRed'), legend=c('sites'), bty='n')


# plot 
pdf()
par(mar=c(1,1,7,1))
plotSPC(x, label = 'generalized_taxonname', print.id = TRUE, id.style = 'side', name = 'hzn_desgn', color='bs82', cex.id=0.7, y.offset=0, col.label='', scaling.factor=1, width=0.15, axis.line.offset = -7)
title('bs82', cex.main=3, line=5)

# make column to store order and fill with values 
x$order <- rep(NA, times=length(x))
x$order <- c(17,16,1,7,15,8,14,12,9,13,10,2,3,4,11,5,6)
head(x)

# plot in order color = BS
par(mar=c(1,1,7,1))
plotSPC(x, label = 'generalized_taxonname', print.id = TRUE, id.style = 'side', name = 'hzn_desgn', color='bs82', cex.id=0.6, y.offset=0.3, col.label='', plot.order=order(x$order))
title('bs82', cex.main=1)

# plot color = pH
par(mar=c(1,1,7,1))
plotSPC(x, label = 'site_id', print.id = TRUE, id.style = 'side', name = 'hzn_desgn', color='ph_h2o', cex.id=0.6, y.offset=0.3, col.label='')
title('pH 1:1 H2O', cex.main=1)

# find missing pH values
x$ph_h2o


# fill in missing pH values with data from NASIS. Use figure to find pedon IDs
x$ph_h2o <- c(6.1, 6.4, 6.4, 6.3, 6.4, 6.6, NA, 6.7, 6.7, 6.1, 6.4, 6.5, 6.5, 6.8, NA, NA, 6.7, 5.7, 5.6, 5.7, 5.8, 5.9, 5.9, 5.7, 6.3, 6.1, 6.1, 6.0, 6.0, 5.9, 6.0, 6.1, 6.0, 5.9, 6.0, 6.0, 6.5, 7.3, 7.4, 7.3, 5.9, 5.9, 5.7, 5.7, 5.7, 5.5, 5.6, 5.5, 5.6, 5.2, 5.3, 5.3, 5.4, 5.6, 5.0, 6.0, 5.5, 5.5, 5.6, 5.5, 6.2, 4.8, 6.0, 6.3)

# plot color = pH
par(mar=c(1,1,7,1))
plotSPC(x, label = 'generalized_taxonname', print.id = TRUE, id.style = 'side', name = 'hzn_desgn', color='ph_h2o', cex.id=0.6, y.offset=0.3, col.label='', plot.order=order(x$order))
title('pH 1:1 H2O', cex.main=1)

# xyplot
xyplot(bs82 ~ ph_h2o, data=horizons(x), type= c("p","r","g"), main='title') + latticeExtra::layer(panel.abline(v=6, h=35, col='red', lty=2))

# robust linear model, not discriminating enough
l <- rlm(bs82 ~ ph_h2o, data=horizons(x))

xyplot(bs82 ~ ph_h2o, data=horizons(x), type= c("p","r","g"), main='title') + latticeExtra::layer(panel.abline(v=6, h=35, col='red', lty=2)) + latticeExtra::layer(panel.abline(coef(l), col='DarkGreen', lty=2))

## rms package
library(rms)
# include non-linearity into a linear model framework
h <- horizons(x)[, c('bs82', 'ph_h2o')]
dd <- datadist(h)
options(datadist="dd")
(l <- ols(bs82 ~ ph_h2o, data=h))
xyplot(bs82 ~ ph_h2o, data=horizons(x)) + plot(Predict(l))

# use RCS to deal with non-linearity
(l.2 <- ols(bs82 ~ rcs(ph_h2o), data=h))
xyplot(bs82 ~ ph_h2o, data=horizons(x)) + plot(Predict(l.2))

plot(nomogram(l.2))
