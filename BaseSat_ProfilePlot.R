library(soilDB)
library(aqp)
library(plyr)
library(lattice)
library(maps)
library(raster)
library(rgdal)
library(sharpshootR)
library(RColorBrewer)
library(classInt)
library(sp)
library(cluster)
library(ape)
library(latticeExtra)
library(stats)
library(plotrix)
library(dismo)

# load extent data and colors
cols <-brewer.pal('Set1', n=3)
jocal <- seriesExtent('jocal')
josephine <- seriesExtent('josephine')
#margins
par(mar=c(0,0,0,0))
#map
map('county', 'California', xlim=c(-122, -119), ylim=c(37, 40))
#add extent data
plot(jocal, col=cols[1], add=TRUE)
plot(josephine, col=cols[2], add=TRUE)
#draw box
box()
map.axes()
#add legend
legend('topright', col=cols, pch=15, legend= c('Jocal', 'Josephine'))


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

plot(jocal, label='pedon_id', color='bs82', name='hzn_desgn')
# set up colors for BS comparison of x.22
plotclr <- brewer.pal(3, "PuBuGn")
plotvar <- jocal$bs82
class <- classIntervals(plotvar, 3, style = "fixed", fixedBreaks = c(0, 35, 75, 100))
colcode <- findColours(class, plotclr)
jocal$colcode <- colcode
class.names <- names(attr(colcode, "table"))
plot(jocal, label='pedon_id', color='colcode', name='hzn_desgn')

plot(josephine, label='pedon_id', color='bs82', name='hzn_desgn')
#subset points in MLRA 22A
josephine2<- subsetProfiles(josephine, s="mlra == '22A'")
plot(josephine2, label='pedon_id', color='bs82', name='hzn_desgn')
# set up colors for BS comparison of x.22
plotclr2 <- brewer.pal(3, "PuBuGn")
plotvar2 <- josephine2$bs82
class2 <- classIntervals(plotvar2, 3, style = "fixed", fixedBreaks = c(0, 35, 75, 100))
colcode2 <- findColours(class2, plotclr2)
josephine2$colcode <- colcode2
class.names2 <- names(attr(colcode, "table"))
plot(josephine2, label='pedon_id', color='colcode', name='hzn_desgn')

plot(schaad, label='pedon_id', color='bs82', name='hzn_desgn')
# set up colors for BS comparison of x.22
plotclr3 <- brewer.pal(3, "PuBuGn")
plotvar3 <- schaad$bs82
class3 <- classIntervals(plotvar3, 3, style = "fixed", fixedBreaks = c(0, 35, 75, 100))
colcode3 <- findColours(class3, plotclr3)
schaad$colcode <- colcode3
class.names3 <- names(attr(colcode, "table"))
plot(schaad, label='pedon_id', color='colcode', name='hzn_desgn')

plot(nedsgulch, label='pedon_id', color='bs82', name='hzn_desgn')
# set up colors for BS comparison of x.22
plotclr4 <- brewer.pal(3, "PuBuGn")
plotvar4 <- nedsgulch$bs82
class4 <- classIntervals(plotvar4, 3, style = "fixed", fixedBreaks = c(0, 35, 75, 100))
colcode4 <- findColours(class4, plotclr4)
nedsgulch$colcode <- colcode4
class.names4 <- names(attr(colcode, "table"))
plot(nedsgulch, label='pedon_id', color='colcode', name='hzn_desgn')

#combine data
x <- rbind(jocal, josephine2, schaad, nedsgulch)

x.22 <- subsetProfiles(x, s="mlra == '22A'")

# Read CSV file with climate attributes
u <- read.csv(file = 'C:/Paolucci/GIS/Jocal22_sitedata6.csv')
# Add climate data to site data
site(x.22) <- u

# Plot Grouped Profile plot
par(mar=c(1,1,7,1), mfrow=c(1,1))
groupedProfilePlot(x.22, groups='Current_Taxoname', group.name.offset=-20, label='Current_Order', color='colcode', 
                   name='hzn_desgn', id.style ='top', abbr=TRUE, abbr.cutoff=3, y.offset =5, x.idx.offset = 0, shrink=TRUE, shrink.cutoff=3)
title('8.2 Base Saturation Profile Plot ', cex.main=1)
legend("bottom", horiz=TRUE, legend = c("0-35%", "35-75%"), fill=plotclr, cex=1)

new.order <- order(x.22$Precipitation, decreasing =FALSE)

labelsmm  <- c(842.82, 847.40, 851.21, 900.14, 908.49, 913.63, 954.45, 959.07, 983.79, 1000.60, 1008.58, 1075.61, 1156.98, 1181.54, 1185.69, 1340.28, 1949.22)
labelsin <- c(33.187, 33.36, 33.51, 35.43, 35.77, 35.97, 37.58, 37.76, 38.73, 39.39, 39.71, 42.34, 45.55, 46.52, 46.68, 52.77, 76.74)
labeltax <- c('Alfisol','Alfisol','Ultisol','Ultisol','Ultisol','Alfisol', 'Ultisol','Ultisol','Alfisol','Alfisol','Alfisol','Ultisol','Alfisol','Inceptisol','Alfisol','Ultisol','Ultisol')

# Plot
par(mar=c(1,1,3,1))
plot(x.22, plot.order=new.order, label='Current_Order', color='colcode', name='hzn_desgn', y.offset=17, id.style = 'top', abbr=TRUE, abbr.cutoff=3, cex.id=0.75, cex.names=0.5, scaling.factor=0.75, axis.line.offset=-4, cex.depth.axis=0.75)
axis(1, line=-5, at=1:length(x.22), labels=round(labelsin, 0), cex.axis=0.75, font=2, col='black', col.axis='black',lwd=0.5, lty=1, tcl=0.5)
mtext('Mean Annual Precipitation (In)', side=1, line=-2, font=1, col='black', cex=1)
title('Base Saturation (By Sum of Cations) ', cex.main=1.5, font=1, line=1)
legend("top", horiz=TRUE, legend = c("0 - 35 %", "35 - 75 %"), fill=plotclr, cex=1.1, bty='n')



