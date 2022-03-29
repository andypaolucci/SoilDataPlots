# setup selected set in NASIS

# load required packages
library(lattice)
library(soilDB)
library(sharpshootR)
library(aqp)
library(rgdal)
library(plyr)
library(sp)
library(cluster)
library(ape)
library(RColorBrewer)
library(latticeExtra)
library(plotrix)
library(vegan)

# setup plotting style for later
tps <- list(superpose.symbol=list(pch=16, col=c('Blue', 'Red', 'Green', 'Purple', 'Yellow', 'Grey', 'Orange', 'Brown', 'Black'), lwd=2, cex=1.5))

# get pedon data
x <- fetchNASIS()
str(x)

# keep only those pedons with real coordinates
good.idx <- which(!is.na(x$x_std))
x <- x[good.idx, ]
head(x)

#initalize spatial coordinates
coordinates(x) <- ~ x_std + y_std

# get/set spatial reference system
proj4string(x) <- '+proj=longlat +datum=WGS84'

# extract spatial data + site level attributes
x.sub <- as(x, 'SpatialPointsDataFrame')

# Get map unit data
mu <- readOGR(dsn = 'L:/CA630/FG_CA630_OFFICIAL.gdb', layer = 'ca630_a', stringsAsFactors=FALSE)

# Transform coor from UTM to geographic coordinates
mu.gcs <- spTransform(mu, CRS('+proj=longlat +datum=WGS84'))

# subset seleted map unit
s <- grep('8174', mu.gcs$MUSYM, ignore.case = TRUE)
s.sub <- mu.gcs[s , ]

# Checks
# Make coordinate systems the same
proj4string(x.sub) 
proj4string(s.sub)
# Graph
plot(s.sub)
points(x.sub, col=2, cex=0.5)

# check class
class(s.sub)
class(x.sub)

#extract pedons within map unit  
x.sub <- over(s.sub, x.sub)

#remove NAs
x.sp <- which(!is.na(x.sub$pedon_id))
y <- x.sub[x.sp, ]

# get profiles in x whose peiid are in subset y
idx <- which(profile_id(x) %in% y$peiid)
x.mu <- x[idx, ]
groupedProfilePlot(x.mu, groups = 'taxonname', label='pedon_id', color='clay')
sort(table(x.mu$taxonname), decreasing = TRUE)

# make a new column to store generalized taxonname, and fill with NA
x.mu$generalized_taxname <- rep(NA, times=length(x.mu))

# Sort pedons with intrusive parent materials
x.mu$generalized_taxname[x.mu$taxonname %in% c('Nedsgulch','Jocal')] <- 'Nedsgulch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Sheepranch')] <- 'Sheepranch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Wallyhill')] <- 'Wallyhill'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Fricot')] <- 'Fricot'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Loamy-skeletal, mixed, super, acid, thermic Lithic Xerorthen')] <- 'Xerorthents'


# subset only those pedons with generalized taxname
x.mu <- subsetProfiles(x.mu, s="!is.na(generalized_taxname)")

#Get soil depth class
sdc <- getSoilDepthClass(x.mu)
site(x.mu) <- sdc

# check
sort(table(x.mu$generalized_taxname), decreasing = TRUE)

par(mar=c(1,1,7,1))
plotSPC(x.mu, plot.order=order(x.mu$depth), color='clay', col.label='Clay (%)', label='pedon_id', alt.label = 'generalized_taxname', alt.label.col='black', cex.names = 0.3)  
groupedProfilePlot(x.mu, groups = 'generalized_taxname', break.style='line', id.style='side', label='pedon_id', group.name.offset =-10, group.line.col=3, group.line.lwd=2, group.line.lty=5, width= 0.15, y.offset=1, cex.names=0.4, group.name.cex=0.4)


#Profile Dissimilarity 
#set up colors
cols <- brewer.pal(7, 'Set1')

# convert SPC -> DF
x.mu.df <- as(x.mu, 'data.frame')

# normalize taxonnames, and convert to factor
x.mu.df$generalized_taxname <- factor(tolower(x.mu.df$generalized_taxname))

# inspect variables used to determine dissimilarity
#3D plot
cloud(clay ~ phfield + total_frags_pct, groups=generalized_taxname, data=x.mu.df, 
      auto.key=list(space='right', points=TRUE, lines=FALSE), 
      par.settings=list(superpose.symbol=list(pch=16, col=cols, cex=1.5)))

#2D plot
xyplot(clay ~ total_frags_pct, groups=generalized_taxname, data=x.mu.df, 
      auto.key=list(space='right', points=TRUE, lines=FALSE), 
      par.settings=list(superpose.symbol=list(pch=16, col=cols, cex=1.5)))

# particle size control section data
# get PSCS depths
f.psc.prop <- function(x, prop) {
  # these are accessed from @site
  sv <- c(x$psctopdepth, x$pscbotdepth)
  
  # test for missing PCS data
  if(any(is.na(sv)))
    return(NA) 
  
  # check to make sure that the lower PSC boundary is shallower than the depth
  if(sv[2] > max(x))
    return(NA)
  
  # create formula from named property
  fm <- as.formula(paste0(' ~ ', prop))
  
  # return just the (weighted) mean, accessed from @horizons
  s <- slab(x, fm, slab.structure=sv, slab.fun=mean, na.rm=TRUE)$value
  return(s)
}

# calculating PCS weighted means
x.mu$pcs.clay <- profileApply(x.mu, f.psc.prop, prop='clay')
x.mu$pcs.frags <- profileApply(x.mu, f.psc.prop, prop='total_frags_pct')
x.mu$soil.depth <- profileApply(x.mu, estimateSoilDepth)

# normalize taxonname
x.mu$generalized_taxname <- tolower(x.mu$generalized_taxname)

#plot PCS means
xyplot(pcs.clay ~ pcs.frags, groups=generalized_taxname, data=site(x.mu), 
       auto.key=list(space='right', points=TRUE, lines=FALSE), 
       par.settings=list(superpose.symbol=list(pch=16, col=cols, cex=1.5)))

# compute between-profile dissimilarity, using pedon-level data
# use functions from cluster package

# distance matrix, requires FULL dataset = no NA
# simpler to work with just site data
s <- site(x.mu)

# check univariate group-wise patterns
bwplot(generalized_taxname ~ pcs.clay, data=s)
bwplot(generalized_taxname ~ pcs.frags, data=s)
bwplot(generalized_taxname ~ soil.depth, data=s)

d <- daisy(s[, c('pcs.clay', 'pcs.frags', 'soil.depth')], stand = TRUE)

# group via divisive hierarchical clustering
d.diana <- diana(d)

# convert to hclust object and copy labels
h <- as.hclust(d.diana)
h$labels <- paste0(s$generalized_taxname, ' ', s$pedon_id)

# convert classes, for better plotting
d.phylo <- as.phylo(h)

# plot dendrogram in next panel
plot(d.phylo, label.offset=0.1, cex=0.75, no.margin=TRUE, panel=2)
plot(d.phylo)

xyplot(soil.depth ~ pcs.frags, groups=generalized_taxname, data=site(x.mu), 
       auto.key=list(space='right', points=TRUE, lines=FALSE), 
       par.settings=list(superpose.symbol=list(pch=16, col=cols, cex=1.5)))


cloud(soil.depth ~ pcs.frags + pcs.clay, groups=generalized_taxname, data=site(x.mu), 
      auto.key=list(space='right', points=TRUE, lines=FALSE), 
      par.settings=tps)



