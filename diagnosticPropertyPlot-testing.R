library(aqp)
library(soilDB)
library(sharpshootR)

# get pedon data
x.mu <- fetchNASIS()
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
mu <- readOGR(dsn = 'L:/NRCS/MLRAShared/CA630/FG_CA630_OFFICIAL.gdb', layer = 'ca630_a', stringsAsFactors=FALSE)

# Transform coor from UTM to geographic coordinates
mu.gcs <- spTransform(mu, CRS('+proj=longlat +datum=WGS84'))

# subset seleted map unit
s <- grep('3058', mu.gcs$MUSYM, ignore.case = TRUE)
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
head(y)

# get profiles in x whose peiid are in subset y
idx <- which(profile_id(x) %in% y$peiid)
x.mu <- x[idx, ]
plot(x.mu, label='pedon_id')
groupedProfilePlot(x.mu, groups = 'taxonname')
sort(table(x.mu$taxonname), decreasing = TRUE)

# make a new column to store generalized taxonname, and fill with NA
x.mu$generalized_taxname <- rep(NA, times=length(x.mu))

# Sort pedons with intrusive parent materials
x.mu$generalized_taxname[x.mu$taxonname %in% c('Nedsgulch', 'NEDSGULCH')] <- 'Nedsgulch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Sheepranch', 'SHEEPRANCH')] <- 'Sheepranch'
x.mu$generalized_taxname[x.mu$taxonname %in% c('Wallyhill', 'WallyHill','WALLYHILL', 'Millvilla', 'haploxeralfs with redox')] <- 'Wallyhill'
x.mu$generalized_taxname[x.mu$taxonname %in% c('FRICOT')] <- 'Fricot'
x.mu$generalized_taxname[x.mu$taxonname %in% c('MARJONMINE')] <- 'Marjonmine'

# subset only those pedons with generalized bedrock
x.mu <- subsetProfiles(x.mu, s="!is.na(generalized_taxname)")

# get depth class
sdc <- getSoilDepthClass(x.mu)
site(x.mu) <- sdc

# convert PSC to factor
x.mu$part_size_class <- factor(x.mu$part_size_class)

# work-around: new function in sharpshootR
hp <- multinominal2logical(x.mu, 'part_size_class')
site(x.mu) <- hp

# diagnostic properties to consider, no need to convert to factors
v <- c( 'argillic.horizon', 'ochric.epipedon', 'mollic.epipedon', 'very.shallow',
       'shallow', 'mod.deep', 'deep', 'very.deep', 'fine', 'fine.loamy', 'loamy.skeletal','loamy')


x <- diagnosticPropertyPlot(x.mu, v, k=7, grid.label='bedrock_kind', dend.label = 'generalized_taxname')
x <- diagnosticPropertyPlot(x.mu, v, k=5, grid.label='pedon_id', dend.label = 'generalized_taxname')

x <- diagnosticPropertyPlot2(x.mu, v, k=3, grid.label='pedon_id')
x <- diagnosticPropertyPlot2(x.mu, v, k=3, grid.label='pmorigin')


## binary + multinominal categories


# work-around: new function in sharpshootR
hp <- multinominal2logical(loafercreek, 'hillslope_pos')
site(loafercreek) <- hp

# init variable names
v <- c('lithic.contact', 'paralithic.contact', 'argillic.horizon', 
       'cambic.horizon', 'ochric.epipedon', 'mollic.epipedon', 'very.shallow',
       'shallow', 'mod.deep', 'deep', 'very.deep', levels(loafercreek$hillslope_pos))

# should work as before
x <- diagnosticPropertyPlot(loafercreek, v, k=5, grid.label='bedrock_kind', dend.label = 'taxonname')
x <- diagnosticPropertyPlot2(loafercreek, v, k=5, grid.label='bedrock_kind')

loafercreek$hillslope_pos[is.na(loafercreek$hillslope_pos)] <- 'missing'
groupedProfilePlot(loafercreek, groups = 'hillslope_pos')
