
fetchNASISDiagnostics <- function() {
  if(!requireNamespace('RODBC')) stop('please install the `RODBC` package', call.=FALSE)
  q <- "SELECT peiidref as peiid, dfk.ChoiceName as diag_kind, featdept, featdepb\nFROM pediagfeatures_View_1 \n\tLEFT OUTER JOIN (SELECT * FROM MetadataDomainDetail WHERE DomainID = 147) AS dfk ON pediagfeatures_View_1.featkind = dfk.ChoiceValue\n\tORDER BY pediagfeatures_View_1.peiidref, pediagfeatures_View_1.featdept;"
  channel <- RODBC::odbcDriverConnect(connection="DSN=nasis_local;UID=NasisSqlRO;PWD=nasisRe@d0n1y")
  d <- RODBC::sqlQuery(channel, q, stringsAsFactors=FALSE)
  RODBC::odbcClose(channel)
  return(d)  
}

f <- fetchNASIS()
foo <- fetchNASISDiagnostics()
#calc thickness
foo$thickness <- foo$featdepb-foo$featdept
#filter to desired diagnostic
foo2 <- foo[which(foo$diag_kind == "mollic epipedon"),]
foo3 <- data.frame(peiid=foo2$peiid,mollic_thickness=foo2$thickness)
site(f) <- merge(site(f), foo3, by="peiid")


#check
head(f)

#rename x so we can save to shapefile

x <- f 
# setup plotting style for later
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen'), lwd=2))

#Load raster layers:
lastfrost <- raster('E:/OR644_Geodata/last_spring_frost_mean_800m.tif')
temperature <- raster('E:/OR644_Geodata/final_MAAT_800m.tif')
precip <- raster('E:/OR644_Geodata/final_MAP_mm_800m.tif')
elevation <- raster('C:/Paolucci/OR644/dem_10m_clip')
slope <- raster('C:/Paolucci/OR644/slope_10m')
frostfree <- raster('E:/OR644_Geodata/ffd_mean_800m.tif')


# keep only those pedons with real coordinates
good.idx <- which(!is.na(x$x_std))
x <- x[good.idx, ]

#iniitalize spatial coordinates
coordinates(x) <- ~ x_std + y_std

# get/set spatial reference system
proj4string(x) <- '+proj=longlat +datum=WGS84'

# extract spatial data + site level attributes
x.sp <- as(x, 'SpatialPointsDataFrame')

# sample raster using series data points
x.frost <- extract(lastfrost, x.sp)
x.temp <- extract(temperature, x.sp)
x.precip <- extract(precip, x.sp)
x.elev <- extract(elevation, x.sp)
x.slope <- extract(slope, x.sp)
x.frostfree <- extract(frostfree, x.sp)


# dave to dataframe
d <- data.frame(site_id= x.sp$site_id, series = x.sp$taxonname, mollic=x.sp$mollic.epipedon, mollic_thickness= x.sp$mollic_thickness, last_frost = x.frost, temp = x.temp, precip = x.precip, elev = x.elev, slope = x.slope, frostfree = x.frostfree)

bwplot(precip ~ mollic2, data = d)

table(d$mollic)


xyplot(mollic_thickness ~ precip, data=d, par.settings =list(superpose.symbol = list(pch=17, cex=2, col=c( "purple"))

d$mollic2 <- as.character(d$mollic)                                                            
                            
head(d)                                                            


#Add data to spatial points data frame 
x.sp$frost <- x.frost 
x.sp$temp <- x.temp
x.sp$precip <- x.precip
x.sp$elev <- x.elev 
x.sp$slope <- x.slope 
x.sp$frostfree <- x.frostfree

head(x.sp)

#Remove Unused Columns 
x.sp2 <- x.sp[, c('site_id', 'x', 'y', 'utmnorthing', 'utmeasting',  'class_type', 
                    'taxonname', 'elev_field', 'mollic.epipedon', 'mollic_thickness')]

#Write to shapefile
writeOGR(x.sp2, dsn='C:/Paolucci/OR644/Mollic', layer='OR644_Mollics_8_10_2017AJP', driver='ESRI Shapefile', overwrite_layer=TRUE)
