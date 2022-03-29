# setup selected set in NASIS

# load required packages
library(lattice)
library(soilDB)
library(sharpshootR)

# setup plotting style for later
tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkRed', 'DarkGreen'), lwd=3))

# get pedon data
x <- fetchNASIS()

# make a new column to store generalized taxonname, and fill with NA
x$generalized_taxname <- rep(NA, times=length(x))

# generalize taxonnames to group soils with slightly different names. H
x$generalized_taxname[x$taxonname %in% c('Merrimac', 'merrimac')] <- 'Merrimac'
x$generalized_taxname[x$taxonname %in% c('Enfield', 'enfield')] <- 'Enfield'

# subset only those pedons with generalized taxname
x.sub <- subsetProfiles(x, s="!is.na(generalized_taxname)")
plot(x.sub)

# tabulate number of pedons / generalized taxonname
(taxname.tab <- table(x.sub$generalized_taxname))

# custom 'slab' function, returning mean +/- 1SD
mean.and.sd <- function(values) {
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  upper <- m + s
  lower <- m - s
  res <- c(mean=m, lower=lower, upper=upper)
  return(res)
}

# aggregate several variables at once, within 'group'
a <- slab(x.sub, generalized_taxname ~ clay + sand + phfield, slab.fun=mean.and.sd)

# append the number of pedons to the taxonname label
a$generalized_taxname <- factor(a$generalized_taxname, levels=names(taxname.tab), labels=paste(names(taxname.tab), ' (', taxname.tab, ')', sep=''))


tps <- list(superpose.line=list(col=c('RoyalBlue', 'DarkGreen'), lwd=3))

# make plot with overlapping data 
xyplot(top ~ mean | variable, groups=generalized_taxname, data=a, ylab=list(label='Depth', cex=1.5),
       xlab=list(label='Mean bounded by +/- 1 Standard Deviation', cex=1.5),
       lower=a$lower, upper=a$upper, ylim=c(100,-5), 
       panel=panel.depth_function, alpha=0.5, sync.colors=TRUE,
       prepanel=prepanel.depth_function,
       cf=a$contributing_fraction,
       strip=strip.custom(bg=grey(0.85)), par.strip.text=list(cex=1.5),
       layout=c(3,1), scales=list(x=list(alternating=1, relation='free', cex=1), y=list(alternating=2, cex=1)), par.settings=tps, auto.key=list(columns=2, lines=TRUE, points=FALSE, cex=1.5))
