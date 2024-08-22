#--------------------------------------------------
# Appendix 5. LIVE_DEAD SEAGRASS COMPARISONS
# 
# Last updated: August 21 2024
# Written by M. Kowalewski (kowalewski@ufl.edu)
#
# files archived at https://github.com/MJKowalewski/Seagrass-Fidelity-Project
#
#--------------------------------------------------

# Initial steps ====
#         Upload libraries, define functions, set parameter values, 
#         initial data summary and basic quality control
outPDF <- F # "T" outputs PDF figures, "F" prints figures in Rstudio console 
options(scipen=100)
pdf.options(paper="special", onefile=TRUE, family='Helvetica',
            pointsize=10, encoding="ISOLatin1.enc")
outPDF <- TRUE # set to 'T' or TRUE to output pdf files with figures
pdf.options(paper="special", onefile=TRUE, family='Helvetica',
            pointsize=10, encoding="ISOLatin1.enc")
library(plyr)
library(vegan)
# to install PaleoFidelity package run the next line
# devtools::install_github('mjkowalewski/PaleoFidelity', build_vignettes = TRUE)
library(PaleoFidelity)

# define custom functions ====
# a function to provide quick summary of a dataset
data.summary <- function(x, name) {
  y <- x[,-(1:loc.vars)]
  y2 <- rbind(total.samples = nrow(y),
              smallest.sample = min(rowSums(y)),
              mean.sample.size = mean(rowSums(y)),
              largest.sample = max(rowSums(y)),
              least.diverse.sample = min(rowSums(y>0)),
              mean.sample.diversity = mean(rowSums(y>0)),
              most.diverse.sample = max(rowSums(y>0)),
              least.abundant.taxon = min(colSums(y)),
              mean.taxon.abundance = mean(colSums(y)),
              most.abundant.taxon = max(colSums(y)),
              total.taxa = sum(colSums(y>0) > 0),
              missing.taxa = sum(colSums(y>0) == 0),
              number.of.singletons = sum(colSums(y>0)==1),
              mean.number.of.occurrences = mean(colSums(y>0)),
              most.frequently.occurring.taxon = max(colSums(y>0)))
  colnames(y2) <- name
  return(round(y2,0)) # change rounding as needed
}

# datsets and arguments ====
fd <- read.csv('appendix1.csv') # upload raw data
n.min.q <- 30 # minimum number of specimens per quadrant
loc.vars <- 14 # number of environmental and locality info variables 
fid.times <- 1000 # number of iterations, fidelity analyses (time consuming, set to 100 for trial runs)
sum(colSums(fd[,-(1:loc.vars)]) == 0) # confirm that all taxa have specimens
# map-related  datasets
env <- read.csv('appendix2.csv', stringsAsFactors = T) # long and lat coordinates of sites
map1 <- read.csv('appendix3.csv', na.strings=c('',NA,'.'), header=F) # NOAA map coordinates (study area)
map2 <- read.csv('appendix4.csv', header=F, na.strings=c(NA, '.', '')) # NOAA map coordinates (Florida outline)
plgs <- which(is.na(map1[2,]))
plgs2 <- which(is.na(map2[2,]))

# FIGURE 1 study area map ====
# NOTE: coastal coordinates downloaded from:
# https://gnome.orr.noaa.gov/goods/tools/GSHHS/coast_subset
# color definitions for Figure 1
water.col <- 'gray100'
land.col <- 'gray80'
coast.col <- 'black'
site.col <- 'black'
water.col2 <- 'gray100'
land.col2 <- 'gray70'
coast.col2 <- 'black'
# Figure 1 - plot the map
if(outPDF) pdf(paste0("Fig 1 map seagrass Florida", runif(1), ".pdf"),
               height=6, width=5)
plot.lim <- t(map1[,2:5]) # full map
plot(plot.lim, type='n', xlab='', ylab='', las=1)
 polygon(t(map1[,2:5]), col=water.col, border=NA)
 for (i in 2:length(plgs)) {
  if (i < length(plgs)) {
    a1 <- plgs[i]+1
    a2 <- plgs[i+1]-1
    polygon(t(map1[,a1:a2]), col=land.col, border=NA)
  }
 }
 points(env[,c(3,2)], pch=21, cex=1, bg='white', col=site.col)
 text(-83.31, 29.66,  'Steinchatchee', cex=0.7, col='black')
 text(-82.71, 29.18,  'Waccasassa', cex=0.7, col='black')
 text(-82.84, 28.88,  'Crystal River', cex=0.7, col='black')
 text(-82.84, 28.75,  'Homosassa', cex=0.7, col='black')
 text(-82.81, 28.66,  'Chassahowitzka', cex=0.7, col='black')
 text(-82.82, 28.55,  'Weeki Wachee', cex=0.7, col='black')
 points(c(-83.7, -83.6), c(28.55, 28.55), type='l', lwd=2, lend=3)
 text(-83.65, 28.6, '10 km', cex=0.8)
 box(lwd=1.5)
 mtext(side=2,'latitude', line=3, cex=1.2)
 mtext(side=1,'longitude', line=2.5, cex=1.2)
 polygon(t(map1[,2:5]), col=NA, border='black')
tempar <- par(new=T, fig=c(0.59, 0.913, 0.67, 0.866), mar=c(0,0,0,0))
 plot.lim <- t(map2[,2:5]) # full map
 plot(plot.lim, type='n', xlab='', ylab='', las=1, axes=F)
 polygon(t(map2[,2:5]), col=water.col2, border=NA)
 for (i in 2:length(plgs2)) {
  if (i < length(plgs2)) {
    a1 <- plgs2[i]+1
    a2 <- plgs2[i+1]-1
    polygon(t(map2[,a1:a2]), col=land.col2, border='black')
  }
 }
 rect(-83.55, 28.5, -82.6, 29.75, border='black', lwd=1,
     col=adjustcolor('gray30', .2))
 text(-81.5, 27.5, 'Florida', srt=280, cex=0.7)
 text(-84.8, 27, 'Gulf of Mexico', srt=0, cex=0.7)
 polygon(t(map2[,2:5]), col=NA, border='black')
par(tempar)
if(outPDF) dev.off() 

# FIDELITY ANALYSES ====
# FRACTION LEVEL ANALYSES ====
# create a factor for grouping samples by quadrat and fraction
sam <- as.factor(paste(fd$system, fd$site_num, fd$station, fd$fraction))
# Find live-dead pairs of samples that meet sample size requirements
Frac <- NULL
for (i in levels(sam)) {
  sam.set <- which(sam == i)
  if (length(sam.set) == 2) {
    if (min(rowSums(fd[sam.set,-(1:loc.vars)])) >= n.min.q) {
       Frac <- rbind(Frac, fd[sam.set,])
    }
  }
}

paste('samples retained =', nrow(Frac)/2,
      'min sample =', min(rowSums(Frac[,-(1:loc.vars)])))
( frac_filter <- data.summary(Frac, 'Fractions Filtered: Total') )
FracL <- droplevels(Frac[which(Frac$type == 'live'),])
FracD <- droplevels(Frac[-which(Frac$type == 'live'),])
# summary of live fraction data
( frac_l_sum <- data.summary(FracL, 'Fractions: Live') )
# summary of dead fraction data
( frac_d_sum <- data.summary(FracD, 'Fractions: Dead') )
# output fraction dataset summary
# write.csv(cbind(frac_l_sum, frac_d_sum), 'fractions.csv')
# check how many levels (groups) remain present
res1 <- FidelityEst(as.matrix(FracL[,-(1:loc.vars)]), as.matrix(FracD[,-(1:loc.vars)]),
                    t.filters = 1, report=F, iter = fid.times, tfsd='wisconsin',
                    iter2 = fid.times, sim.measure = 'bray', cor.measure = 'spearman')

# QUADRAT LEVEL ANALYSES ====,
# (group by system, site, quadrat, and type)
quad <- as.factor(paste(fd$system, fd$site_num, fd$station, fd$type))
quad.list <- split(fd, quad)
quad.out.sp <- t(sapply(quad.list, function(x) colSums(x[,-(1:loc.vars)])))
quad.out.loc <- t(sapply(quad.list, function(x) x[,1:loc.vars][1,]))
qdat <- droplevels(data.frame(quad.out.loc, quad.out.sp))
# check basic sample info
( q.sum <- data.summary(qdat, 'quadrants') )
# output data summary table
# write.csv(q.sum, 'quadrat data summary combined live_dead.csv')
quad.sam <- as.factor(paste(qdat$system, qdat$site, qdat$station))
levels(quad.sam)
quad <- NULL
for (i in levels(quad.sam)) {
  sam.set <- which(quad.sam == i)
  if (length(sam.set) == 2) {
    if (min(rowSums(qdat[sam.set,-(1:loc.vars)])) >= n.min.q) {
      quad <- rbind(quad, qdat[sam.set,])
    }
  }
}
# check basic sample info
( q.sum.filt <- data.summary(quad, 'Stations filtered') )
# write.csv(q.sum.filt, 'filtered quadrant data summary combined live_dead.csv')
QuadL <- droplevels(quad[which(quad$type == 'live'),])
QuadD <- droplevels(quad[-which(quad$type == 'live'),])
# summary of live fraction data
( qdlive <- data.summary(QuadL, 'Stations: Live') )
# summary of dead fraction data
( qddead <- data.summary(QuadD, 'Stations: Dead') )
# output quadrant dataset summary
# write.csv(cbind(qdlive, qddead), 'quadrants live and dead separately.csv')
res2 <- FidelityEst(as.matrix(QuadL[,-(1:loc.vars)]), sim.measure = 'bray',
                    as.matrix(QuadD[,-(1:loc.vars)]),
                    t.filters = 1, report=F, tfsd='wisconsin', iter = fid.times,
                    iter2 = fid.times)


# SITE LEVEL ANALYSES ====
# (group by system, site, and type)
site <- as.factor(paste(fd$system, fd$site_num, fd$type))
site.list <- split(fd, site)
site.out.sp <- t(sapply(site.list, function(x) colSums(x[,-(1:loc.vars)])))
site.out.loc <- t(sapply(site.list, function(x) x[,1:loc.vars][1,]))
sitedat <- droplevels(data.frame(site.out.loc, site.out.sp))

# check basic sample info
( site_sum <- data.summary(sitedat, 'sites') )
site.sam <- as.factor(paste(sitedat$system, sitedat$site))
levels(site.sam)
siteout <- NULL
for (i in levels(site.sam)) {
  site.sam.set <- which(site.sam == i)
  if (length(site.sam.set) == 2) {
    if (min(rowSums(sitedat[site.sam.set,-(1:loc.vars)])) >= n.min.q) {
      siteout <- rbind(siteout, sitedat[site.sam.set,])
    }
  }
}
# check basic sample info after filtering
( site_sum_filt <- data.summary(siteout, 'sites filtered') )
SiteL <- droplevels(siteout[which(siteout$type == 'live'),])
SiteD <- droplevels(siteout[-which(siteout$type == 'live'),])
# summary of live fraction data
( sitelsum <- data.summary(SiteL, 'Sites: Live') )
# summary of dead fraction data
( sitedsum <- data.summary(SiteD, 'Sites: Dead') )

res3 <- FidelityEst(as.matrix(SiteL[,-(1:loc.vars)]), sim.measure = 'bray',
                    as.matrix(SiteD[,-(1:loc.vars)]),
                    t.filters = 1, report=F, tfsd='wisconsin', iter = fid.times,
                    iter2 = fid.times)

# ESTUARY LEVEL ANALYSES ====
# (group by system and type)
estu <- as.factor(paste(fd$system, fd$type))
estu.list <- split(fd, estu)
estu.out.sp <- t(sapply(estu.list, function(x) colSums(x[,-(1:loc.vars)])))
estu.out.loc <- t(sapply(estu.list, function(x) x[,1:loc.vars][1,]))
estudat <- droplevels(data.frame(estu.out.loc, estu.out.sp))
# check basic sample info
( estu_sum <- data.summary(estudat, 'estuaries') )
estu.sam <- as.factor(paste(estudat$system))
levels(estu.sam)
estuout <- NULL
for (i in levels(estu.sam)) {
  estu.sam.set <- which(estu.sam == i)
  if (length(estu.sam.set) == 2) {
    if (min(rowSums(estudat[estu.sam.set,-(1:loc.vars)])) >= n.min.q) {
      estuout <- rbind(estuout, estudat[estu.sam.set,])
    }
  }
}
# check basic sample info after filtering
( estu_sum_filt <- data.summary(estuout, 'estuaries filtered') )
EstuL <- droplevels(estuout[which(estuout$type == 'live'),])
EstuD <- droplevels(estuout[-which(estuout$type == 'live'),])
# summary of live fraction data
( estulsum <- data.summary(EstuL, 'Estuaries: Live') )
# summary of dead fraction data
( estudsum <- data.summary(EstuD, 'Estuaries: Dead') )

res4 <- FidelityEst(as.matrix(EstuL[,-(1:loc.vars)]), sim.measure = 'bray',
                    as.matrix(EstuD[,-(1:loc.vars)]),
                    t.filters = 1, tfsd='wisconsin', report=F, iter = fid.times,
                    iter2 = fid.times)

res5 <- FidelityEst(live=rbind(colSums(fd[fd$type=='live',-(1:loc.vars)])),
                    sim.measure = 'bray',
                    dead=rbind(colSums(fd[fd$type=='dead',-(1:loc.vars)])),
                    t.filters = 1, report=F, tfsd='wisconsin', iter = fid.times,
                    iter2 = fid.times)


# Tables and data summaries ====
# basic info about the data
( dsumall <- data.summary(fd, 'all data') )

# Summary of all datasets ====
( data.summary.all <- cbind(dsumall, frac_filter, frac_l_sum,
                            frac_d_sum, q.sum.filt, qdlive,
                            qddead, site_sum_filt, sitelsum,
                            sitedsum, estu_sum_filt,
                            estulsum, estudsum) )

# Fidelity Summary ====
fidstatsF <- function(x, name) {
  y <- cbind(c(length(x$x), min(x$x), max(x$x), x$x.stats, median(x$x),
               min(x$xc), max(x$xc),x$xc.stats, min(x$xc), max(x$xc), 
               x$xs.stats, x$values$min.sam,
               min(x$x.pf.dist), max(x$x.pf.dist),
               mean(x$x.pf.dist), median(x$x.pf.dist)))
  rownames(y) <- c('number of samples', 'minimum [observed]',
                   'maximum [observed]', 'mean [observed]',
                   'median [observed]', 'minimum [corr]',
                   'maximum [corr]','mean [corr]',
                   'lower CL [corr]', 'median [corr]', 'upper CL [corr]',
                   'minimum [samsdt]', 'maximum [samstd]', 'mean [samstd]',
                   'lower CL [samstd]', 'median [samstd]', 'upper CL [samstd]', 'samstd n',
                   'minimum [perfect]', 'maximum [perfect]',
                   'mean [expected perfect]', 'median [expected perfect]')
  colnames(y) <- name
 return(y)  
}


#---------------------------------------------
# Diversity Analyses standardized at n = 30) ====
dtstDIV <- list("fractions live" = FracL[,-(1:loc.vars)],
                "fractions dead" = FracD[,-(1:loc.vars)],
                "quadrants live" = QuadL[,-(1:loc.vars)],
                "quadrants dead" = QuadD[,-(1:loc.vars)],
                "sites live" = SiteL[,-(1:loc.vars)],
                "sites dead" = SiteD[,-(1:loc.vars)],
                "estuaries live" = EstuL[,-(1:loc.vars)],
                "estuaries dead" = EstuD[,-(1:loc.vars)])
pie.f <- function(x) (sum(x > 0) / (sum(x > 0) - 1)) * (1 - sum((x / sum(x)) ^ 2))
standR <- rbind(sapply(dtstDIV, function(x) mean(rarefy(x, 30))),
                sapply(dtstDIV, function(x) quantile(rarefy(x, 30))))
evenn <- rbind(sapply(dtstDIV, function(x) mean(apply(rrarefy(x, 30), 1, pie.f), na.rm=T)),
               sapply(dtstDIV, function(x) quantile(apply(rrarefy(x, 30), 1, pie.f), na.rm=T)),
               sapply(dtstDIV, nrow))

# Diversity Offsets
repdiv1 <- FidelityDiv(as.matrix(FracL[,-(1:loc.vars)]),
                       as.matrix(FracD[,-(1:loc.vars)]),
                       t.filters = 1, iter = 1000, CImean = 0.95)
repdiv2 <- FidelityDiv(as.matrix(QuadL[,-(1:loc.vars)]),
                       as.matrix(QuadD[,-(1:loc.vars)]),
                       t.filters = 1, iter = 1000, CImean = 0.95)
repdiv3 <- FidelityDiv(as.matrix(SiteL[,-(1:loc.vars)]),
                       as.matrix(SiteD[,-(1:loc.vars)]),
                       t.filters = 1, iter = 1000, CImean = 0.95)
repdiv4 <- FidelityDiv(as.matrix(EstuL[,-(1:loc.vars)]),
                       as.matrix(EstuD[,-(1:loc.vars)]),
                       t.filters = 1, iter = 1000, CImean = 0.95)

rdivlist <- list(fractions=repdiv1, quadrats=repdiv2, sites=repdiv3,
                 estuaries=repdiv4)
divoff <- t(rbind(sapply(rdivlist, function(x) x$xmean),
                sapply(rdivlist, function(x) x$p.values[1]),
                sapply(rdivlist, function(x) x$ymean),
                sapply(rdivlist, function(x) x$p.values[2]),
                sapply(rdivlist, function(x) nrow(x$x))))
colnames(divoff) <- c('S', 'S2.5', 'S97.5', 'pS', 'E', 'E2.5',
                      'E97.5', 'pE', 'n') 

# Table 1====
Table1 <- cbind('number of samples' = data.summary.all[1,c(3,6,9,12)])
rownames(Table1) <- c('fractions', 'quadrats', 'sites', 'estuaries')
Table1
write.csv(Table1, 'table1.csv')

dim(data.summary.all)
# Table 2====
Table2 <- data.summary.all[,c(1,3,4,6,7,9,10,12,13)]
Table2
write.csv(Table2, 'table2.csv')

# Table 3====
Table3 <- cbind(fidstatsF(res1, 'fractions'), fidstatsF(res2, 'quadrats'),
                fidstatsF(res3, 'sites'), fidstatsF(res4, 'estuaries'))
Table3 <- Table3[c(1,4,8,9,11,14,15,17,21,18),]
Table3
write.csv(Table3, 'table3.csv')

# Table 4====
div_summ <- cbind(t(standR), t(evenn))
colnames(div_summ) <- c('mean S', 'min S', 'q25 S', 'median S', 'q75 S',
                        'max S', 'mean E', 'min E', 'q25 E', 'median E',
                        'q75 E', 'max E', 'n')
Table4 <- div_summ[,c(13, 1,3,5,7,9,11)]
Table4
write.csv(Table4, 'table4.csv')

# Table 5====
table5 <- as.data.frame(divoff)
table5
write.csv(table5, 'Table5.csv')

# FIGURES ====
#---------------------------------------------
# Comparison of sample sizes between live & dead sites & quadrants ====
siteliveN <- rowSums(SiteL[,-(1:loc.vars)])
sitedeadN <- rowSums(SiteD[,-(1:loc.vars)])
quadliveN <- rowSums(QuadL[,-(1:loc.vars)])
quaddeadN <- rowSums(QuadD[,-(1:loc.vars)])
myxlabs <- c(10, expression(10^2), expression(10^3),
             expression(10^4), expression(10^5))
quick.sum.F <- function(x) {
  c(n=length(x), mean=mean(x), median=median(x),
    std.dev=sd(x), min=min(x), max=max(x))
}

##### FIGURE 2 live dead sample size comparison ====
if (outPDF) pdf('Figure 2 live dead sample size comparison.pdf',
                width = 5, height = 5)
tempar <- par(mfrow=c(2,1), mar=c(2,0,1,0), oma=c(3,4,0,1))
plot(siteliveN, sitedeadN, log='xy', xlim=c(10,100000),
     ylim=c(10,100000), type='n', axes=F,
     xlab='', ylab='')
axis(1, at=c(10,100,1000,10^4,10^5), labels=myxlabs)
axis(2, at=c(10,100,1000,10^4,10^5), labels=myxlabs, las=1)
box()
abline(a = 0, b = 1, lwd=2)
points(c(1,10^6), c(10,10^7), type='l')
points(c(1,10^6), c(100,10^8), type='l')
points(c(1,10^6), c(1000,10^9), type='l')
points(siteliveN, sitedeadN, pch=21, col='black',
       bg=adjustcolor('gray30', .5), cex=1.2)
mtext(side=3, line=-1, adj=0.01, 'A')
plot(quadliveN, quaddeadN, log='xy', xlim=c(10,100000),
     ylim=c(10,100000), type='n', axes=F)
axis(1, at=c(10,100,1000,10^4,10^5), labels=myxlabs)
axis(2, at=c(10,100,1000,10^4,10^5), labels=myxlabs, las=1)
box()
abline(a = 0, b = 1, lwd=2)
points(c(1,10^6), c(10,10^7), type='l')
points(c(1,10^6), c(100,10^8), type='l')
points(c(1,10^6), c(1000,10^9), type='l')
points(quadliveN, quaddeadN, pch=21, col='black',
       bg=adjustcolor('gray30', .5), cex=1.5)
mtext(side=3, line=-1, adj=0.01, 'B')
mtext(outer=T, side=2, line=2.5, 'number of dead specimens')
mtext(outer=T, side=1, line=0.5, 'number of live specimens')
par(tempar)
if(outPDF) dev.off()

#### DIVERSITY/EVENNESS OFFSET FIGURES ====
col.sym.pt <- 'darkgray'
bg.sym.pt <- 'gray'
col.cf.col <- 'gray'
col.colmean <- 'black'
my.transp <- 0.7
max.x <- max(abs(c(repdiv1$x[,2],repdiv2$x[,2],repdiv3$x[,2],repdiv4$x[,2])))
max.y <- max(abs(c(repdiv1$y[,2],repdiv2$y[,2],repdiv3$y[,2],repdiv4$y[,2])))

######### Figure 3 with Spearman Plots
histlist <- list(res4$x, apply(res4$x.pf.dist, 2, mean),
                 res3$x, apply(res3$x.pf.dist, 2, mean),
                 res2$x, apply(res2$x.pf.dist, 2, mean),
                 res1$x, apply(res1$x.pf.dist, 2, mean))
gpnames <- c("estuaries", "sites", "quadrats", "quadrat-fractions")
bar.col <- 'white'
border.col <- 'black'
bar.col.mod <- 'gray70'
border.col.mod <- 'gray40'

if(outPDF) pdf(paste0('Figure 3 spearman plots_',runif(1),'.pdf'),
               width=5, height=7)
tempar <- par(mfrow=c(4,1), mar=c(0,0,2,2), oma=c(5,5,0,0))
j <- 0
for (i in seq(1, 7, 2)) {
  j <- j + 1
  if (j == 1) my.ylim <- 5
  if (j == 2) my.ylim <- 10
  if (j == 3) my.ylim <- 20
  if (j == 4) my.ylim <- 30
  hist(histlist[[i+1]], breaks=seq(-1,1,0.05), main='', axes=F,
       ylab='', xlab='', col=adjustcolor(bar.col.mod, 0.3),
       border=NA, ylim=c(0,my.ylim)) 
  points(mean(histlist[[i+1]]), 0, pch=16, col=bar.col.mod, xpd=NA, cex=3)
  points(mean(histlist[[i]]), 0, pch=16, col=border.col, xpd=NA, cex=3)
  hist(histlist[[i]], breaks=seq(-1,1,0.05), main='', axes=F,
       ylab='', xlab='', col=adjustcolor(bar.col, 0),
       border=adjustcolor(border.col, 1),
       ylim=c(0,30), add=T) 
  axis(2, las=2, at=seq(0,30,5), labels=seq(0,30,5))
  axis(1, at=seq(-1,1,0.2), labels=seq(-1,1,0.2))
  if(j > 3) mtext(side=1, line=3, bquote("spearman" ~ italic(rho)),
                  cex=1)
  mtext(side=3, line=-2, adj=0.05, LETTERS[j], cex=1.1)
  mtext(side=3, line=-4, adj=0.05, gpnames[j], cex=0.7)
  mtext(side=3, line=-5, adj=0.05, cex=0.7,
        paste('n =', length(histlist[[i]])))
  meanrho <- round(mean(histlist[[i]]),2)
  mtext(side=3, line=-6, adj=0.05, cex=0.7,
        bquote(italic(rho)[mean]==.(meanrho)))
  
}
mtext(side=2, line=3, 'number of samples', cex=1, outer=T)
par(tempar)
if(outPDF) dev.off()

# Figure 4 fidelity for all data pooled
if (outPDF) pdf('Figure 4 fidelity_pooled.pdf', width = 5, height = 5)
tmp.par <- par(mar=c(3,8,2,8))
rep4 <- LDPlot(colSums(fd[fd$type=="live" , -(1:loc.vars)]),
               colSums(fd[fd$type=="dead" , -(1:loc.vars)]),
               tax.names = colnames(fd)[-(1:loc.vars)],
               toplimit = 30, cex.names = 0.6, report=T) 
par(tmp.par)
if (outPDF) dev.off()

#### Figure 5 - FOUR PANEL PLOT OF DIVERSITY OFFSETS ====
if(outPDF) pdf(paste0('Figure 5 diversity_offset_multipanel',runif(1),'.pdf'),
               width=4.5, height=7)
tempar <- par(mfrow=c(4,1), oma=c(4,1,0.5,0.5), mar=c(0,5,0,0))
AlphaPlot(repdiv1, transp=my.transp, xlab='', col=col.sym.pt, cf.col=col.cf.col,
          bgpt=bg.sym.pt, colmean = col.colmean, axes=F, xmax=max.x, ymax=max.y)
axis(2, las=1)
mtext(side=3, line=-1.5, adj=0.05, 'A. Fractions')
AlphaPlot(repdiv2, transp=my.transp, xlab='', col=col.sym.pt, cf.col=col.cf.col,
          bgpt=bg.sym.pt, colmean = col.colmean, axes=F, xmax=max.x, ymax=max.y)
axis(2, las=1)
mtext(side=3, line=-1.5, adj=0.05, 'B. Quadrats')
AlphaPlot(repdiv3, transp=my.transp, xlab='', col=col.sym.pt, cf.col=col.cf.col,
          bgpt=bg.sym.pt, colmean = col.colmean, axes=F, xmax=max.x, ymax=max.y)
axis(2, las=1)
mtext(side=3, line=-1.5, adj=0.05, 'C. Sites')
AlphaPlot(repdiv4, transp=my.transp, xlab='', col=col.sym.pt, cf.col=col.cf.col,
          bgpt=bg.sym.pt, colmean = col.colmean, axes=F, xmax=max.x, ymax=max.y)
mtext(side=3, line=-1.5, adj=0.05, 'D. Estuaries')
mtext(side=1, line=2.5, bquote(Delta[S]))
axis(2, las=1)
axis(1)
par(tempar)
if(outPDF) dev.off()

#### Figure 6 - SUMMARY OF DIVERSITY OFFSETS ====
if(outPDF) pdf(paste0('Figure 6 diversity_offset_',runif(1),'.pdf'),
               width=7, height=4.5)
tempar <- par(mar=c(5,5,0.5,5))
plot(0,0, xlim=c(0,10), ylim=c(-0.05,0.25), type='n', axes=F,
     xlab='', ylab='')
 axis(1, at=c(1:4,6:9), labels=rep(rownames(divoff),2),
      las=2)
 axis(2, las = 1, col.axis = 'gray20')
 axis(4, las = 1, col.axis = 'gray60', col.ticks = 'gray60')
 abline(h=0, col='black', lty=3, lwd=0.5)
 box()
for (i in 1:nrow(divoff)) {
    points(c(i,i), c(divoff[i,2], divoff[i,3]), type='l', col='gray20')
    points(i, divoff[i,1], pch=16, col='gray20')
}
for (i in 1:nrow(divoff)) {
  points(c(i+5,i+5), c(divoff[i,6], divoff[i,7]), type='l', col='gray60')
  points(i+5, divoff[i,5], pch=16, col='gray60')
 }
mtext(side=2, line=3.2, bquote(Delta[S]), cex=1.3) 
mtext(side=4, line=3.2, bquote(Delta[PIE]), col='gray60', cex=1.3) 
par(tempar)
if(outPDF) dev.off()

#----------------------------------------------
# PART 5: Beta Gradient Analysis (units: quadrants) ====
# in terms of environmental distance 

out.comp <- NULL
for (i in 1:nrow(QuadL)) {
  for (j in 2:nrow(QuadL)) {
    if (j > i) {
      samesys <- sum(as.factor(QuadL$system[i]) == as.factor(QuadL$system[j]))
      samesite <- sum(as.factor(QuadL$site[i]) == as.factor(QuadL$site[j]))
      out.comp <- rbind(out.comp, cbind(samesys, samesite))
  }
 }
}
distL <- 1 - vegdist(QuadL[,-(1:loc.vars)])
distD <- 1 - vegdist(QuadD[,-(1:loc.vars)])
distF <- cbind(apply(out.comp[,1:2], 1, sum), distL, distD)
n.gp <- table(distF[,1])

# randomized means
r.mean.d <- NULL
for (i in 1:10000) {
rdistd <- sample(distF[,3])
r.gp <- c(rep(0, n.gp[1]), rep(1, n.gp[2]), rep(2,n.gp[3]))
r.mean.d <- rbind(r.mean.d, tapply(rdistd, r.gp, mean))
}

r.mean.l <- NULL
for (i in 1:10000) {
  rdistl <- sample(distF[,2])
  r.gp <- c(rep(0, n.gp[1]), rep(1, n.gp[2]), rep(2,n.gp[3]))
  r.mean.l <- rbind(r.mean.l, tapply(rdistl, r.gp, mean))
}

rmd.sim <- apply(r.mean.d, 2, function(x) c(mean(x), quantile(x, prob=0.025),
                               quantile(x, prob=0.975)))
rml.sim <- apply(r.mean.l, 2, function(x) c(mean(x), quantile(x, prob=0.025),
                                 quantile(x, prob=0.975)))

if (outPDF) pdf('Figure 7 pairwise comparisons.pdf')
tempar <- par(mfrow=c(3,2), mar=c(2,2,0,0), oma=c(3,3,0.5, 0.5))
k <- 0
for (i in 3:1) {
  for (j in 2:3) {
    k <- k + 1
    hist(distF[distF[,1]==i-1,j], breaks=seq(0,1,0.025), main='',
         xlab='', ylab='', axes=F)
    mtext(side=3, line=-1, LETTERS[k], adj=0.8)
    meanBC <- mean(distF[distF[,1]==i-1,j])
    points(meanBC, 0, pch=16, cex=2, col='black', xpd=NA)
    if (j == 2) points(rml.sim[i,], rep(0, 3), col='red', lend=3, lwd=4, type='l')
    if (j == 3) points(rmd.sim[i,], rep(0, 3), col='red', lend=3, lwd=4, type='l')
    if (j == 2) mtext(side=2, line=3, cex=0.8, 'number of comparisons')
    axis(2, las=1)
    if (i == 1) axis(1, line=1)
    if (i == 1) mtext(side=1, line=3.5, 'Bray-Curtis similarity')
}
}
par(tempar)
if (outPDF) dev.off()

### Figure 8 NMDS ====
quad.data <- rbind(QuadL[,-(1:loc.vars)], QuadD[,-(1:loc.vars)])
quad.d <- quad.data[,-which(colSums(quad.data > 0) < 2)]
res5 <- metaMDS(wisconsin(quad.d), autotransform=F, k=3, try = 50)
qtype <- as.factor(c(unlist(QuadL$type), unlist(QuadD$type)))
qsite <- as.factor(c(unlist(QuadL$env.site), unlist(QuadD$env.site)))
qsite.col <- c('yellow3', 'yellow3', 'yellow3', 'blue1', 'blue1',
               'blue1', 'blue1', 'blue1', 'green4', 'green4', 'green4',
               'red1', 'red1', 'red1', 'red1', 'black', 'black', 'black',
               'black', 'darkgray', 'darkgray')
qsite.text <- c(unlist(QuadL$site_num),
                unlist(QuadL$site_num))
pchd <- 16
pchl <- 21 
cex.sym <- 0.9
dhf <- nrow(res5$points)/2
dhf2 <- dhf + 1
dhf3 <- 2*dhf
########## FIGURE NMDS ====
if (outPDF) pdf(paste0('Figure 8 NMDS_quadrats',runif(1),'.pdf'), width = 5, height = 7)
tempar <- par(mfrow=c(3,1), mar=c(0,0,1,0), oma=c(5,5,0,1))
plot(res5$points[,1], res5$points[,2], pch=c(pchd,pchl)[qtype],
     col=qsite.col[qsite], cex=cex.sym, ylim=c(-1.5,1), xlim=c(-1,1.5),
     axes=F)
  box()
  axis(2, las=1)
  mtext(side=3, line=-1.5, adj=0.95, 'A')
  text(rep(1.45,6), seq(-0.25,-1.5,-0.25), pos=4, levels(as.factor(fd$system)),
       cex=0.7, col=unique(qsite.col))
  points(rep(1.45,6), seq(-0.25,-1.5,-0.25), pch=16, cex=cex.sym, col=unique(qsite.col))
  text(1.2, -1, 'live', pos=4, cex=cex.sym)
  text(1.2, -1.3, 'dead', pos=4, cex=cex.sym)
  points(1.2, -1, pch=21, cex=cex.sym)
  points(1.2, -1.3, pch=16, cex=cex.sym)
  text(-0.85, 0.9, paste('stress =', round(res5$stress, 3)), cex=cex.sym)
  mtext(side=2, line=3, 'NMDS 2')
plot(res5$points[1:dhf,1], res5$points[1:dhf,2], pch=pchl,
     col=qsite.col[qsite], ylim=c(-1.5,1), xlim=c(-1,1.5),
     axes=F, type='n')
  text(res5$points[1:dhf,1], res5$points[1:dhf,2],
       col=qsite.col[qsite], qsite.text, cex=cex.sym)
  box()
  axis(2, las=1)
  #axis(1, at=seq(-0.5,1,0.5), labels=seq(-0.5,1,0.5))
  mtext(side=3, line=-1.5, adj=0.95, 'B')
  mtext(side=2, line=3, 'NMDS 2')
plot(res5$points[dhf2:dhf3,1], res5$points[dhf2:dhf3,2], pch=pchd,
     col=qsite.col[qsite], ylim=c(-1.5,1), xlim=c(-1,1.5),
     axes=F, type='n')
  text(res5$points[dhf2:dhf3,1], res5$points[dhf2:dhf3,2],
     col=qsite.col[qsite], qsite.text, cex=cex.sym)
  box()
  axis(2, las=1)
  axis(1, at=seq(-1,1.5,0.2), labels=seq(-1,1.5,0.2))
  mtext(side=3, line=-1.5, adj=0.95, 'C')
  mtext(side=1, line=3, 'NMDS 1')
  mtext(side=2, line=3, 'NMDS 2')
par(tempar)
if(outPDF) dev.off()
