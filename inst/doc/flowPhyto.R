### R code from vignette source 'flowPhyto.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("flowPhyto")


###################################################
### code chunk number 2: loadlib
###################################################
library(flowPhyto)


###################################################
### code chunk number 3: read
###################################################
CHANNEL.CLMNS


###################################################
### code chunk number 4: read
###################################################
evt.file.path <- system.file("extdata","seaflow_cruise","2011_001", "2.evt",
				package="flowPhyto")
evt <- readSeaflow(evt.file.path)


###################################################
### code chunk number 5: filter
###################################################
opp <- filter(evt, notch=1.1)


###################################################
### code chunk number 6: classify1
###################################################
opp.path <- system.file("extdata","seaflow_cruise","2011_001", "2.evt.opp", 
				package="flowPhyto")
pop.def.path <- system.file("extdata","seaflow_cruise","pop.def.tab", 
				package="flowPhyto")
opp <- readSeaflow(opp.path)
def <- readPopDef(pop.def.path)
def


###################################################
### code chunk number 7: classify1
###################################################
pop <- classify(x=opp, pop.def= def, func=2 )

table(pop$pop)


###################################################
### code chunk number 8: plotCytogram
###################################################
plotCytogram(pop, "fsc_small","chl_small", pop.def= def, add.legend=TRUE, cex=1)



###################################################
### code chunk number 9: fig1
###################################################
plotCytogram(pop, "fsc_small","chl_small", pop.def= def, add.legend=TRUE, cex=1)



###################################################
### code chunk number 10: consensus
###################################################
vct.paths <- sapply(c(1,439,440), function(i) 
		system.file("extdata","seaflow_cruise","2011_001", 
			paste("1.evt.opp.",i,'-class.vct',sep=''), 
				package="flowPhyto"))
mat <- do.call(cbind,lapply(vct.paths, read.delim))
consen.df <- consensus(mtrx=mat)
table(consen.df$pop)
aggregate(consen.df$support,list(consen.df$pop), mean)


###################################################
### code chunk number 11: census
###################################################
census(v=pop$pop, pop.def=def)


###################################################
### code chunk number 12: summarize
###################################################
filter.df <- readSeaflow(opp.path, add.yearday.file=TRUE)
classed <- cbind.data.frame(filter.df, consen.df)
names(opp.path) <- getFileNumber(opp.path)
class.jn <- joinSDS(classed, opp.path)
nrow.opp <- sapply(opp.path, function(p) readSeaflow(              p , count.only=TRUE))
nrow.evt <- sapply(opp.path, function(p) readSeaflow(sub('.opp','',p), count.only=TRUE))
class.jn$opp <- rep(nrow.opp, times=nrow.opp)
class.jn$evt <- rep(nrow.evt, times=nrow.opp)

summarize(class.jn,  opp.paths.str=opp.path)


###################################################
### code chunk number 13: plotStatMap
###################################################
stat.tab <- system.file("extdata","seaflow_cruise","stats.tab", 
					package="flowPhyto")
stats <- read.delim(stat.tab)
plotStatMap(df=stats, pop='ultra', z.param='conc', margin=0.2, zlab=expression(paste('Cell concentration, ',10^6 * cells/L)), 
main="Cell concentration of Ultra-plankton population")
mtext(line=1, side=4, "cell concentration 10^6 cells / L")



###################################################
### code chunk number 14: fig2
###################################################
stat.tab <- system.file("extdata","seaflow_cruise","stats.tab", 
					package="flowPhyto")
stats <- read.delim(stat.tab)
plotStatMap(df=stats, pop='ultra', z.param='conc', margin=0.2, zlab=expression(paste('Cell concentration, ',10^6 * cells/L)), 
main="Cell concentration of Ultra-plankton population")
mtext(line=1, side=4, "cell concentration 10^6 cells / L")



###################################################
### code chunk number 15: validateSDS
###################################################
path <- system.file("extdata","seaflow_cruise",package="flowPhyto")
sds <- combineSdsFiles(path)
plot(sds$LON, sds$LAT)


###################################################
### code chunk number 16: validatePopDef
###################################################
validatePopDef(readPopDef(pop.def.path))


###################################################
### code chunk number 17: pipeline
###################################################
repository.dir <- '.' 
output.path <- paste(repository.dir,'/','seaflow_cruise',sep='')
seaflow.path <- system.file("extdata", 'seaflow_cruise', package="flowPhyto")
file.copy(from=seaflow.path, to=repository.dir, recursive=TRUE)
pipeline(repo= repository.dir, cruise.name='seaflow_cruise', steps=4, parallel=FALSE, submit.cmd='qsub -l walltime=00:20:00') 
unlink(output.path, recursive=TRUE)


