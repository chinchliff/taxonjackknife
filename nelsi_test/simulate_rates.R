library(NELSI)
library(phangorn)

plot <- function(s, x, y) {
	ntrees <- length(s)
	ncols <- round(sqrt(ntrees * x / y))
	nrows <- round(sqrt(ntrees * y / x))
	nrows <- ifelse(nrows*ncols >= ntrees, nrows, nrows+1)
	dev.new(width=x, height=y)
	plot.new()
	layout(matrix(c(1:ntrees, rep(0, (nrows * ncols - ntrees))), nrows, ncols, byrow=T))
	layout.show(ntrees)
	sapply(1:ntrees, function(i) {
		par(mar=c(1,1,1,1))
		plot.phylo(s[[i]][[1]],show.tip.label=F,edge.width=0.2) })	
}

simulate.rates <- function(i) {

	simtype <- d[i,1]
	ntrees <- as.integer(d[i,2])
	v <- type.convert(d[i,3])

	t <- read.tree(file.path(basedir, 'source_trees', paste(simtype,'.tre', sep='')))
	s <- lapply(1:ntrees, function(i) {
		simulate.autocor.kishino(t, params = list(initial.rate = 0.2, v = v)) })
	
	plot(s, 8, 10)

	outdir <- file.path(basedir, 'scaled_trees', paste(simtype, ntrees, sep="_"))
	dir.create(outdir)
	
	sapply(1:ntrees, function(i) {
		write.tree(s[[i]][[1]], file=file.path(outdir, paste(i, '.tre', sep=''))) })
	
}

basedir <- '/Users/cody/Dropbox/Projects_current/Jackknife_test/nelsi_test'

r = 0.5 # heritability (autocorrelation coefficient)
d <- rbind(
#	c('balanced_10000', 1, r),
	c('balanced_1000', 5, r),
#	c('balanced_500', 10, r),
#	c('balanced_100', 50, r),
#	c('balanced_50', 100, r),
	c('balanced_100', 100, r),
	c('pectinate_1000', 5, r / log(1000, 2)),
#	c('pectinate_500', 10, r / log(500, 2)),
#	c('pectinate_100', 50, r / log(100, 2)),
#	c('pectinate_50', 100, r / log(50, 2)),
	c('balanced_100', 100, r))

simulate.rates(8)