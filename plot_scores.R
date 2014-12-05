plot.ica.vs.node.depth <- function(d, write.x.label=FALSE, write.y.axis=FALSE) {	

	x1 <- d[,1]
	x2 <- d[,2]
	df <- data.frame(x1,x2)
	
	# Use densCols() output to get density at each point
	x <- densCols(x1, x2, colramp=colorRampPalette(c("black", "white")))
	df$dens <- col2rgb(x)[1,] + 1L
	
	# Map densities to colors
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$dens]

#	plot(d[,2] ~ d[,1], pch=20, cex=0.5, col=rgb(0,0,0,0.15), ylim=c(-1,1), axes=FALSE)
	plot(x2~x1, data=df[order(df$dens),], pch=20, cex=0.5, col=col, ylim=c(-1,1), axes=FALSE)

	if (write.y.axis) {
		axis(2, cex.axis=0.75, hadj=0.5, padj=1.5, tcl=-0.2)
		mtext(side=2, cex=0.75, text="ica", line=1)
	}
	axis(1,cex.axis=0.75, hadj=0.5, padj=-1.75, tcl=-0.2)
	if (write.x.label) {
		mtext(side=1, cex=0.75, text="node depth", line=1)
	}
	box()
	y <- lm(d[,2] ~ d[,1])
#	abline(y,col=col)
}

plot.j.ica.vs.b.ica.density <- function(a, write.x.label=FALSE) {
	## 3. Create a color density plot to compare jackknife ica to bootstrap ica for nodes not in the original topology
	## from http://stackoverflow.com/questions/17093935/r-scatter-plot-symbol-color-represents-number-of-overlapping-points

	x1 <- a[,"b_ica"]
	x2 <- a[,"j_ica"]
	df <- data.frame(x1,x2)
	
	# Use densCols() output to get density at each point
	x <- densCols(a[,"b_ica"],a[,"j_ica"], colramp=colorRampPalette(c("black", "white")))
	df$dens <- col2rgb(x)[1,] + 1L
	
	# Map densities to colors
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$dens]
	
	# Plot it, reordering rows so that densest points are plotted on top
	plot(x2~x1, data=df[order(df$dens),], pch=20, cex=0.5, col=col, ylim=c(-1,1), xlim=c(-1,1),axes=FALSE)
#	mtext(side=2, cex=0.75, text="quartet jackknife ica", line=1)
#	axis(2,cex.axis=0.75,hadj=0.5,padj=1.5,tcl=-0.2)
	axis(1, cex.axis=0.75, hadj=0.5, padj=-1.75, tcl=-0.2)
	if (write.x.label) {
		mtext(side=1, cex=0.75, text="bootstrap ica", line=1)
	}
	box()
	abline(0,1,lty=3)
}

c1="#0099ff"
c2="#ff6600"

w <- 8
h <- 10
dev.new(width=w,height=h)
close.screen(all=TRUE)
plot.new()

# margins
om <- 0.005 # 0.01 		# outer
lm <- 0.1 				# left
rm <- 1 - om		 	# right
tm <- 1 - om * w / h	# top
bm <- 0.03 * w / h		# bottom

bh <- 0.03				# border - horizontal
bv <- bh - om 			# border - vertical	-- use scaled (larger) border for vertical to allow for additional axes
t <- 2*bh
sh <- (rm-lm-3*bh) / 4
sv <- sh * w / h

top <- c(
	1+bm-t-0.5*bv-0*sv-om,
	1+bm-t-1.5*bv-1*sv-om,
	1+bm-t-2.5*bv-2*sv-om,
	1+bm-t-3.5*bv-3*sv-om,
	1+bm-t-4.5*bv-4*sv-om)

bottom <- c(
	tm+bm-t-0.5*bv-1*sv,
	tm+bm-t-1.5*bv-2*sv,
	tm+bm-t-2.5*bv-3*sv,
	tm+bm-t-3.5*bv-4*sv,
	tm+bm-t-4.5*bv-5*sv)

left <- c(
	lm+0*sh+0*bh,
	lm+1*sh+1*bh,
	lm+2*sh+2*bh,
	rm-sh)

right <- c(
	lm+1*sh+0*bh,
	lm+2*sh+1*bh,
	lm+3*sh+2*bh,
	rm)

#	left			right			bottom			top
m <- rbind(

# column labels
c(	left[1],		right[1],		tm+bm-t,			tm),
c(	left[2],		right[2],		tm+bm-t,			tm),
c(	left[3],		right[3],		tm+bm-t,			tm),
c(	left[4],		right[4],		tm+bm-t,			tm),

# row labels
c(	om,				lm,				bottom[1],		top[1]),
c(	om,				lm,				bottom[2],		top[2]),
c(	om,				lm,				bottom[3],		top[3]),
c(	om,				lm,				bottom[4],		top[4]),
c(	om,				lm,				bottom[5],		top[5]),

# plots
c(	left[1],		right[1],		bottom[1],		top[1]),
c(	left[2],		right[2],		bottom[1],		top[1]),
c(	left[3],		right[3],		bottom[1],		top[1]),
c(	left[4],		right[4],		bottom[1],		top[1]),
c(	left[1],		right[1],		bottom[2],		top[2]),
c(	left[2],		right[2],		bottom[2],		top[2]),
c(	left[3],		right[3],		bottom[2],		top[2]),
c(	left[4],		right[4],		bottom[2],		top[2]),
c(	left[1],		right[1],		bottom[3],		top[3]),
c(	left[2],		right[2],		bottom[3],		top[3]),
c(	left[3],		right[3],		bottom[3],		top[3]),
c(	left[4],		right[4],		bottom[3],		top[3]),
c(	left[1],		right[1],		bottom[4],		top[4]),
c(	left[2],		right[2],		bottom[4],		top[4]),
c(	left[3],		right[3],		bottom[4],		top[4]),
c(	left[4],		right[4],		bottom[4],		top[4]),
c(	left[1],		right[1],		bottom[5],		top[5]),
c(	left[2],		right[2],		bottom[5],		top[5]),
c(	left[3],		right[3],		bottom[5],		top[5]),
c(	left[4],		right[4],		bottom[5],		top[5]))

split.screen(m)

col.title <- c(
	"True tree nodes: bs ica",
	"True tree nodes: qj ica",
	"True tree nodes\nNOT found by bs: qj ica",
	"Nodes NOT in true tree\nfound by bootstraps")

for (i in 1:length(col.title)) {
	screen(i)
	par(mar = c(0,0,0,0))
	plot(1, axes=FALSE, type="n")
	text(x=1,y=1,col.title[i],cex=0.75)
}

row.title <- c(
	"b-d | random",
	"b-d | constant",
	"pec | random",
	"sym | random | long tips",
	"sym | random | short tips")

for (i in 1:length(row.title)) {
	screen(i+length(col.title))
	par(mar = c(0,0,0,0))
	plot(1, axes=FALSE, type="n")
	text(x=0.7,y=1,row.title[i],cex=0.75,srt=90)
}

## Process data
#input.files <- c(
#	"~/Dropbox/Projects_current/Jackknife_test/data_products/random_rates_model_100_trees_of_50_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/equal_rates_model_100_trees_of_50_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/pectinate_random_model_100_trees_of_50_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_model_100_trees_of_50_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_short_tips_model_100_trees_of_50_tips/all_scores/ALL.scores.csv"
#)

#input.files <- c(
#	"~/Dropbox/Projects_current/Jackknife_test/data_products/random_rates_model_50_trees_of_100_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/equal_rates_model_50_trees_of_100_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/pectinate_random_model_50_trees_of_100_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_model_50_trees_of_100_tips/all_scores/ALL.scores.csv"
#	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_short_tips_model_50_trees_of_100_tips/all_scores/ALL.scores.csv"
#)

input.files <- c(
	"~/Dropbox/Projects_current/Jackknife_test/data_products/random_rates_model_10_trees_of_500_tips/all_scores/ALL.scores.csv"
	,"~/Dropbox/Projects_current/Jackknife_test/data_products/equal_rates_model_10_trees_of_500_tips/all_scores/ALL.scores.csv"
	,"~/Dropbox/Projects_current/Jackknife_test/data_products/pectinate_random_model_10_trees_of_500_tips/all_scores/ALL.scores.csv"
	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_model_10_trees_of_500_tips/all_scores/ALL.scores.csv"
	,"~/Dropbox/Projects_current/Jackknife_test/data_products/balanced_random_short_tips_model_10_trees_of_500_tips/all_scores/ALL.scores.csv"
)

### may want to organize these primarily by treeshape/dataset instead of size. so, all sizes of pectinate/random 100% together, etc.  

for (i in 1:length(input.files)) {
	
	xx <- read.csv(input.files[i])
	
	# nodes not in true tree (but in bootstraps)
	aa <- xx[which(xx[,"in_true_tree"] == FALSE),]
	
	# nodes in true tree
	bb <- xx[which(xx[,"in_true_tree"] == TRUE),]
	
	# nodes in true tree but not in bootstraps
	cc <- bb[which(is.na(bb[,"b_ica"])),]
	
	# nodes in bootstraps and in true tree
	dd <- bb[which(!is.na(bb[,"b_ica"])),]

	nc <- length(col.title)
	nr <- length(row.title)

	screen(nr + nc*i + 1)
	par(mar=rep(bh,4))
	plot.ica.vs.node.depth(data.frame(dd[,"depth"],dd[,"b_ica"]), write.x.label=ifelse(i==length(input.files),TRUE,FALSE), write.y.axis=TRUE)

	screen(nr + nc*i + 2)
	par(mar=rep(bh,4))
	plot.ica.vs.node.depth(data.frame(dd[,"depth"],dd[,"j_ica"]), write.x.label=ifelse(i==length(input.files),TRUE,FALSE))
	
	screen(nr + nc*i + 3)
	par(mar=rep(bh,4))
	plot.ica.vs.node.depth(data.frame(cc[,"depth"],cc[,"j_ica"]), write.x.label=ifelse(i==length(input.files),TRUE,FALSE))
	
	screen(nr + nc*i + 4)
	par(mar=rep(bh,4))
	plot.j.ica.vs.b.ica.density(aa, ifelse(i==length(input.files),TRUE,FALSE))

}

### use to make a table showing when the j ica is greater than the b ica for nodes not in the true tree
#j.greater.than.b <- sum(a[,"j_ica"] > a[,"b_ica"]) / nrow(a)
#j.less.than.b <- sum(a[,"j_ica"] < a[,"b_ica"]) / nrow(a)
#text(paste("P(j < b) =",format(j.less.than.b,digits=3)), x=0.8, y=-1)
#text(paste("P(j > b) =",format(j.greater.than.b,digits=3)), x=-0.8, y=.95)
