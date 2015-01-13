plot.proportions <- function(d, breaks, cutoff, write.x.labels=FALSE, write.y.labels=FALSE, label='') {
	
	yticks = c(0, 0.2, 0.4, 0.6, 0.8, 1);
	xticks = c(-1, -0.5, 0, 0.5, 1);
	
	plot.new();
	plot.window(c(-1, 1), c(0, 1.1));
	if (write.x.labels) {
		axis(1, at=xticks, cex.axis=0.75, hadj=0.5, padj=-1.75, tcl=-0.2)
		mtext(side=1, cex=0.75, text=paste(label,"ICA",sep=" "), line=1.5)
	} else {
		axis(1, at=xticks, labels=NA, cex.axis=0.75, hadj=0.5, padj=-1.75, tcl=-0.2)
	}
	
	if (write.y.labels) {
		axis(2, at=yticks, cex.axis=0.75, hadj=0.5, padj=1.5, tcl=-0.2)
		mtext(side=2, cex=0.75, text="frequency", line=1.5)
	} else {
		axis(2, at=yticks, labels=NA, cex.axis=0.75, hadj=0.5, padj=1.5, tcl=-0.2)
	}
	
#	m <- max(d[1,] + d[2,])
	m <- sum(d[1,] + d[2,])
	p = d / m

	w <- 1 / ncol(p);
	z <- w*3/4
	l <- breaks - z;
	r <- breaks + z;
	
	t1 <- p[1,]
	t2 <- p[1,] + p[2,]
	rect(l,0,r,t1,col=c1,border=NA);
	rect(l,t1,r,t2,col=c2,border=NA);
	box();

	ps <- p[1,] + p[2,]
	v <- p[2,] / ps > cutoff;
	x <- p[1,] / ps > cutoff;

	if (sum(v[!is.na(v)]) > 0) {
		text(l[v], t2[v], adj=c(-z*5,-z*3), labels="*", col=c2);
	}
	
	if (sum(x[!is.na(x)]) > 0) {
		text(l[x], t2[x], adj=c(-z*5,-z*3), labels="*", col=c1)
	}
}

draw.plot <- function(stat, label='') {
	
	dev.new(width=w, height=h)
	close.screen(all=TRUE)
	plot.new()
	split.screen(m)
	
	for (i in 1:length(col.title)) {
		screen(i)
		par(mar = c(0,0,0,0))
		plot(1, axes=FALSE, type="n")
		text(x=1,y=1,col.title[i],cex=0.75)
	}
	
	for (i in 1:length(row.title)) {
		screen(i+length(col.title))
		par(mar = c(0,0,0,0))
		plot(1, axes=FALSE, type="n")
		text(x=0.7,y=1,row.title[i],cex=0.75,srt=90)
	}

	for (i in 1:length(input.files)) {
	
		xx <- read.csv(input.files[i])
		
		# nodes not in true tree (but in bootstraps)
		a <- xx[which(xx[,"in_true_tree"] == FALSE),][,stat]
		a <- a[which(!is.na(a))]
		
		# nodes in true tree
		b <- xx[which(xx[,"in_true_tree"] == TRUE),][,stat]
		b <- b[which(!is.na(b))]

		d <- rbind(sapply(1:(length(breaks)-1), function(j) {sum(a[a > breaks[j]] <= breaks[j+1])}),
				   sapply(1:(length(breaks)-1), function(j) {sum(b[b > breaks[j]] <= breaks[j+1])}));
	
		screen(k + i)
		par(mar=rep(bh,4))
		plot.proportions(d, breaks[2:length(breaks)], cutoff, write.x.labels=(i %% nrow == 0), write.y.labels=(i <= ncol), label=label);
	}
}

w <- 8
h <- 7.5
ncol <- 3
nrow <- 3

# margins
om <- 0.005 # 0.01 						# outer
lm <- 0.1 								# left
rm <- 1 - om		 					# right
tm <- 1 - om * w / h					# top
bm <- 0.04 * w / h						# bottom

bh <- 0.03								# border - horizontal
bv <- bh - om 							# border - vertical
t <- 2*bh								# top header height
cw <- (rm - lm - (ncol-1)*bh) / ncol			# plot cell height
ch <- (tm - bm - t - (nrow-1)*bv) / nrow		# plot cell width

top <- c(
	1+bm-t-0.5*bv-0*ch-om,
	1+bm-t-1.5*bv-1*ch-om,
	1+bm-t-2.5*bv-2*ch-om)

bottom <- c(
	tm+bm-t-0.5*bv-1*ch,
	tm+bm-t-1.5*bv-2*ch,
	tm+bm-t-2.5*bv-3*ch)

left <- c(
	lm+0*cw+0*bh,
	lm+1*cw+1*bh,
	rm-cw)

right <- c(
	lm+1*cw+0*bh,
	lm+2*cw+1*bh,
	rm)

#	left			right			bottom			top
m <- rbind(

# column labels
c(	left[1],		right[1],		tm+bm-t,			tm),
c(	left[2],		right[2],		tm+bm-t,			tm),
c(	left[3],		right[3],		tm+bm-t,			tm),

# row labels
c(	om,				lm,				bottom[1],		top[1]),
c(	om,				lm,				bottom[2],		top[2]),
c(	om,				lm,				bottom[3],		top[3]),

# plots
c(	left[1],		right[1],		bottom[1],		top[1]),
c(	left[1],		right[1],		bottom[2],		top[2]),
c(	left[1],		right[1],		bottom[3],		top[3]),
c(	left[2],		right[2],		bottom[1],		top[1]),
c(	left[2],		right[2],		bottom[2],		top[2]),
c(	left[2],		right[2],		bottom[3],		top[3]),
c(	left[3],		right[3],		bottom[1],		top[1]),
c(	left[3],		right[3],		bottom[2],		top[2]),
c(	left[3],		right[3],		bottom[3],		top[3]))

col.title <- c(
	"complete",
	"\"phylogenetic\"",
	"uniform 50%")

row.title <- c(
	"b-d | random",
	"sym | random | long tips",
	"sym | random | short tips")

k <- length(col.title) + length(row.title)

basedir <- "~/Dropbox/Projects_current/Jackknife_test/data_products/"
input.files <- c(
	# complete sampling
	paste(basedir,"random_rates_model_5_trees_of_1000_tips/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_model_5_trees_of_1000_tips/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_short_tips_model_5_trees_of_1000_tips/all_scores/ALL.scores.csv",sep="")
	
	# beta subsampling
	,paste(basedir,"random_rates_model_5_trees_of_1000_tips_SUBSAMPLED_beta/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_model_5_trees_of_1000_tips_SUBSAMPLED_beta/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_short_tips_model_5_trees_of_1000_tips_SUBSAMPLED_beta/all_scores/ALL.scores.csv",sep="")

	# uniform subsampling
	,paste(basedir,"random_rates_model_5_trees_of_1000_tips_SUBSAMPLED_uniform/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_model_5_trees_of_1000_tips_SUBSAMPLED_uniform/all_scores/ALL.scores.csv",sep="")
	,paste(basedir,"balanced_random_short_tips_model_5_trees_of_1000_tips_SUBSAMPLED_uniform/all_scores/ALL.scores.csv",sep="")
);

breaks <- c(-1.01,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0);
cutoff <- 0.95;

# for ica
c1 <- '#ffcc77'
c2 <- '#77ccff'
draw.plot("j_ica", "QJ")

# for bootstraps
c1 <- '#ff55cc'
c2 <- '#99dd44'
draw.plot("b_ica", "BS")