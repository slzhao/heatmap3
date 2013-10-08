###############################################################################
# Data: 2012-09-19
# Author: Shilin Zhao (zhaoshilin@gmail.com)
###############################################################################

##' heatmap3
##' 
##' The function heatmap3 is completely compatible with the original R function heatmap, and provides more new features.
##' 
##' 
##' 
##' @inheritParams stats::heatmap
##' @param legendfun function used to generate legend in top left of the figure. If not specified, the color bar will be plotted. You can use any plot functions to generate your own legend. Or a function \code{\link{showLegend}} is also provided as a example.
##' @param method the agglomeration method to be used by \code{\link{hclust}} function. This should be (an unambiguous abbreviation of) one of "ward", "single", "complete", "average", "mcquitty", "median" or "centroid".
##' @param balanceColor logical indicating if the colors need to be balanced so that the median color will represent the 0 value. The default value is F.
##' @param ColAxisColors integer indicating which coloum of ColSideColors will be used as colors for labels in coloum axis. The default value is 0, which means all coloum labels will be in black color.
##' @param RowAxisColors integer indicating which coloum of RowSideColors will be used as colors for labels in row axis. The default value is 0, which means all row labels will be in black color.
##' @param showColDendro logical indicating if the coloum dendrogram should be plotted (when Colv isn't NA).
##' @param showRowDendro logical indicating if the row dendrogram should be plotted (when Rowv isn't NA).
##' @param col specifying the colors, used in \code{\link{image}} function.
##' @param cexRow,cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
##' @param labRow,labCol character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
##' @param main,xlab,ylab main, x- and y-axis titles; defaults to none.
##' @param ... additional arguments passed on to \code{\link{image}}.
##' @export
##' @return The same return value as \code{\link{hclust}} function.
##' @examples #gererate data
##' set.seed(123456789)
##' rnormData<-matrix(rnorm(1000), 40, 25)
##' rnormData[1:15, seq(6, 25, 2)] = rnormData[1:15, seq(6, 25, 2)] + 2
##' rnormData[16:40, seq(7, 25, 2)] = rnormData[16:40, seq(7, 25, 2)] + 4
##' colnames(rnormData)<-c(paste("Control", 1:5, sep = ""), paste(c("TrtA", "TrtB"), rep(1:10,each=2), sep = ""))
##' rownames(rnormData)<-paste("Probe", 1:40, sep = "")
##' ColSideColors<-c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10))
##' #A simple example
##' heatmap3(rnormData,ColSideColors=ColSideColors,showRowDendro=FALSE)
##' #A more detail example
##' ColSideColors<-cbind(Method1=c(rep("steelblue2",5), rep(c("lightgoldenrod"),20)),Method2=c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10)))
##' RowSideColors<-colorRampPalette(c("chartreuse4", "white", "firebrick"))(40)
##' heatmap3(rnormData,ColSideColors=ColSideColors,RowSideColors=RowSideColors,col=colorRampPalette(c("green", "black", "red"))(1024),ColAxisColors=1,RowAxisColors=1,legendfun=function() showLegend(legend=c("Control","Treatment","TrtA Treatment","TrtB Treatment"),col=c("steelblue2","lightgoldenrod","brown1","mediumpurple2")))
heatmap3<-function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
		distfun = function(x) as.dist(1 - cor(t(x),use="pa")),balanceColor=F, ColSideLabs,RowSideLabs,showColDendro=T,showRowDendro=T,col=colorRampPalette(c("navy", "white", "firebrick3"))(1024),legendfun,method="complete",ColAxisColors=0,RowAxisColors=0, hclustfun = hclust, reorderfun = function(d, 
				w) reorder(d, w), add.expr,symm = FALSE, revC = identical(Colv, 
				"Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
		margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
				1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
		labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
		verbose = getOption("verbose"),...) 
{
	scale <- if (missing(scale)) 
				"none"
			else match.arg(scale)
	if (is.data.frame(x)) {x<-as.matrix(x)}
	if (!missing(ColSideColors)) {
		if (is.vector(ColSideColors)) {
			ColSideColors<-cbind(ColSideColors)
		}
	}
	if (!missing(RowSideColors)) {
		if (is.vector(RowSideColors)) {
			RowSideColors<-cbind(RowSideColors)
		}
	}
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("'x' must be a numeric matrix")
	nr <- di[1L]
	nc <- di[2L]
	if (nr <= 1 || nc <= 1) 
		stop("'x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2L) 
		stop("'margins' must be a numeric vector of length 2")
	doRdend <- !identical(Rowv, NA)
	doCdend <- !identical(Colv, NA)
	if (!doRdend && identical(Colv, "Rowv")) 
		doCdend <- FALSE
	if (is.null(Rowv)) 
		Rowv <- rowMeans(x, na.rm = na.rm)
	if (is.null(Colv)) 
		Colv <- colMeans(x, na.rm = na.rm)
	if (doRdend) {
		if (inherits(Rowv, "dendrogram")) 
			ddr <- Rowv
		else {
			hcr <- hclustfun(distfun(x),method=method)
			ddr <- as.dendrogram(hcr)
			if (!is.logical(Rowv) || Rowv) 
				ddr <- reorderfun(ddr, Rowv)
		}
		if (nr != length(rowInd <- order.dendrogram(ddr))) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else rowInd <- 1L:nr
	if (doCdend) {
		if (inherits(Colv, "dendrogram")) 
			ddc <- Colv
		else if (identical(Colv, "Rowv")) {
			if (nr != nc) 
				stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
			ddc <- ddr
		}
		else {
			hcc <- hclustfun(distfun(if (symm) 
										x
									else t(x)))
			ddc <- as.dendrogram(hcc)
			if (!is.logical(Colv) || Colv) 
				ddc <- reorderfun(ddc, Colv)
		}
		if (nc != length(colInd <- order.dendrogram(ddc))) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else colInd <- 1L:nc
	x <- x[rowInd, colInd]
	labRow <- if (is.null(labRow)) 
				if (is.null(rownames(x))) 
					(1L:nr)[rowInd]
				else rownames(x)
			else labRow[rowInd]
	labCol <- if (is.null(labCol)) 
				if (is.null(colnames(x))) 
					(1L:nc)[colInd]
				else colnames(x)
			else labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 1L, sd, na.rm = na.rm)
		x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
	}
	else if (scale == "column") {
		x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
		sx <- apply(x, 2L, sd, na.rm = na.rm)
		x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
	}
	lmat <- rbind(c(NA, 3), 2:1)
#	lwid <- c(if (doRdend) 1 else 0.05, 4)
	lwid <- c(1, 4)
#	lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
#			4)
	lhei <- c( 1 + if (!is.null(main)) 0.2 else 0,4)
	if (!missing(ColSideColors)) {
		if (!is.character(ColSideColors) || nrow(ColSideColors) != 
				nc) 
			stop("'ColSideColors' must be a character vector or matrix of length ncol(x)")
		lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
		lhei <- c(lhei[1L], 0.2, lhei[2L])
	}
	if (!missing(RowSideColors)) {
		if (!is.character(RowSideColors) || nrow(RowSideColors) != 
				nr) 
			stop("'RowSideColors' must be a character vector or matrix of length nrow(x)")
		lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
						1), lmat[, 2] + 1)
		lwid <- c(lwid[1L], 0.2, lwid[2L])
	}
	lmat<-lmat+1
	lmat[is.na(lmat)] <- 0
	lmat[1,1]<-1
	if (verbose) {
		cat("layout: widths = ", lwid, ", heights = ", lhei, 
				"; lmat=\n")
		print(lmat)
	}
	dev.hold()
	on.exit(dev.flush())
	op <- par(no.readonly = TRUE)
	on.exit(par(op), add = TRUE)
	
	#balanceColor
	if (balanceColor) {
		if (abs(max(x,na.rm=T))>=abs(min(x,na.rm=T))) {
			cut.off<-round(quantile(1:length(col),probs=1-(abs(max(x,na.rm=T))+abs(min(x,na.rm=T)))/(2*abs(max(x,na.rm=T)))))
			col<-col[cut.off:length(col)]
		} else {
			cut.off<-round(quantile(1:length(col),probs=(abs(max(x,na.rm=T))+abs(min(x,na.rm=T)))/(2*abs(min(x,na.rm=T)))))
			col<-col[1:cut.off]
		}
	}
	
	layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
	if (!missing(legendfun)) {
		par(mar = c(0, 0, 0, 0))
		legendfun()
	} else {
		par(mar = c(5, 1, 1, 0))
		dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), 
				length = length(col))
		dummy.z <- matrix(dummy.x, ncol = 1)
		image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", 
				col = col,cex.axis=cexCol,xlab="")
	}
	if (!missing(RowSideColors)) {
		par(mar = c(margins[1L], 0, 0, 0.5))
		if (revC) {
			rsc = RowSideColors[rev(rowInd),,drop=F]
		} else {
			rsc = RowSideColors[rowInd,,drop=F]
		}
		rsc.colors = matrix()
		rsc.names = names(table(rsc))
		rsc.i = 1
		for (rsc.name in rsc.names) {
			rsc.colors[rsc.i] = rsc.name
			rsc[rsc == rsc.name] = rsc.i
			rsc.i = rsc.i + 1
		}
		rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
		image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
		
		if (missing(RowSideLabs)) {
			RowSideLabs<-colnames(RowSideColors)
		}
		if (dim(rsc)[2]==1) {
			axis(1, 0, RowSideLabs, las = 2, tick = FALSE)
		} else {
			axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), RowSideLabs, 
					las = 2, tick = FALSE)
		}
	}
	if (!missing(ColSideColors)) {
		par(mar = c(0.5, 0, 0, margins[2L]))
		csc = ColSideColors[colInd,,drop=F]
		csc.colors = matrix()
		csc.names = names(table(csc))
		csc.i = 1
		for (csc.name in csc.names) {
			csc.colors[csc.i] = csc.name
			csc[csc == csc.name] = csc.i
			csc.i = csc.i + 1
		}
		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
		image(csc, col = as.vector(csc.colors), axes = FALSE)
		if (missing(ColSideLabs)) {
			ColSideLabs<-colnames(ColSideColors)
		}
		if (dim(csc)[2]==1) {
			axis(4, 0, ColSideLabs,las = 2, tick = FALSE)
		} else {
			axis(4, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), ColSideLabs,las = 2, tick = FALSE)
		}
	}
	par(mar = c(margins[1L], 0, 0, margins[2L]))
	if (!symm || scale != "none") 
		x <- t(x)
	if (revC) {
		iy <- nr:1
		if (doRdend)
			ddr <- rev(ddr)
		x <- x[, iy]
	}
	else iy <- 1L:nr
	image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
					c(0, nr), axes = FALSE, xlab = "", ylab = "", col=col,...)
	if (!missing(ColSideColors) & ColAxisColors!=0) {
		mtext(1, at=1L:nc, text = labCol, las = 2, line = 0.5,cex = cexCol,col=ColSideColors[colInd,ColAxisColors])
	} else {
		axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0,cex.axis = cexCol)
	}
	
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1L] - 1.25)
	if (!missing(RowSideColors) & RowAxisColors!=0) {
		mtext(4, at=iy, text = labRow, las = 2, line = 0.5,cex = cexRow,col=RowSideColors[rowInd,RowAxisColors])
	} else {
		axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow)
	}
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2L] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	par(mar = c(margins[1L], 0, 0, 0))
	if (doRdend & showRowDendro) 
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	else frame()
	par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
	if (doCdend & showColDendro) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}	else if (!is.null(main)) 
		frame()
	if (!is.null(main)) {
		par(xpd = NA)
		title(main, cex.main = 1.5 * op[["cex.main"]])
	}
	invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
							doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}

##' showLegend
##' 
##' The function showLegend is an example for generating legend in the figure of heatmap3 function. You can use your any plot functions to generate your own legend.
##' 
##' 
##' 
##' @inheritParams graphics::legend
##' @param lwd the line widths for lines appearing in the legend.
##' @param ... additional arguments passed on to \code{\link{legend}}
##' @export
##' @examples RowSideColors<-rep("steelblue2",nrow(mtcars))
##' RowSideColors[c(4:6,15:17,22:26,29)]<-"lightgoldenrod"
##' RowSideColors[c(1:3,19:21)]<-"brown1"
##' heatmap3(mtcars,scale="col",margins=c(2,10),RowSideColors=RowSideColors,legendfun=function() showLegend(legend=c("European","American","Japanese"),col=c("steelblue2","lightgoldenrod","brown1"),cex=1.5))
showLegend<-function(legend=c("Group A","Group B"),lwd=3,cex=1.1,col=c("red","blue"),...) {
	plot(0,xaxt="n",bty="n",yaxt="n",type="n",xlab="",ylab="")
	legend("topleft",legend=legend,lwd=lwd,col=col,bty="n",cex=cex,...)
}
