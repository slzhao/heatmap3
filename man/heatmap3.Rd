\name{heatmap3}
\alias{heatmap3}
\title{heatmap3}
\usage{
  heatmap3(x, Rowv = NULL,
    Colv = if (symm) "Rowv" else NULL,
    distfun = function(x) as.dist(1 - cor(t(x), use = "pa")),
    balanceColor = F, ColSideLabs, RowSideLabs,
    showColDendro = T, showRowDendro = T,
    col = colorRampPalette(c("navy", "white", "firebrick3"))(1024),
    legendfun, method = "complete", ColAxisColors = 0,
    RowAxisColors = 0, hclustfun = hclust,
    reorderfun = function(d, w) reorder(d, w), add.expr,
    symm = FALSE, revC = identical(Colv, "Rowv"),
    scale = c("row", "column", "none"), na.rm = TRUE,
    margins = c(5, 5), ColSideColors, RowSideColors,
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc),
    labRow = NULL, labCol = NULL, main = NULL, xlab = NULL,
    ylab = NULL, keep.dendro = FALSE,
    verbose = getOption("verbose"), ...)
}
\arguments{
  \item{legendfun}{function used to generate legend in top
  left of the figure. If not specified, the color bar will
  be plotted. You can use any plot functions to generate
  your own legend. Or a function \code{\link{showLegend}}
  is also provided as a example.}

  \item{method}{the agglomeration method to be used by
  \code{\link{hclust}} function. This should be (an
  unambiguous abbreviation of) one of "ward", "single",
  "complete", "average", "mcquitty", "median" or
  "centroid".}

  \item{balanceColor}{logical indicating if the colors need
  to be balanced so that the median color will represent
  the 0 value. The default value is F.}

  \item{ColAxisColors}{integer indicating which coloum of
  ColSideColors will be used as colors for labels in coloum
  axis. The default value is 0, which means all coloum
  labels will be in black color.}

  \item{RowAxisColors}{integer indicating which coloum of
  RowSideColors will be used as colors for labels in row
  axis. The default value is 0, which means all row labels
  will be in black color.}

  \item{showColDendro}{logical indicating if the coloum
  dendrogram should be plotted (when Colv isn't NA).}

  \item{showRowDendro}{logical indicating if the row
  dendrogram should be plotted (when Rowv isn't NA).}

  \item{col}{specifying the colors, used in
  \code{\link{image}} function.}

  \item{cexRow,cexCol}{positive numbers, used as cex.axis
  in for the row or column axis labeling. The defaults
  currently only use number of rows or columns,
  respectively.}

  \item{labRow,labCol}{character vectors with row and
  column labels to use; these default to rownames(x) or
  colnames(x), respectively.}

  \item{main,xlab,ylab}{main, x- and y-axis titles;
  defaults to none.}

  \item{...}{additional arguments passed on to
  \code{\link{image}}.}

  \item{x}{numeric matrix of the values to be plotted. }

  \item{Rowv}{determines if and how the \emph{row}
  dendrogram should be computed and reordered.  Either a
  \code{\link{dendrogram}} or a vector of values used to
  reorder the row dendrogram or \code{\link{NA}} to
  suppress any row dendrogram (and reordering) or by
  default, \code{\link{NULL}}, see \sQuote{Details} below.}

  \item{Colv}{determines if and how the \emph{column}
  dendrogram should be reordered.  Has the same options as
  the \code{Rowv} argument above and \emph{additionally}
  when \code{x} is a square matrix, \code{Colv = "Rowv"}
  means that columns should be treated identically to the
  rows (and so if there is to be no row dendrogram there
  will not be a column one either).}

  \item{distfun}{function used to compute the distance
  (dissimilarity) between both rows and columns.  Defaults
  to \code{\link{dist}}.}

  \item{hclustfun}{function used to compute the
  hierarchical clustering when \code{Rowv} or \code{Colv}
  are not dendrograms.  Defaults to \code{\link{hclust}}.
  Should take as argument a result of \code{distfun} and
  return an object to which \code{\link{as.dendrogram}} can
  be applied.}

  \item{reorderfun}{\code{function(d, w)} of dendrogram and
  weights for reordering the row and column dendrograms.
  The default uses \code{\link{reorder.dendrogram}}.}

  \item{add.expr}{expression that will be evaluated after
  the call to \code{image}.  Can be used to add components
  to the plot.}

  \item{symm}{logical indicating if \code{x} should be
  treated \bold{symm}etrically; can only be true when
  \code{x} is a square matrix.}

  \item{revC}{logical indicating if the column order should
  be \code{\link{rev}}ersed for plotting, such that e.g.,
  for the symmetric case, the symmetry axis is as usual.}

  \item{scale}{character indicating if the values should be
  centered and scaled in either the row direction or the
  column direction, or none.  The default is \code{"row"}
  if \code{symm} false, and \code{"none"} otherwise.}

  \item{na.rm}{logical indicating whether \code{NA}'s
  should be removed.}

  \item{margins}{numeric vector of length 2 containing the
  margins (see \code{\link{par}(mar = *)}) for column and
  row names, respectively.}

  \item{ColSideColors}{(optional) character vector of
  length \code{ncol(x)} containing the color names for a
  horizontal side bar that may be used to annotate the
  columns of \code{x}.}

  \item{RowSideColors}{(optional) character vector of
  length \code{nrow(x)} containing the color names for a
  vertical side bar that may be used to annotate the rows
  of \code{x}.}

  \item{keep.dendro}{logical indicating if the
  dendrogram(s) should be kept as part of the result (when
  \code{Rowv} and/or \code{Colv} are not NA).}

  \item{verbose}{logical indicating if information should
  be printed.}
}
\value{
  The same return value as \code{\link{hclust}} function.
}
\description{
  The function heatmap3 is completely compatible with the
  original R function heatmap, and provides more new
  features.
}
\examples{
#gererate data
set.seed(123456789)
rnormData<-matrix(rnorm(1000), 40, 25)
rnormData[1:15, seq(6, 25, 2)] = rnormData[1:15, seq(6, 25, 2)] + 2
rnormData[16:40, seq(7, 25, 2)] = rnormData[16:40, seq(7, 25, 2)] + 4
colnames(rnormData)<-c(paste("Control", 1:5, sep = ""), paste(c("TrtA", "TrtB"), rep(1:10,each=2), sep = ""))
rownames(rnormData)<-paste("Probe", 1:40, sep = "")
ColSideColors<-c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10))
#A simple example
heatmap3(rnormData,ColSideColors=ColSideColors,showRowDendro=FALSE)
#A more detail example
ColSideColors<-cbind(Method1=c(rep("steelblue2",5), rep(c("lightgoldenrod"),20)),Method2=c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10)))
RowSideColors<-colorRampPalette(c("chartreuse4", "white", "firebrick"))(40)
heatmap3(rnormData,ColSideColors=ColSideColors,RowSideColors=RowSideColors,col=colorRampPalette(c("green", "black", "red"))(1024),ColAxisColors=1,RowAxisColors=1,legendfun=function() showLegend(legend=c("Control","Treatment","TrtA Treatment","TrtB Treatment"),col=c("steelblue2","lightgoldenrod","brown1","mediumpurple2")))
}

