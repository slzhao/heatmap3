heatmap3
============
* [Introduction](#Introduction)
* [Download and install](#download)
* [Example](#example)
* [Reproduce the figures](#reproduce)

<a name="Introduction"/>
# Introduction #
Heat map and clustering are used frequently in expression analysis studies for data visualization and quality control. Simple clustering and heat map can be produced from the “heatmap” function package in R. However, the “heatmap” function lacks functionality and customizability, preventing it from generating advanced heat maps and dendrograms. To tackle the limitations of “heatmap” function, we have developed an R package “heatmap3” which significantly improves the original “heatmap” function by adding several more powerful and convenient features. “heatmap3” packages allows user to produce highly customizable state of art heatmap and dendrogram. The heatmap3 package is developed based on the “heatmap” function in R language and it is completely compatible with it. And the new features of “heatmap3” include highly customizable legend and side annotation, a wider range of color selections, new labeling features which allow user to define multiple layers of phenotype variables and automatically conduct association test based on the phenotypes provided. Additional features such as different agglomeration methods for estimating distance between two samples are also added for clustering.    

<a name="download"/>
# Download and install #
You can install heatmap3 package in R from [github](https://github.com/slzhao/heatmap3/) by following codes:

	library(devtools)
	install_github("heatmap3", user="slzhao")

<a name="example"/>
# Example #
After you have installed heatmap3 package. You can enter R and use following R codes to see the examples for it.
	
	#Load MultiRankSeq package
	library(heatmap3)
	#View help files
	?heatmap3
	#The examples of heatmap3 and other functions
	example(heatmap3)
	example(showLegend)
	example(showAnn)
	example(colByValue)

<a name="reproduce"/>
# Reproduce the figures #

We used the TCGA BRCA data to generate two example figures, which were used in our paper.

First you will need to download the count data. We generate the count data based on the TCGA BRCA samples. You can use the following codes in Linux system to download the count data. Or you can download it in [github release page](https://github.com/slzhao/heatmap3/releases/tag/example).

	wget https://github.com/slzhao/heatmap3/releases/download/example/allSample_edgeR_result.csv
	wget https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples.csv
	wget https://github.com/slzhao/heatmap3/releases/download/example/BRCA_30Samples_clinic.csv

Here are the R codes to generate the figures.

	#Assume you have already installed heatmap3 package. And the 3 example file were download in current working directory. You can use the following codes in R to generate the figures.

    #Prepare expression data
    counts<-read.csv("BRCA_30Samples.csv",header=T,row.names=1)
    #Prepare column side annotation
    clinic<-read.csv("BRCA_30Samples_clinic.csv",header=T,row.names=1)
    #Prepare row side color bar annotation
    edgeR_result<-read.csv("allSample_edgeR_result.csv",header=T,row.names=1)
    temp1<-(edgeR_result$logFC)
    temp2<--log10(edgeR_result$FDR)
    temp1<-colByValue(as.matrix(temp1),range=c(-4,4),col=colorRampPalette(c('chartreuse4','white','firebrick'))(1024))
    temp2<-colByValue(as.matrix(temp2),range=c(0,5),col=heat.colors(1024))
    colGene<-cbind(temp1,temp2)
    row.names(colGene)<-row.names(edgeR_result)
    colnames(colGene)<-c("log2FC","-Log10P")
    
    #Generate Figure1
    #counts, colGene and clinic were read throught the csv file
    ##Assume counts has counts information for each gene, colGene has the colors for each gene, clinic has the clinic information for each sample
    temp<-apply(counts,1,sd)
    selectedGenes<-rev(order(temp))[1:500]
    heatmap3(counts[selectedGenes,],labRow="",margin=c(7,0),RowSideColors=colGene[selectedGenes,],ColSideCut=0.85,ColSideAnn=clinic,ColSideFun=function(x) showAnn(x),ColSideWidth=1.2,balanceColor=T)

    #Generate Figure2
    heatmap3(counts,topN=c(500,3000,nrow(counts)),Rowv=NA,labRow="",margin=c(7,0),RowSideColors=colGene,ColSideCut=0.85,ColSideAnn=clinic,ColSideFun=function(x) showAnn(x),ColSideWidth=1.2,balanceColor=T)

Here is the environment (including the version of packages).

    R version 3.0.2 (2013-09-25)
    Platform: x86_64-w64-mingw32/x64 (64-bit)

    locale:
    [1] LC_COLLATE=English_United States.1252 
    [2] LC_CTYPE=English_United States.1252   
    [3] LC_MONETARY=English_United States.1252
    [4] LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.1252    

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] heatmap3_1.0.3

    loaded via a namespace (and not attached):
    [1] fastcluster_1.1.13 tools_3.0.2       
