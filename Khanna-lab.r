##SINGLE-CELL-SEQ########
########################################
# set working directory
chmod -R 775 as18818/
cd /gpfs/scratch/as18818/Khanna-lab/2025
# get on a screenscreen -S khanna

# get on the compute node
srun -c2 --partition=a100_short -t4000:00 --mem=70000 --pty /bin/bash
sinfo
srun --x11 --partition=a100_short --nodes=1 --mem=330GB -t1:00:00  --pty /bin/bash

module purge
module load r/4.0.3 
R# Install from CRAN 
install.packages("iCellR")
setwd("/your/download/directory")
library("iCellR")

##AGRREGATE-THE-DATA##
s1 <- load10x("RNAs_CMO302_1922/")
s2 <- load10x("RNAs_CMO310_2448/")

my.data <- data.aggregation(samples = c("s1","s2"),
	condition.names = c("P1","P2"))

head(my.data)
my.obj <- make.obj(my.data)
my.obj


"""> my.obj
###################################
,--. ,-----.       ,--.,--.,------. 
`--''  .--./ ,---. |  ||  ||  .--. ' 
,--.|  |    | .-. :|  ||  ||  '--'.' 
|  |'  '--'\   --. |  ||  ||  |   
`--' `-----' `----'`--'`--'`--' '--' 
###################################
An object of class iCellR version: 1.6.7
Raw/original data dimentions (rows,columns): 32285,4370
Data conditions in raw data: P1,P2 (1922,2448)
Row names: A030001D20Rik,A030003K21Rik,A030005K14Rik ...
Columns names: P1_AAACCCAAGACGGAAA.1,P1_AAACCCACAGCAGTGA.1,P1_AAACCCACAGCGTATT.1 ...
###################################
   QC stats performed:FALSE, PCA performed:FALSE
   Clustering performed:FALSE, Number of clusters:0
   tSNE performed:FALSE, UMAP performed:FALSE, DiffMap performed:FALSE
   Main data dimensions (rows,columns): 0,0
   Normalization factors:,... 
   Imputed data dimensions (rows,columns):0,0
############## scVDJ-seq ###########
VDJ data dimentions (rows,columns):0,0
############## CITE-seq ############
   ADT raw data  dimensions (rows,columns):0,0
   ADT main data  dimensions (rows,columns):0,0
   ADT columns names:... 
   ADT row names:... 
############## scATAC-seq ############
   ATAC raw data  dimensions (rows,columns):0,0
   ATAC main data  dimensions (rows,columns):0,0
   ATAC columns names:... 
   ATAC row names:... 
############## Spatial ###########
Spatial data dimentions (rows,columns):0,0
########### iCellR object ##########
"""
##PERFORM-SOME-QC##
my.obj <- qc.stats(my.obj)
stats.plot(my.obj,
	plot.type = "three.in.one",
	out.name = "UMI-plot",
	interactive = FALSE,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green")
dev.off()
##FILTER-CELLS##
my.obj <- cell.filter(my.obj,
	min.mito = 0,
	max.mito = 0.05,
	min.genes = 500,
	max.genes = 4000,
	min.umis = 0,
	max.umis = Inf)

# check to see how many cells are left.  
dim(my.obj@main.data)

'''> dim(my.obj@main.data)
[1] 32285  3900'''

##optional-down-sampling:
#for having the same number of cells for each condition.
#my.obj <- down.sample(my.obj)

##Normalize data##
my.obj <- norm.data(my.obj, 
     norm.method = "ranked.glsf",
     top.rank = 500) # best for scRNA-Seq

##Scale data (optional)
#my.obj <- data.scale(my.obj)

##Make a gene model for clustering
##This function will help you find a good number of genes to use for running PCA.
my.obj <- qc.stats(my.obj,which.data = "main.data")

#plot after filtering:
stats.plot(my.obj,
	plot.type = "three.in.one",
	out.name = "UMI-plot",
	interactive = FALSE,
	cell.color = "slategray3", 
	cell.size = 1, 
	cell.transparency = 0.5,
	box.color = "red",
	box.line.col = "green")
dev.off()

##FILTER-
# See model plot 
my.obj <- gene.stats(my.obj, which.data = "main.data")
head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])

make.gene.model(my.obj, my.out.put = "plot",
	dispersion.limit = 1.5, 
	base.mean.rank = 1500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	out.name = "gene.model")
head(my.obj@main.data,15)[1:4]
'''> head(my.obj@gene.data[order(my.obj@gene.data$numberOfCells, decreasing = T),])
        genes numberOfCells totalNumberOfCells percentOfCells  meanExp      SDs
27654  Tmsb4x          3893               3900       99.82051 73.93348 40.39624
594      Actb          3892               3900       99.79487 47.97247 28.52823
18763  Malat1          3892               3900       99.79487 83.74273 53.72605
6624     Fth1          3857               3900       98.89744 31.55277 47.72071
1814      B2m          3855               3900       98.84615 20.88300 17.38422
12777 Gm42418          3851               3900       98.74359 35.18469 62.57193
      condition
27654       all
594         all
18763       all
6624        all
1814        all
12777       all
'''
"SDs – The standard deviation of gene expression across cells, which indicates variability."

# Write the gene model data into the object
my.obj <- make.gene.model(my.obj, my.out.put = "data",
	dispersion.limit = 1.5, 
	base.mean.rank = 1500, 
	no.mito.model = T, 
	mark.mito = T, 
	interactive = F,
	out.name = "gene.model")

head(my.obj@gene.model)
"""The make.gene.model() function in iCellR is used to identify highly variable genes (HVGs), which are critical for clustering, dimensionality reduction, and downstream analysis. 
It filters genes based on expression levels and variability, creating a model of the most informative genes for analysis.
A gene is considered highly variable if its dispersion (variability) is greater than the threshold set here.
"""
#Perform Principal component analysis (PCA):
##Note:skip this step if you plan to do batch correction. For batch correction (sample alignment/harmonization/integration)
# If you run PCA (run.pca) there would be no batch alignment but if you run CPCA (using iba function) this would perform batch alignment and PCA after batch alignment. Example for batch alignment using iba function: 
# my.obj <- iba(my.obj,dims = 1:30, k = 10,ba.method = "CPCA", method = "gene.model", gene.list = my.obj@gene.model)

# run PCA in case no batch alignment is necessary
my.obj <- run.pca(my.obj, method = "gene.model", gene.list = my.obj@gene.model,data.type = "main")

opt.pcs.plot(my.obj)

length(my.obj@gene.model)
##Perform other dimensionality reductions (tSNE, UMAP, KNetL, PHATE, destiny, diffusion maps)
# KNetL is fundamentally more powerful.
# tSNE
my.obj <- run.pc.tsne(my.obj, dims = 1:10)

# UMAP
my.obj <- run.umap(my.obj, dims = 1:10)

# KNetL (for lager than 5000 cell use a zoom of about 400) 
# Because knetl has a very high resolution it's best to use a dim of 20 (this usually works best for most data)
my.obj <- run.knetl(my.obj, dims = 1:20, zoom = 300) # (Important note!) don't forget to set the zoom in the right range  
##try run zoom 300/110--do cluster using Knetl again


########################## IMPORTANT DISCLAIMER NOTE ###########################
            *** KNetL map is very dynamic with zoom and dims! ***
                 *** Therefore it needs to be adjusted! ***
# For data with less than 1000 cells use a zoom of about 5-50.
# For data with 1000-5000 cells use a zoom of about 50-200.
# For data with 5000-10000 cells use a zoom of about 100-300.
# For data with 10000-30000 cells use a zoom of about 200-500.
# For data with more than 30000 cells use a zoom of about 400-600.
# zoom 400 is usually good for big data but adjust for intended resolution.
# Lower number for zoom in and higher for zoom out (its reverse).
# dims = 1:20 is generally good for most data.
# other parameters are best as default.

#### Just like a microscope, you need to zoom to see the intended amount of details. 
#### Here we use a zoom of 100 or 110 but this might not be ideal for your data.
#### example: # my.obj <- run.knetl(my.obj, dims = 1:20, zoom = 400)
#### Because knetl has a very high resolution it's best to use a dim of 20 (this usually works best for most data)
###################################################################################

###Visualizing the results of dimensionality reductions before clustering
A= cluster.plot(my.obj,plot.type = "pca",interactive = F)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F)

library(gridExtra)
name <- "All-before-clustering.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
grid.arrange(A,B,C,D)
dev.off()

#####Clustering#####

##since the knel plot looks the best, we will proceed doing clustring based on Knetl.
##Generally peoplw use pca

# clustering based on KNetL

my.obj <- iclust(my.obj, sensitivity = 150, data.type = "knetl") 

#Visualize data clustering results
A= cluster.plot(my.obj,plot.type = "pca",interactive = F)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F)

name <- "grid-plot-after-clustering-zoom=300.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
grid.arrange(A,B,C,D)
dev.off()

#Look at the conditions:
# conditions 
A <- cluster.plot(my.obj,plot.type = "pca",col.by = "conditions",interactive = F,cell.size = 0.5)
B <- cluster.plot(my.obj,plot.type = "umap",col.by = "conditions",interactive = F,cell.size = 0.5)
C <- cluster.plot(my.obj,plot.type = "tsne",col.by = "conditions",interactive = F,cell.size = 0.5)
D <- cluster.plot(my.obj,plot.type = "knetl",col.by = "conditions",interactive = F,cell.size = 0.5)

name <- "gridPlot-conditions-zoom=300.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
grid.arrange(A,B,C,D)
dev.off()

###important

name <- "gridPlot.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
grid.arrange(A,B,C,D)
dev.off()

##save the object
save(my.obj, file = "my.obj.Robj")

##Renumbering the clusters, helps in better visualisation of heatmap

my.obj <- clust.ord(my.obj,top.rank = 500, how.to.order = "distance")
#my.obj <- clust.ord(my.obj,top.rank = 500, how.to.order = "random")

my.obj <- change.clust(my.obj, change.clust = 9, to.clust = 10)

#A= cluster.plot(my.obj,plot.type = "pca",interactive = F,cell.size = 0.5,cell.transparency = 1, anno.clust=T)
B= cluster.plot(my.obj,plot.type = "umap",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
C= cluster.plot(my.obj,plot.type = "tsne",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
D= cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.5,cell.transparency = 1,anno.clust=T)
name <- "AllClusts3.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
grid.arrange(B,B,C,D)
dev.off()

my.obj <- clust.ord(my.obj,top.rank = 500, how.to.order = "distance")


##@@@@@@@@@Find marker genes@@@@@@@@##

marker.genes <- findMarkers(my.obj,
	fold.change = 2,
	padjval = 0.1)

dim(marker.genes)
"""[1] 7785   20"""
head(marker.genes)
##Heatmap:

# find top genes
MyGenes <- top.markers(marker.genes, topde = 10, min.base.mean = 0.2,filt.ambig = F)
MyGenes <- unique(MyGenes)

# main data 
name <- "heat-map.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", conds.to.plot = NULL)

# sort cells and plot only one condition
name <- "heat-map-conditions.png"
png(name,width = 10, height = 10, units = 'in', res = 300)
heatmap.gg.plot(my.obj, gene = MyGenes, interactive = F, cluster.by = "clusters", data.type = "imputed", cell.sort = TRUE, conds.to.plot = c("p1"))

##Bubble Heatmap
png('heatmap_bubble_gg_genes.png', width = 10, height = 20, units = 'in', res = 300)
bubble.gg.plot(my.obj, gene = MyGenes, interactive = F, conds.to.plot = NULL, size = "Percent.Expressed",colour = "Expression")
dev.off()

#genes:
#genelist = c("CD3E","CD4","CD8A", "FOXP3","TNFRSF9", "NKG7","GNLY","GZMH","FCGR3A","CD19","MS4A1","CD69","CD169","LYZ","LY6G","CD15","CD68","CD1C","CTLA4","CCR7")
genelist = c("Ccr7","Cd14","Cd19","Cd3e","Cd4","Cd68","Cd8a","Ctla4","Epcam","Itgam","Fcgr3","Foxp3","Gnly","Kras","Fcgr2b","Ms4a1","Nkg7","Ptprc","S100a9","Tnfrsf9","Gzmb","Gzma","Gnly")
rm(list = ls(pattern="PL_"))

for(i in genelist){
####
    MyPlot <- gene.plot(my.obj, gene = i,
        interactive = F,
        cell.size = 0.1,
        plot.data.type = "knetl",
        data.type = "main",
        scaleValue = T,
        min.scale = 0,max.scale = 2.0,
        cell.transparency = 1)
####
    NameCol=paste("PL",i,sep="_")
    eval(call("<-", as.name(NameCol), MyPlot))
}

library(cowplot)
filenames <- ls(pattern="PL_")

B <- cluster.plot(my.obj,plot.type = "knetl",interactive = F,cell.size = 0.1,cell.transparency = 1,anno.clust=T)
filenames <- c("B",filenames)

png('genes_KNetL.png',width = 15, height = 12, units = 'in', res = 300)
plot_grid(plotlist=mget(filenames))
dev.off()

#Mice
#c("Ccr7","Cd14","Cd19","Cd3e","Cd4","Cd68","Cd8a","Cd1c","Ctla4","Epcam","Itgam","Fcgr3","Foxp3","Gnly","Kras","Fcgr2b","Ms4a1","Nkg7","Ptprc","S100a9","Tnfrsf9","Gzmb","Gzma","Gnly")
#HUMANS
#c("CD3E","CD4","CD8A", "FOXP3","TNFRSF9", "NKG7","GNLY","GZMH","FCGR3A","CD19","MS4A1","CD69","CD169","LYZ","LY6G","CD15","CD68","CD1C","CTLA4","CCR7")