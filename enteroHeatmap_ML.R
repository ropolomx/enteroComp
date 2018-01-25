# Installing and loading packages -----------------------------------------

# You can install packages by typing:
# install.packages("gplots")

library(gplots)
library(heatmaply)

# Read data ---------------------------------------------------------------

# You have provided multiple SNP distance matrices
# We want to load all of them into R

# This is how you would load only one file
# You assign it to a variable and then use the read.csv function to parse it.

groel_ungapped <- read.csv('SNP table groEL ungapped jan122018.csv')

row.names(groel_ungapped) <- groel_ungapped$X

groel_ungapped <- groel_ungapped[,2:ncol(groel_ungapped)]

groel_ungapped_dist <- as.dist(groel_ungapped)

heatmaply(as.matrix(groel_ungapped_dist))

# Using the class function, you will get the following information:

# class(groel_ungapped)
# [1] "data.frame"

# Nice! We have read in a SNP matrix as a data frame

# But how about loading multiple files at a time?

# Let's use some wizardry

allSNPfiles <- Sys.glob(file.path("*.csv"))

allSNPdf <- lapply(allSNPfiles, function(x){
  df <- read.csv(x)
  df
})

allSNPdf <- lapply(allSNPdf, function(x){
  row.names(x) <- x$X
  x
})

# Preparing data for generating heatmap -----------------------------------

groel_ungapped <- groel_ungapped[,2:ncol(groel_ungapped)]

groel_ungapped[is.na(groel_ungapped)] <- 0

dist_groel_ungapped <- as.dist(groel_ungapped_rep_na)

hc_dist_groel_ungapped <- hclust(dist_groel_ungapped, method = "complete")

dend_hc <- as.dendrogram(hc_dist_groel_ungapped)

heatmap.2(as.matrix(as.disdist_groel_ungapped),
dendrogram = "both", # plot both dendrograms for rows and columns
#scale = "row", # Centers and scales data and calculates z-scores in either the row or column direction (or none)
key.title = "SNP Distance",
srtCol = 90, # Angle of column labels
cexCol = 1, # Modify to change font size of column labels
cexRow = 0.4, # Modify to change font size of column labels
Rowv = dend_hc, # Plot dendrogram object for the rows
#Colv = dend_c, # Plot dendrogram object for the columns
trace = "none" # Good that you removed this! 
#RowSideColors = colour_blocks
#lmat=rbind(c(4, 3), c(2, 1)), lhei=c(0.3, 2), lwid=c(1.5, 4)
# These are options to modify the layout. We can leave them commented out for now.
)

# Gapped

groel_gapped <- read.csv('SNP table groel gapped jan122018.csv')

row.names(groel_gapped) <- groel_gapped$X

groel_gapped <- groel_gapped[,2:ncol(groel_gapped)]

groel_gapped_rep_na <- lapply(groel_gapped, function(x){
  x <- str_replace_na(x, replacement = 0)
  x
})

groel_gapped_rep_na <- do.call("rbind", groel_gapped_rep_na)

dist_groel_gapped <- dist(groel_gapped_rep_na)

hc_dist_groel_gapped <- hclust(dist_groel_gapped, method = "complete")

dend_hc <- as.dendrogram(hc_dist_groel_gapped)

# Heatmap with ggplot2 ----------------------------------------------------

heatmap.2(as.matrix(dist_groel_gapped),
dendrogram = "both", # plot both dendrograms for rows and columns
#scale = "row", # Centers and scales data and calculates z-scores in either the row or column direction (or none)
key.title = "SNP Distance",
srtCol = 90, # Angle of column labels
cexCol = 0.4, # Modify to change font size of column labels
cexRow = 0.3, # Modify to change font size of column labels
Rowv = dend_hc, # Plot dendrogram object for the rows
#Colv = dend_c, # Plot dendrogram object for the columns
trace = "none" # Good that you removed this! 
#RowSideColors = colour_blocks
#lmat=rbind(c(4, 3), c(2, 1)), lhei=c(0.3, 2), lwid=c(1.5, 4)
# These are options to modify the layout. We can leave them commented out for now.
)

# Heatmap with heatmaply --------------------------------------------------

# Use heatmaply

allSNPdfDist <- lapply(allSNPdf, function(x){
  as.dist(x)
})
