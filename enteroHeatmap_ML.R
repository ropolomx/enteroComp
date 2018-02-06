# Installing and loading packages -----------------------------------------

# You can install packages by typing:
# install.packages("gplots")

library(gplots)
library(dplyr)
library(heatmaply)

# Read data ---------------------------------------------------------------

# You have provided multiple SNP distance matrices
# We want to load all of them into R

# This is how you would load only one file
# You assign it to a variable and then use the read.csv function to parse it and
# load it into R.

groel_ungapped <- read.csv('./data/SNP_table_groEL_ungapped_jan122018.csv')

# Nice! We have read in a SNP matrix as a data frame
# The data frame is one of the most important data structures in R
# A data frame is the most common way to store data in R
# A data frame has similarities to an Excel spreadsheet, or a table, but those 
# are very different things under the hood

# How do I know that this is a data frame?

# Using the class function, you will get the following information:

# class(groel_ungapped)
# [1] "data.frame"

# Let's now read the metadata file
# In this case, we are using a different function called read.table.
# This function reads files that are separated by a user-provided separator character
# In this case, that character is the tab (\t)

metadata <- read.table('data/sample_metadata.txt', 
                       sep="\t",
                       header=TRUE)

# Preprocessing data frames -----------------------------------------------

# In this particular format from Geneious, there are points in the data frame 
# that have NA values when they are read into R. This means that the values are
# empty, and R assigns them a missing value indicator ('Not available')
# We can replace those NA values by 0 for this particular format
# This is because the NA were generated for identical strains (and thus zero difference)

groel_ungapped[is.na(groel_ungapped)] <- 0

# Now, let's rename the first column of the groel_ungapped_meta

groel_ungapped <- groel_ungapped %>%
  rename(Iso = X)

# The reason for this is that we are going to use that column to merge it with the
# metadata file using that column

groel_ungapped_meta <- left_join(groel_ungapped, metadata, by="Iso")

# Let's select certain columns to keep for the metadata

groel_ungapped_meta <- groel_ungapped_meta %>%
  dplyr::select(Final.ID, everything()) %>%
  dplyr::select(-dplyr::contains("Iso")) %>%
  dplyr::select(-dplyr::contains("Location"))
                  

# Generate interactive heatmap --------------------------------------------

groel_ungapped_heatmap <- heatmaply(groel_ungapped_meta)

# Working in batch --------------------------------------------------------

# But how about loading multiple files at a time?

# Let's use some wizardry

# First, let's create a list of the names of the CSV files in this folder

allSNPfiles <- Sys.glob(file.path("*.csv"))

allSN1Pnames <- map(allSNPfiles, 
                   function(x) str_split(string = x, pattern = "_")) %>% 
  flatten() %>% 
  map(function(x) paste0(x[3],"_",x[4]))

allSNPdf <- lapply(allSNPfiles, function(x){
  df <- read.csv(x)
  df
})

allSNPdf <- lapply(allSNPdf, function(x){
  row.names(x) <- x$X
  x
})

allSNPdfDist <- lapply(allSNPdf, function(x){
  x <- x[,2:ncol(x)]
  dists <- as.dist(x)
  dists
})

allSNPdfHeatmaps <- lapply(allSNPdfDist, function(x){
  hm <- heatmaply(as.matrix(x))
  hm
})

# Generating static heatmaps with gplots ----------------------------------

entero_metadata <- read.csv()

groel_ungapped <- groel_ungapped[,2:ncol(groel_ungapped)]

groel_ungapped[is.na(groel_ungapped)] <- 0

dist_groel_ungapped <- as.dist(groel_ungapped_rep_na)

hc_dist_groel_ungapped <- hclust(dist_groel_ungapped, method = "complete")

dend_hc <- as.dendrogram(hc_dist_groel_ungapped)

heatmap.2(as.matrix(dist_groel_ungapped),
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

