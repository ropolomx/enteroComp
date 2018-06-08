# Installing and loading packages -----------------------------------------

# You can install packages by typing:
# install.packages("gplots")

library(gplots)
library(tidyverse) # We will be using dplyr, purrr, and stringr
library(heatmaply)
library(here)

# Read data ---------------------------------------------------------------

# You have provided multiple SNP distance matrices
# We want to load all of them into R

# This is how you would load only one file
# You assign it to a variable and then use the read.csv function to parse it and
# load it into R.

groel_gapped <- read.csv(
  here('data', 'SNP_table_groel_gapped_jan122018.csv'),
  stringsAsFactors = FALSE)

groel_gapped$X <- make.names(groel_gapped$X)

# Nice! We have read in a SNP matrix as a data frame
# The data frame is one of the most important data structures in R
# A data frame is the most common way to store data in R
# A data frame has similarities to an Excel spreadsheet, or a table, but those 
# are very different things under the hood

# How do I know that this is a data frame?

# Using the class function, you will get the following information:

# class(groel_gapped)
# [1] "data.frame"

# Let's now read the metadata file
# In this case, we are using a different function called read.table.
# This function reads files that are separated by a user-provided separator character
# In this case, that character is the tab (\t)

metadata <- read_csv('data/sample_metadata.csv')

metadata$`Final ID` <- na_if(metadata$`Final ID`, "")

metadata$Iso <- make.names(metadata$Iso)

# Pre-processing one dataframe -----------------------------------------------

# In this particular format from Geneious, there are points in the data frame 
# that have NA values when they are read into R. This means that the values are
# empty, and R assigns them a missing value indicator ('Not available')
# We can replace those NA values by 0 for this particular format
# This is because the NA were generated for identical strains (and thus zero difference)

groel_gapped[is.na(groel_gapped)] <- 0

# Now, let's rename the first column of the groel_ungapped_meta

groel_gapped <- groel_gapped %>%
  rename(Iso = X)

# The reason for this is that we are going to use that column to merge it with the
# metadata file using that column

groel_gapped_meta <- left_join(groel_gapped, metadata, by="Iso")

meta_in_groel <- intersect(metadata$Iso,groel_gapped_meta$Iso)

meta_in_groel <-
  metadata %>%
  filter(Iso %in% meta_in_groel)

# Let's select certain columns to keep for the metadata

groel_gapped_meta <- 
  groel_gapped_meta %>%
  dplyr::select(-Iso,-Location) %>%
  dplyr::select(`Final ID`, everything())
                  
# Generate one interactive heatmap ----------------------------------------

groel_gapped_heatmap <- heatmaply(groel_gapped_meta)

# Preliminary steps for generating heatmaps in batch ----------------------

# How about processing all the SNP distance files?
# Let's use some wizardry

# First, let's create a list of the names of the CSV files in this folder

allSNPfiles <- Sys.glob(file.path("./data/SNP*.csv"))

# Now, let's create a list with all the filename prefixes only (e.g. 16SatpApheS_gapped)

allSNPnames <- map(allSNPfiles, function(x) str_split(string = x, pattern = "_")) %>% 
  flatten() %>% 
  map(function(x) paste0(x[3],"_",x[4]))

# The code below will read all the CSV files in the allSNPfiles list and load them 
# into a list in the R Global environment

allSNPdf <- map(allSNPfiles, function(x){
  df <- read.csv(x)
  df
}) %>%
  set_names(nm = allSNPnames)


# Extract isolate names ---------------------------------------------------

# Let's fix the names so that we only have the isolate names

allSNPdf <- map(allSNPdf, function(x){
  df <- x %>%
    rename(Iso = X) %>%
    mutate(Iso = str_replace(Iso, "(_|\\s|Assembly|\\.).*$", "")) %>%
    mutate(Iso = str_replace(Iso, "\\.", "-"))
  names(df) <- str_replace(names(df), "(_|\\s|Assembly|\\.).*$", "")
  names(df) <- str_replace(names(df), "\\.", "-")
  df
})

# In this function we are replacing the NA values with zero

allSNPdf <- map(allSNPdf, function(x){
  x[is.na(x)] <- 0
  x
})

# We might have data duplicates 
# Let's check for those

duplicate_positions <- map(allSNPdf, ~ which(duplicated(.x$Iso))) %>%
  discard( ~ length(.x) == 0)

duplicate_positions

# $atpa_gapped
# [1] 26
# 
# $pheS_gapped
# [1] 176

allSNPdf$atpa_gapped <- allSNPdf$atpa_gapped[-26,-27]

allSNPdf$pheS_gapped <- allSNPdf$pheS_gapped[-176,-177]

allSNPdf <- map(allSNPdf, function(x) {
  df <- x %>%
    mutate(Iso = make.names(Iso))
  df
})


# Merge isolate metadata with SNP distances -------------------------------

allSNPdfMeta <- map(allSNPdf, function(x){
  meta <- left_join(x, metadata, by="Iso")
  meta
})

allSNPdfMeta <- map(allSNPdfMeta, function(x){
  meta <- dplyr::select(x, `Final ID`, everything()) %>%
    dplyr::select(-dplyr::contains("Iso")) %>%
    dplyr::select(-dplyr::contains("Location")) %>%
    rename(Species=`Final ID`)
})

allSNPdfMetaSub <- discard(allSNPdfMeta, str_detect(names(allSNPdfMeta), "ungapped"))

# Generate all heatmaps ---------------------------------------------------

allSNPdfHeatmaps <- imap(allSNPdf_MetaSub, function(x,y){
  hm <- heatmaply(x, 
                  dendrogram = "both",
                  plot_method = "plotly",
                  main = y,
                  # margins = c(50,32,NA,11),
                  fontsize_row = 8,
                  fontsize_col = 11,
                  column_text_angle = 45,
                  key.title = "SNP distance",
                  # colorbar_xpos = -1,
                  colorbar_ypos = 2,
                  showticklabels = FALSE
         )
  hm
})

# Export heatmaps as PNG --------------------------------------------------

allSNPdfHeatmaps %>%
  iwalk( ~ export(
    p = .x,
    file = paste0(here("snp_heatmaps/"),
      .y,
      ".png"),
    vwidth = 1246,
    vheight = 814
  ))


# Export heatmaps as HTML -------------------------------------------------

allSNPdfHeatmaps %>%
  iwalk(
     ~ htmlwidgets::saveWidget(as_widget(.x), paste0(here("snp_heatmaps/"),.y,".html"))
  )


# Cluster analysis for Comparing Partitions: one example ------------------

# Purpose: extract clusters at all thresholds for analysis with the Comparing 
# Partitions tool (http://www.comparingpartitions.info/)

# Example with one dataset

groel_gapped <- groel_gapped[,2:ncol(groel_gapped)]

groel_gapped_rep_na <- lapply(groel_gapped, function(x){
  x <- str_replace_na(x, replacement = 0)
  x
})

groel_gapped_rep_na <- do.call("rbind", groel_gapped_rep_na)

dist_groel_gapped <- as.dist(groel_gapped_rep_na)

all_dists <- unique(as.vector(dist_groel_gapped))

hc_dist_groel_gapped <- hclust(dist_groel_gapped, method = "complete")

# clusters_groel_strains <- metadata[match(hc_dist_groel_gapped$labels, metadata$Iso),]

# write.csv(clusters_groel_strains, 'reference_clusters.csv')

# hc_dist_groel_gapped$labels <- hc_dist_groel_gapped$labels[hc_dist_groel_gapped$order]

clusters_groel <- cutree(hc_dist_groel_gapped, h = all_dists)

write.csv(clusters_groel, file = here('cluster_analysis', 'clusters_at_all_SNP_thresholds.csv'))

clusters_groel_species <- cutree(hc_dist_groel_gapped, k = 7)

# dend_hc <- as.dendrogram(hc_dist_groel_gapped)

clusters_groel_sum <- apply(clusters_groel, 2, max)

clusters_groel_sum <- data.frame(
  distances = as.numeric(names(clusters_groel_sum)),
  clusters = clusters_groel_sum
) %>%
  arrange(distances)

# Cluster analysis for Comparing Partitions: in batch ---------------------

all_dists <-
  allSNPdfMetaSub %>%
  map( ~ (.x[, 2:ncol(.x)] %>%
      as.dist(.)))

all_clusts <-
  all_dists %>%
  map( ~ (hclust(.x, method = "complete")))

all_flat_clusters <-
  all_clusts %>%
  map( ~ (cutree(.x, k = 7))) # Check for high and low thresholds

all_flat_clusters_df <-
  all_flat_clusters %>%
  map( ~ enframe(.x, name = "Iso", value = "Cluster"))

# all_flat_clusters_df <-
#   all_flat_clusters_df %>%
#   map(~ .x[match(metadata$Iso, .x$Iso), ] %>%
#       na.omit())

all_flat_clusters_df <-
  all_flat_clusters_df %>%
  map_dfr(~ left_join(.x, metadata,by="Iso"), .id = 'Scheme')


butter_spread <- spread(all_flat_clusters_df, key = Scheme, value = Cluster)
butter_spread <- butter_spread[match(metadata$Iso, butter_spread$Iso),]
butter_spread <- butter_spread %>% na.omit()
write.csv(butter_spread, here('cluster_analysis', 'all_seven_clusters.csv', row.names = FALSE, quote = FALSE)
