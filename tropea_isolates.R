# Script for analyzing the presence/absence data provided by 
# Erica Tropea to Haley Sanderson


# Loading packages --------------------------------------------------------

library(tidyverse)
library(readxl)
library(pheatmap)
library(here)


# Read Excel file ---------------------------------------------------------

isolate_path <- here('data','Erica_Tropea_isolate_information.xlsx')

binary_profiles <- isolate_path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = isolate_path)



# Extracting numerical columns for building matrices ----------------------

# Let's take advantage of the fact that the columns
# with the presence/absence data are of a different format 
# than the rest of the columns (i.e. double/numeric)
# R read the binary profile column as character, so that also helps

# > head(binary_profiles$Fecal)
# # A tibble: 6 x 10
# Isolates     iutA  ccdB clpXET  hra1   phd  traT Phylogroup Fingerprint Source 
#   <chr>       <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <chr>      <chr>       <chr>  
# 1 Cat 1_1         0     0      0     0     0     1 B1         000001      Cat    
# 2 Cat 1_2         0     0      0     0     0     1 B1         000001      Cat    
# 3 Cat 1_3         0     0      0     0     0     1 B1         000001      Cat    
# 4 Cat 3_1         0     0      0     0     0     0 B1         000000      Cat    
# 5 Chicken 3_1     0     0      1     0     1     1 B1         001011      Chicken
# 6 Chicken 3_2     0     0      1     0     1     1 B1         001011      Chicken

binary_matrices <-
  binary_profiles %>% 
  map( ~ .x %>% select_if(is.double) %>% as.matrix(.))


# Preparing metadata ------------------------------------------------------


# Assign isolate names as rownames of matrices

rownames(binary_matrices$Fecal) <- binary_profiles$Fecal$Isolates
rownames(binary_matrices$`Well water`) <- binary_profiles$`Well water`$Isolates

# Create annotation row data frames for adding metadata to the heatmap

annotation_rows <-
  binary_profiles %>%
  map( ~ .x %>% select(Source) %>% as.data.frame(.))

rownames(annotation_rows$Fecal) <- binary_profiles$Fecal$Isolates


# Plotting heatmaps -------------------------------------------------------

pheatmap(binary_matrices$Fecal,
         color = viridis::cividis(n=2),
         legend = FALSE,
         border_color = NA,
         cellwidth = 90,
         cellheight = 6.5,
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         clustering_method = "average",
         # annotation_row = annotation_rows$Fecal,
         fontsize_row = 8
         )

# Change annotation colorscheme

# First, investigate how many unique sources we have

length(unique(annotation_rows$Fecal$Source))
# 15

ann_colors <- list(
  
)

