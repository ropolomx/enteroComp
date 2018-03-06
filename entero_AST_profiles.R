
# Load packages -----------------------------------------------------------

library(readxl)
library(dplyr)
library(purrr)
library(heatmaply)
library(UpSetR)

# Read data ---------------------------------------------------------------

ast_path <- 'data/Isolates_with_Binary_AMR phenotype_Jan292018.xlsx'

# Read all sheets in Excel file

ast_profiles <- ast_path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = ast_path)

# Pre-process data --------------------------------------------------------

# Drop first data frame from list
ast_profiles <- ast_profiles[-1]

# Remove all rows with NA values
ast_profiles <- lapply(ast_profiles, function(x){
  na.omit(x)
})


# Tabulate binary profiles ------------------------------------------------

ast_binary_freqs <- ast_profiles %>%
  map(~ sort(table(.x[16]), decreasing = TRUE))

# Calculate Hamming distance matrices -------------------------------------

# Create function for calculating Hamming distance of binary AST data
# This function uses matrix multiplication

hamming <- function(x){
  D <- (1-x) %*% t(x)
  D + t(D)
}

# Calculating the Hamming distance for each AST binary profile matrix

ast_hamming <- ast_profiles %>%
  map(~ hamming(as.matrix(.x[4:15])))

# ast_hamming_e1071 <- ast_profiles[2:5] %>%
#   map(~ hamming.distance(as.matrix(.x[4:15])))

# Create list with Species names for heatmaps

ast_hamming_species <- map2(ast_profiles, ast_hamming, ~ data.frame(.x$Species,.y))

# Generate heatmaps --------------------------------------------------------

# Presence/absence heatmaps (binary patterns only)

ast_heatmaps <- map(ast_profiles, function(x){
  hm <- heatmaply(x[c(2,4:15)],
                  dendrogram="both",
                  labRow=x$Iso
                  )
  hm
})

# Heatmaps of Hamming distance matrices

ast_hamming_heatmaps <- map2(ast_hamming_species, ast_profiles, function(x,y){
  hm <- heatmaply(x,
                  dendrogram="both",
                  hclust_method = "mcquitty",
                  seriate = "mean",
                  plot_method="plotly",
                  margins = c(30,NA,NA,NA),
                  
                  # scale_fill_gradient_fun = scale_fill_gradient(low="black", high="red")
                  labRow = y$Iso
                  )
  hm
})


# Generate UpSet figures --------------------------------------------------

binaries <- map(ast_profiles, function(x){
  bin_profile <- x$Binary
  bin_profile
})

names(binaries) <- c("BAF FE", "BAF PE", "CAS FE", "CAS PE", "Biomass")
  
upset(fromList(binaries),
      order.by = "freq",
      mainbar.y.label = "Shared Binary Profiles",
      # empty.intersections = TRUE,
      point.size = 3.2,
      text.scale=c(2,2,2,1.5,2))





