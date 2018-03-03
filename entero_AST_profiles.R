
# Load packages -----------------------------------------------------------

library(readxl)
library(dplyr)
library(heatmaply)

# Read data ---------------------------------------------------------------

ast_path <- 'data/Isolates_with_Binary_AMR phenotype_Jan292018.xlsx'

ast_profiles <- ast_path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = ast_path)

# Pre-process data --------------------------------------------------------

# Remove all rows with NA values
ast_profiles <- lapply(ast_profiles, function(x) {
  na.omit(x)
})

# Generate plots ----------------------------------------------------------

ast_heatmaps <- map(ast_profiles[2:5], function(x) {
  hm <- heatmaply(x[c(2,4:15)],
                  dendrogram="both",
                  labRow=x$Iso)
  hm
})

