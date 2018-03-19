
# Load packages -----------------------------------------------------------

library(readxl)
library(dplyr)
library(purrr)
library(stringr)
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

# Merge E.cass and E.gall

ast_profiles <- lapply(ast_profiles, function(x) {
  x <- x %>%
    mutate(Species = str_replace(
      Species,
      "^E. (casseliflavus$|gallinarum$|casseliflvaus$)",
      "E. casseliflavus\\/E. gallinarum"
    ))
  x
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

# Hamming datasets with merged data ---------------------------------------

ast_hamming_BAF <- list(
  ast_hamming_species$`BAF FE Species binary`,
  ast_hamming_species$`BAF PE Species binary`
)

ast_hamming_BAF <- do.call("rbind", ast_hamming_BAF)

ast_hamming_CAS <- list(
  ast_hamming_species$`CAS FE Species binary`,
  ast_hamming_species$`CAS PE Species binary`
)

ast_hamming_CAS

# Generate heatmaps --------------------------------------------------------

# Preliminary step: order Species columns the same way for all datasets
# This is to get the same colour scheme in all the plots

ast_profiles$`BAF FE Species binary`$Species <- 
  factor(as.factor(ast_profiles$`BAF FE Species binary`$Species),
         levels = c("E. faecium", 
              "E. faecalis", 
              "E. hirae",
              "E. casseliflavus/E. gallinarum"
              )
  )

ast_profiles$`BAF PE Species binary`$Species <-
  factor(as.factor(ast_profiles$`BAF PE Species binary`$Species),
         levels = c("E. faecium",
                    "E. faecalis",
                    "E. hirae",
                    "E. casseliflavus/E. gallinarum",
                    "E. mundtii",
                    "E. saccharolyticus"
                    )
         )

ast_profiles$`CAS FE Species binary`$Species <-
  factor(as.factor(ast_profiles$`CAS FE Species binary`$Species),
         levels = c("E. faecium",
                    "E. faecalis",
                    "E. hirae",
                    "E. casseliflavus/E. gallinarum"
                    )
         )


ast_profiles$`CAS PE Species binary`$Species <-
  factor(as.factor(ast_profiles$`CAS PE Species binary`$Species),
         levels = c("E. faecium",
                    "E. faecalis",
                    "E. hirae",
                    "E. casseliflavus/E. gallinarum"
                    )
         )

# Presence/absence heatmaps (binary patterns only)

# Will attempt to generate separately due to color schemes
# TODO: customize rowside colors

custom_heatmaply <- function(profile, title, colScale){
  heatmaply(profile[c(2,4:15)],
            dendrogram="both",
            seriate = "mean",
            plot_method = "plotly",
            hclust_method = "average",
            margins = c(45,32,NA,11),
            fontsize_row = 8,
            fontsize_col = 13,
            main = title,
            RowSideColors = colScale,
            hide_colorbar = TRUE,
            labRow=profile$Iso
  )}

ast_heatmap_baf_fe <- custom_heatmaply(
  ast_profiles$`BAF FE Species binary`,
  "BAF FE Species",
  c("light blue",
    "light green",
    "green",
    "pink")
  )

# Map function to generate All heatmaps

ast_heatmaps <- imap(ast_profiles, function(x,y){
  hm <- heatmaply(x[c(2,4:15)],
                  dendrogram="both",
                  seriate = "mean",
                  plot_method = "plotly",
                  hclust_method = "average",
                  main = y,
                  margins = c(45,32,NA,11),
                  fontsize_row = 7,
                  fontsize_col = 10,
                  # row_side_palette = c(
                  hide_colorbar = TRUE,
                  labRow=x$Iso
                  )
  hm
})

# Saving directly to the 

iwalk(ast_profiles, function(x,y){
  heatmaply(x[c(2,4:15)],
            dendrogram="both",
            seriate = "mean",
            plot_method = "plotly",
            hclust_method = "average",
            main = y,
            margins = c(50,35,NA,11),
            hide_colorbar = TRUE,
            labRow=x$Iso,
            file = paste("heatmaps/",y,".png"),
            width = 1200,
            height = 1000
            )
  # browseURL(paste("heatmaps/",y,".png"))
})

# Heatmaps of Hamming distance matrices

ast_hamming_heatmaps <- map2(ast_hamming_species, ast_profiles, function(x,y){
  hm <- heatmaply(x,
                  dendrogram="both",
                  hclust_method = "mcquitty",
                  seriate = "mean",
                  plot_method="plotly",
                  margins = c(50,NA,NA,NA),
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





