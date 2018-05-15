
# Load packages -----------------------------------------------------------

library(readxl)
library(tidyverse)
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

# Generate two lists, one for BAF data and the other for CAS data

ast_profiles_BAF <- list(
  ast_profiles$`BAF FE Species binary`,
  ast_profiles$`BAF PE Species binary`
) %>%
  set_names(nm = c("FE", "PE"))

ast_profiles_BAF$FE <- 
  ast_profiles_BAF$FE %>%
  rename(AMPI = AMPR)

ast_profiles_BAF <- do.call("rbind", ast_profiles_BAF)

ast_profiles_BAF <- 
  ast_profiles_BAF %>%
  mutate(Location = str_replace(Location, "BAF ",""))

ast_profiles_CAS <- list(
  ast_profiles$`CAS FE Species binary`,
  ast_profiles$`CAS PE Species binary`
) %>%
  set_names(nm = c ("FE", "PE"))

ast_profiles_CAS$FE <- 
  ast_profiles_CAS$FE %>%
  rename(TEIC = TERC)

ast_profiles_CAS <- do.call("rbind", ast_profiles_CAS)

ast_profiles_CAS <- 
  ast_profiles_CAS %>%
  mutate(Location = str_replace(Location, "CAS ",""))

# Calculate Hamming distances on BAF and CAS lists

# Retrieve Isolate names, then manually split by sewage treatment plant type

# There has to be a much better way to extract the isolate character vectors

isolate_names <- 
  ast_profiles %>%
  map(
    ~ select(.x,Iso) %>% 
      as_vector(.) %>% 
      unname(.)
  )

isolate_names_BAF <- list(
  isolate_names$`BAF FE Species binary`,
  isolate_names$`BAF PE Species binary`
) %>%
  set_names(nm = c("FE", "PE"))

isolate_names_CAS <- list(
  isolate_names$`CAS FE Species binary`,
  isolate_names$`CAS PE Species binary`
) %>%
  set_names(nm = c("FE", "PE"))


# Tabulate binary profiles ------------------------------------------------

ast_binary_freqs <- 
  ast_profiles %>%
  map(~ sort(table(.x[16]), decreasing = TRUE))

ast_binary_freqs_BAF <- sort(table(ast_profiles_BAF$Binary), decreasing = TRUE)

ast_binary_freqs_CAS <- sort(table(ast_profiles_CAS$Binary), decreasing = TRUE)

# Calculate Hamming distance matrices -------------------------------------

# Create function for calculating Hamming distance of binary AST data
# This function uses matrix multiplication

hamming <- function(x){
  D <- (1-x) %*% t(x) # Matrix multiplication
  D + t(D)
}

# Calculating the Hamming distance for each AST binary profile matrix

ast_hamming <- ast_profiles %>%
  map(~ hamming(as.matrix(.x[4:15])))

ast_hamming_BAF <- hamming(as.matrix(ast_profiles_BAF[4:15]))

ast_hamming_CAS <- hamming(as.matrix(ast_profiles_CAS[4:15]))
  
# Another option from e1071 package:
# ast_hamming_e1071 <- ast_profiles[2:5] %>%
#   map(~ hamming.distance(as.matrix(.x[4:15])))

# Create list with Species names for heatmaps

ast_hamming_species <- map2(
  ast_profiles, 
  ast_hamming, 
  ~ data.frame(.x$Species,.y)
)

# Function to assign the isolate names of each dataset as 
# column name

assign_names <- function(x,y){
  names(x) <- c("Species","Location",y)
  row.names(x) <- y # Isolate name y-axis,
  # x <- select(x, Species, Isolate, everything())
  x
}

ast_hamming_BAF <- data.frame(ast_profiles_BAF$Species, 
                              ast_profiles_BAF$Location,
                              ast_hamming_BAF)

ast_hamming_BAF <- assign_names(ast_hamming_BAF, ast_profiles_BAF$Iso)

ast_hamming_CAS <- data.frame(ast_profiles_CAS$Species, 
                              ast_profiles_CAS$Location,
                              ast_hamming_CAS)

ast_hamming_CAS <- assign_names(ast_hamming_CAS, ast_profiles_CAS$Iso)


# Merge Hamming datasets using tidy methods -------------------------------

# Might want consider using map2 here (and in other sections of the script)

# Extract BAF datasets

ast_hamming_BAF <- list(
  ast_hamming_species$`BAF FE Species binary`,
  ast_hamming_species$`BAF PE Species binary`
) %>%
  set_names(nm = c("FE","PE"))

# Function to assign the isolate names of each dataset as 
# column name

assign_names <- function(x,y){
  names(x) <- c("Species",y)
  x$Isolate_y <- y # Isolate name y-axis,
  x <- select(x, Species, Isolate_y, everything())
  x
}

ast_hamming_BAF <- map2(
  ast_hamming_BAF,
  isolate_names_BAF,
  ~ assign_names(.x,.y)
)

# Might have to tidy data first, then merge, then untidy again

ast_hamming_BAF_tidy <-
  ast_hamming_BAF %>%
  map_dfr(
    ~ gather(.x,
      key = Isolate_x,
      value = hamming,
      -Species, -Isolate_y
    ),
    .id = "Location"
  )

# Let's seriously consider map2 here to avoid repetition, 
# or to work on a list of lists

ast_hamming_CAS <- list(
  ast_hamming_species$`CAS FE Species binary`,
  ast_hamming_species$`CAS PE Species binary`
) %>%
  set_names(c("FE","PE"))

ast_hamming_CAS <- map2(
  ast_hamming_CAS,
  isolate_names_CAS,
  ~ assign_names(.x,.y)
)

ast_hamming_CAS_tidy <-
  ast_hamming_CAS %>%
  map_dfr(
    ~ gather(.x,
      key = Isolate_x,
      value = hamming,
      -Species, -Isolate_y
    ),
    .id = "Location"
  )
  
# TODO: Address the following warning messages

# Warning messages:
# 1: In bind_rows_(x, .id) : Unequal factor levels: coercing to character
# 2: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector
# 3: In bind_rows_(x, .id) :
#   binding character and factor vector, coercing into character vector


# Converting Hamming datasets back to matrix format -----------------------

# Purpose: requested by Haley to plot as heatmap

ast_hamming_BAF_wide <-
  ast_hamming_BAF_tidy %>%
  spread(key = Isolate_x, value = hamming)

# Result: Many columns with NA values because not all strains were retrieved
# from both environments

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

# Saving directly to file

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

# Heatmaps of Hamming distances -------------------------------------------

ast_hamming_heatmaps <- map2(
  ast_hamming_species, 
  ast_profiles, 
  function(x,y){
  hm <- heatmaply(x,
                  dendrogram="both",
                  hclust_method = "mcquitty",
                  seriate = "mean",
                  plot_method="plotly",
                  margins = c(45,NA,NA,NA),
                  # scale_fill_gradient_fun = scale_fill_gradient(low="black", high="red")
                  labRow = y$Iso,
                  labCol = y$Iso
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
