
# Load packages -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(gplots)
library(heatmaply)
library(UpSetR)
library(here)

# Read data ---------------------------------------------------------------

ast_path <- here('data','Isolates_with_Binary_AMR phenotype_Jan292018.xlsx')

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
  ast_profiles_BAF$FE
  # rename(AMPI = AMPR) # Depends on which version of data 

ast_profiles_BAF <- do.call("rbind", ast_profiles_BAF)

ast_profiles_BAF <- 
  ast_profiles_BAF %>%
  mutate(Location = str_replace(Location, "BAF\\s+",""))

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
  mutate(Location = str_replace(Location, "CAS\\s+",""))

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


# Generate binary heatmaps ------------------------------------------------

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

all_species <- map(ast_profiles, ~ unique(as.character(.x$Species))) %>% 
  flatten_chr(.) %>% 
  unique(.)

colorscale <- data.frame(
  Species = all_species,
  color = c(
    "forest green", 
    "red",
    "royal blue",
    "light green",
    "magenta",
    "orange")
  )


assign_colors <- function(x) {
  species_col <- x %>% select(Species)
  merged_colors <- left_join(species_col, colorscale, by = "Species") %>%
    select(-Species)
  # merged_colors <- as.character(merged_colors$color)
}

all_colors <- 
  ast_profiles %>%
  map(
    ~ assign_colors(.x)
  )

# custom_heatmaply <- function(profile, colScale){
#   heatmaply(profile[c(2,4:15)],
#             dendrogram="both",
#             seriate = "mean",
#             plot_method = "plotly",
#             hclust_method = "average",
#             margins = c(45,32,NA,11),
#             fontsize_row = 8,
#             fontsize_col = 13,
#             main = title,
#             RowSideColors = colScale,
#             hide_colorbar = TRUE,
#             labRow=profile$Iso
#   )}

# ast_heatmap_baf_fe <- custom_heatmaply(
#   ast_profiles$`BAF FE Species binary`,
#   "BAF FE Species",
#   c("light blue",
#     "light green",
#     "green",
#     "pink")
#   )

# Map function to generate All heatmaps

ast_profiles_merged <- map2(
  ast_profiles,
  all_colors,
  ~ cbind(.x,.y) %>% select(color, everything())
)

all_colors <- 
  all_colors %>%
  map(~ as.character(.x$color))

binary_heatmap <- function(profile){
  hm <- heatmaply_na(profile[c(3,5:16)],
    dendrogram="both",
    seriate = "mean",
    plot_method = "ggplot",
    hclust_method = "average",
    margins = c(45,32,NA,11),
    fontsize_row = 7,
    fontsize_col = 10,
    RowSideColors = profile$color,
    hide_colorbar = TRUE,
    labRow=profile$Iso
  )
  hm
}

ast_heatmaps <- map(
  ast_profiles_merged,
  # all_colors,
  ~ binary_heatmap(.x)
)


build_ast_dend <- function(profile){
  profile_mat <- as.matrix(profile[,5:16])
  rownames(profile_mat) <- profile$Iso
  profile_dist <- dist(profile_mat, method = "manhattan")
  profile_hc <- hclust(profile_dist, method = "average")
  profile_dend <- reorder(as.dendrogram(profile_hc), rowMeans(profile_mat))
  profile_dist_trans <- dist(t(profile_mat), method = "manhattan")
  profile_hc_trans <- hclust(profile_dist_trans, method = "average")
  profile_dend_trans <- reorder(as.dendrogram(profile_hc_trans), colMeans(profile_mat))
  return(
    list(
    "mat" = profile_mat,
    "row_dend" = profile_dend, 
    "col_dend" = profile_dend_trans
    )
  )
}

# baf_pe_dends <- build_ast_dend(ast_profiles_merged$`BAF PE Species binary`)
# 
# ast_hm_BAF_PE_mat <- as.matrix(ast_profiles_merged$`BAF PE Species binary`[,5:16])
# rownames(ast_hm_BAF_PE_mat) <- ast_profiles_merged$`BAF PE Species binary`$Iso
# ast_hm_BAF_PE_dist <- dist(ast_hm_BAF_PE_mat, "manhattan")
# ast_hm_BAF_PE_hc <- hclust(ast_hm_BAF_PE_dist, method = "average")
# ast_hm_BAF_PE_dend <- reorder(as.dendrogram(ast_hm_BAF_PE_hc),rowMeans(ast_hm_BAF_PE_mat))
# ast_hm_BAF_PE_dist_trans <- dist(t(ast_hm_BAF_PE_mat), "manhattan")
# ast_hm_BAF_PE_hc_trans <- hclust(ast_hm_BAF_PE_dist_trans, method = "average")
# ast_hm_BAF_PE_dend_trans <- reorder(as.dendrogram(ast_hm_BAF_PE_hc_trans), colMeans(ast_hm_BAF_PE_mat))

# baf_pe_colors <- as.character(ast_profiles_merged$`BAF PE Species binary`$color)

binary_heatmap.2 <- function(mat, row_dend, col_dend, colscale,leg,leg_col) {
  heatmap.2(mat,
    Rowv = row_dend,
    Colv = col_dend,
    trace = "none", 
    col = rev(viridis(n = 2, alpha = 1, begin = 0, end= 1, option = "magma")), 
    RowSideColors = colscale,
    lwid = c(2.5,4.0,1.5),
    lhei = c(1.5,5),
    cexRow = 0.75,
    # margins = c(2,2.5),
    # RowSideColors = as.character(ast_profiles_merged$`BAF PE Species binary`$color),
    key = FALSE
  )
  legend(list(x=0.0003023889,y=0.99999),
    legend = leg,
    fill = leg_col)
}

all_species_vec <- map(ast_profiles, ~ unique(as.character(.x$Species)) %>% sort(.))
colorscale <- colorscale %>% arrange(Species)
color_maps <- map(all_species_vec, ~ colorscale$color[colorscale$Species %in% .x])
color_maps <- map(color_maps, ~ as.character(.x))

ast_binary_mats <- 
  ast_profiles_merged %>%
  map(~ build_ast_dend(.x))


pwalk(
  list(ast_binary_mats, 
  all_colors,
  all_species_vec,
  color_maps),
  function(a,b,c,d) binary_heatmap.2(
    a$mat,
    a$row_dend,
    a$col_dend, 
    b,
    c,
    d
  )
) 

binary_heatmap.2(ast_binary_mats$`BAF FE Species binary`$mat,
  ast_binary_mats$`BAF FE Species binary`$row_dend,
  ast_binary_mats$`BAF FE Species binary`$col_dend,
  all_colors$`BAF FE Species binary`,
  all_species_vec$`BAF FE Species binary`,
  color_maps$`BAF FE Species binary`)

heatmap.2(baf_pe_dends$mat,
  Rowv = baf_pe_dends$row_dend,
  Colv = baf_pe_dends$col_dend,
  trace = "none", 
  col = viridis(n = 256, alpha = 1, begin = 0, end= 1, option = "viridis"), 
  RowSideColors = baf_pe_colors,
  # RowSideColors = as.character(ast_profiles_merged$`BAF PE Species binary`$color),
  key = FALSE
)




ast_hm_BAF_PE <- heatmaply(ast_profiles_merged$`BAF PE Species binary`[c(3,5:16)],
  dendrogram="both",
  seriate = "mean",
  plot_method = "ggplot",
  hclust_method = "average",
  margins = c(45,32,NA,11),
  fontsize_row = 7,
  fontsize_col = 10,
  row_side_colors = as.character(ast_profiles_merged$`BAF PE Species binary`$color),
  # hide_colorbar = TRUE,
  labRow = ast_profiles_merged$`BAF PE Species binary`$Iso
)
  


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

hamming_heatmap <- function(x, row_palette){
  heatmaply(
    x,
    dendrogram="both",
    # hclust_method = "mcquitty",
    # seriate = "mean",
    plot_method="plotly",
    fontsize_row = 8,
    fontsize_col = 8,
    RowSideColors = row_palette,
    margins = c(45,NA,NA,NA)
  )
}

# TODO: generate dataframe with colour assignment

hamming_palettes <- list(
  "BAF" = c(
    "dark orange",
    "light orange",
    "red",
    "pink",
    "green",
    "light green",
    "dark blue",
    "light blue"
    ),
  "CAS" = c(
    "dark orange",
    "light orange",
    "green",
    "light green",
    "blue",
    "light blue"
  )
)

ast_hamming_hm <- map2(
  list(
  "BAF" = ast_hamming_BAF, 
  "CAS" = ast_hamming_CAS
  ),
  hamming_palettes,
  ~ hamming_heatmap(.x, .y)
)

heatmaply(ast_hamming_BAF, 
  plot_method = "plotly", 
  fontsize_col = 6,
  fontsize_row = 6
)

heatmaply(ast_hamming_CAS, 
  plot_method = "plotly", 
  fontsize_col = 6,
  fontsize_row = 6
)

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
