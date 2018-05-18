# Load packages -----------------------------------------------------------

library(broom)

# Read Excel --------------------------------------------------------------

cog_file <- read_excel("COGstatsMay182018.xlsx", sheet = 3)

amr <- cog_file$AMR

# Perform models and obtain summary for all COGs and Virulence ------------

cog_models <- lapply(cog_file[,c(4:25,27)], function(x) lm(x ~ amr))
cog_models_tidy <- lapply(cog_models, tidy)
cog_models_tidy_df <- do.call("rbind", cog_models_tidy)

# Dataframe with all results combined
cog_models_tidy_df$Category <- row.names(cog_models_tidy_df)


# Obtain more information from the models  --------------------------------

# Will return information suchas R-squared, etc.

cog_models_glance <- lapply(cog_models, glance)
cog_models_glance_df <- do.call("rbind", cog_models_glance)
cog_models_glance_df$Category <- row.names(cog_models_glance_df)
