# Title: Analysis of host matrices: traits, phylogeny, and interaction similarity
# Authors: Miguel Vinagre de Frutos, Gabriella Lima Tabet Cruz, Cecilia Siliansky de Andreazzi
# Date: July 07, 2025
# Description: Generates and analyses three matrices among mammal hosts: 
# (1) Gower distance based on traits or shared parasites, 
# (2) Jaccard distance based on traits or shared parasites, 
# (3) phylogenetic distance, and 
# (4) similarity based on shared parasites.

# LOAD LIBRARIES
library(tidyverse)    # Data manipulation and visualization
library(readxl)       # Reading Excel files
library(ggpubr)       # For boxplots with statistics
library(GGally)       # For ggpairs()
library(FactoMineR)   # For mixed factor analysis
library(corrplot)     # Correlation matrix plotting
library(car)          # For vif() function
library(FD)      # Gower distance calculation
library(ape)          # Phylogenetic analysis
library(vegan)        # Mantel tests
library(ecodist)      # Distance-based regression (MRM)

# SET PATHS
localDir <- "."
data.directory <- file.path(localDir, "Data")

# PREPARING DATA

# LOAD INTERACTION DATA
interaction_traits <- read_excel(file.path(data.directory, "IPMO-Portugal-Spain-France-Host-Parasite.xlsx"))
interaction_traits <- interaction_traits %>% 
  filter(Domesticated_Captive_or_Wild == "Wild") %>% 
  filter(ZoonoticStatus == "1") %>% 
  filter(ParStatus == "species") %>% 
  filter(HostSpecies != "not_found") %>% 
  filter(Country %in% c("Spain", "Portugal")) %>% 
  distinct()

# Correct outdated or mismatched species names
interaction_traits$HostSpecies <- gsub("Mustela_vison", "Neovison_vison", interaction_traits$HostSpecies)

str(interaction_traits)

# OPTIONAL: ANALYSIS WITH A DATA SUBSET BY HOST ORDER
# Host taxonomic classification
# wildOrder <- read.csv(file.path(data.directory, "List-species-silvestres-brazil-iberian-peninsula.csv"), stringsAsFactors = TRUE)
# wildOrder$HostSpecies <- gsub(" ", "_", wildOrder$HostSpecies)
# wildOrder$HostOrder <- gsub(" ", "", wildOrder$HostOrder)
# wildOrder %>% distinct(HostOrder, HostSpecies)

# interaction_traits <- interaction_traits %>% left_join(wildOrder, by = "HostSpecies")
# interaction_traits <- interaction_traits %>% filter(HostOrder == "Primates") # Example: filter for any desired order
# interaction_traits <- interaction_traits %>% select(-HostOrder)



# HOST TRAIT MATRIX

Tr_data <- interaction_traits %>% 
  distinct(HostSpecies, BodyMass, MainGuild, ForStrat, Activity) # Include all traits to be used in final matrix

# Remove NAs and unwanted categories
Tr_data <- Tr_data %>% 
  filter(!is.na(BodyMass)) %>%   # Remove rows with missing BodyMass
  filter(!is.na(MainGuild)) %>%  # Remove rows with missing MainGuild
  filter(!is.na(ForStrat)) %>%   # Remove rows with missing ForStrat
  filter(!is.na(Activity)) %>%   # Remove rows with missing Activity
  filter(MainGuild != "not_found") %>%  # Remove invalid categories
  filter(ForStrat != "not_found") %>% 
  filter(Activity != "not_found") %>% 
  filter(HostSpecies != "Ovis_aries") %>% # Remove domestic sheep
  distinct()

# Arrange rows by species name
Tr_data <- Tr_data %>% arrange(HostSpecies)
M_HxTr <- Tr_data %>% select(-HostSpecies)

# Convert variables to appropriate types
M_HxTr$BodyMass <- as.numeric(M_HxTr$BodyMass)
M_HxTr <- M_HxTr %>%
  mutate(across(c(MainGuild, ForStrat, Activity), as.factor))

# Check structure and print
str(M_HxTr)
print(M_HxTr)

# Visual exploration of relationships among variables
ggpairs(M_HxTr)

# ANOVA tests for BodyMass differences by categorical variables
aov_mainGuild <- aov(BodyMass ~ MainGuild, data = M_HxTr)
summary(aov_mainGuild)

aov_forStrat <- aov(BodyMass ~ ForStrat, data = M_HxTr)
summary(aov_forStrat)

aov_activity <- aov(BodyMass ~ Activity, data = M_HxTr)
summary(aov_activity)

# Correlation among categorical factors
# Convert factors to dummy variables
M_HxTr_dummy <- model.matrix(~ MainGuild + ForStrat + Activity - 1, data = M_HxTr)
cor_matrix <- cor(M_HxTr_dummy)
corrplot(cor_matrix, method = "color", tl.cex = 0.8)

# Mixed factor analysis (PCA for mixed data)
res.famd <- FAMD(M_HxTr, ncp = 5, graph = TRUE)
plot(res.famd, choix = "ind")  # Plot individuals
plot(res.famd, choix = "var")  # Plot variables

# Create a formula with all explanatory variables
# Response variable is artificial - only used to fit model
dummy_response <- rnorm(nrow(M_HxTr))  # Random vector

# Fit linear model with all variables
mod_vif <- lm(dummy_response ~ BodyMass + MainGuild + ForStrat + Activity, data = M_HxTr)

# Calculate VIF values to check collinearity
vif_values <- vif(mod_vif)
print(vif_values)

# None of the variables show problematic collinearity.
# GVIF^(1/(2*Df)) values are all below 2, indicating these variables 
# can be used together in a model without risk of coefficient instability.



# HOST–PARASITE INTERACTION MATRIX

HxP_data <- interaction_traits %>% 
  distinct(HostSpecies, ParSpecies)

# Check structure
str(HxP_data)

# Convert to binary matrix: rows = Hosts, columns = Parasites
M_HxP <- HxP_data %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = ParSpecies, values_from = value, values_fill = 0)

# Optionally set row names and remove species column
M_HxP <- M_HxP %>% column_to_rownames("HostSpecies")

# Print matrix
print(M_HxP)



# PHYLOGENETIC MATRIX

# Load phylogenetic tree in Nexus format
nex_data <- read.nexus(file.path(data.directory, "1000_consensus_mammal_tree_edges_length.nex"))

# Extract unique host species from interaction data
hspecies <- HxP_data %>% distinct(HostSpecies)

# Create regex pattern of species names to retain from phylogenetic matrix
subset_species <- paste(hspecies$HostSpecies, collapse = "|")

# COMPUTE PHYLOGENETIC DISTANCE MATRIX
phylodist_matrix <- cophenetic.phylo(nex_data)

# Convert distance matrix to data frame
phylodist_matrix_df <- as.data.frame(phylodist_matrix)
phylodist_matrix_df$species <- rownames(phylodist_matrix_df)

# Subset distance matrix to species of interest
subset_phylodist_matrix <- phylodist_matrix_df %>% 
  filter(grepl(subset_species, species, ignore.case = TRUE)) %>%
  select(matches(subset_species, ignore.case = TRUE))

# Restore rownames and reorder
subset_phylodist_matrix$species <- rownames(subset_phylodist_matrix)
hspecies <- hspecies %>% arrange(HostSpecies)
subset_phylodist_matrix <- subset_phylodist_matrix %>% arrange(species)
rownames(subset_phylodist_matrix) <- hspecies$HostSpecies

# Reorder columns alphabetically
subset_phylodist_matrix <- subset_phylodist_matrix %>% select(species, everything())
subset_phylodist_matrix <- subset_phylodist_matrix %>% select(sort(names(.)))

# Remove 'species' column from matrix
subset_phylodist_matrix <- subset_phylodist_matrix %>% select(-species)

# Rename columns to match host species order
colnames(subset_phylodist_matrix) <- hspecies$HostSpecies



# FILTER MATRICES TO COMMON HOSTS
M_HxP <- M_HxP[Tr_data$HostSpecies, , drop = FALSE]
subset_phylodist_matrix <- subset_phylodist_matrix[Tr_data$HostSpecies, Tr_data$HostSpecies, drop = FALSE]

# NORMALIZE PHYLOGENETIC DISTANCE (optional)
subset_phylodist_matrix <- subset_phylodist_matrix / max(subset_phylodist_matrix)
subset_phylodist_matrix <- as.matrix(subset_phylodist_matrix)

# COMPUTE JACCARD DISTANCE MATRIX BASED ON SHARED PARASITES
D_jaccard_parasites <- vegdist(M_HxP, method = "jaccard", binary = TRUE)
D_jaccard_parasites <- as.matrix(D_jaccard_parasites)

# COMPUTE GOWER DISTANCE MATRIX BASED ON TRAITS
D_gower_traits <- FD::gowdis(M_HxTr)
D_gower_traits <- as.matrix(D_gower_traits)

############## DATA PREPARATION END ################

########### ANALYSIS ###############

# MANTEL TESTS
mantel_mod <- vegan::mantel(D_jaccard_parasites, D_gower_traits, 
                            method = "pearson", 
                            permutations = 999)

partial_mantel_mod <- mantel.partial(D_jaccard_parasites, D_gower_traits, 
                                     subset_phylodist_matrix, method = "pearson", 
                                     permutations = 999)

# PRINT RESULTS
print(mantel_mod)
print(partial_mantel_mod)

# MRM ANALYSES (Multiple Regression on distance Matrices)
reg_traits <- MRM(as.dist(D_jaccard_parasites) ~ as.dist(D_gower_traits))
reg_phylo <- MRM(as.dist(D_jaccard_parasites) ~ as.dist(subset_phylodist_matrix))
reg_both <- MRM(as.dist(D_jaccard_parasites) ~ as.dist(D_gower_traits) + as.dist(subset_phylodist_matrix))

print(reg_traits)
print(reg_phylo)
print(reg_both)



# PLOT CORRELATIONS

corrdat <- data.frame(
  Parasites = as.vector(as.dist(D_jaccard_parasites)),
  Traits    = as.vector(as.dist(D_gower_traits)),
  Phylo     = as.vector(as.dist(subset_phylodist_matrix))
)

# Remove missing values
corrdat <- na.omit(corrdat)

# Plot Traits vs Parasites
ggplot(corrdat, aes(x = Traits, y = Parasites)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm") + 
  labs(title = "Jaccard distance Parasites vs Gower distance Traits") + 
  theme_minimal()

# Plot Phylo vs Parasites
ggplot(corrdat, aes(x = Phylo, y = Parasites)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm") + 
  labs(title = "Jaccard distance Parasites vs Phylogenetic distance") + 
  theme_minimal()

