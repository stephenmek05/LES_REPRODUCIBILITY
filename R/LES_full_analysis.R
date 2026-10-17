##################################################
rm(list = ls())

# ============================================================
# Reproducibility package root
# ============================================================

project_root <- normalizePath(".")

setwd(project_root)

# Output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("environment", recursive = TRUE, showWarnings = FALSE)
dir.create("checksums", recursive = TRUE, showWarnings = FALSE)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(phytools)
library(V.PhyloMaker2)
library(caper)
library(lme4)
library(lmerTest)
library(vegan)
library(performance)
library(broom.mixed)
library(ggrepel)
##################################################
# Load TRY dataset
##################################################
try_data <- fread("data/raw/49929.txt")
##################################################
# Extract core traits
##################################################

core_traits <- try_data %>%
  filter(
    TraitID %in% c(
      14,    # Leaf N
      42,    # Growth form
      47,    # LDMC
      3106,  # Height
      3117   # SLA
    )
  )

rm(try_data)
##################################################
# Verify traits
##################################################

table(core_traits$TraitName)
##################################################
# Extract growth-form information
##################################################

growth_forms <- core_traits %>%
  dplyr::filter(TraitID == 42) %>%
  dplyr::select(SpeciesName, OrigValueStr) %>%
  dplyr::distinct()
growth_forms <- growth_forms %>%
  dplyr::mutate(
    OrigValueStr = tolower(trimws(OrigValueStr))
  )
growth_forms <- growth_forms %>%
  dplyr::mutate(
    growth_form_standardized = dplyr::case_when(
      
      OrigValueStr %in% c(
        "tree","conifer","palm"
      ) ~ "tree",
      
      OrigValueStr %in% c(
        "shrub","dwarf shrub"
      ) ~ "shrub",
      
      OrigValueStr %in% c(
        "herb","herbs","herbaceous",
        "herbaceous dicot",
        "herbaceous monocot",
        "herbaceous forb",
        "herbaceous dicotyl",
        "herbaceous monocotyl",
        "forb",
        "leguminous forb"
      ) ~ "herb",
      
      OrigValueStr %in% c(
        "grass",
        "graminoid",
        "sedge",
        "grasses&sedges"
      ) ~ "graminoid",
      
      OrigValueStr %in% c(
        "vine",
        "liana",
        "climber",
        "lianas",
        "lianas and climbers"
      ) ~ "climber",
      
      OrigValueStr %in% c(
        "fern",
        "fern ally",
        "fernally",
        "clubmoss",
        "horsetail"
      ) ~ "fern",
      
      TRUE ~ NA_character_
    )
  )
table(
  growth_forms$growth_form_standardized,
  useNA = "ifany"
)
##################################################
# Species-level growth form assignment
##################################################

species_growth_form <- growth_forms %>%
  dplyr::filter(!is.na(growth_form_standardized)) %>%
  dplyr::group_by(SpeciesName) %>%
  dplyr::summarise(
    growth_form =
      names(
        sort(
          table(growth_form_standardized),
          decreasing = TRUE
        )
      )[1],
    .groups = "drop"
  )

species_growth_form <- species_growth_form %>%
  dplyr::filter(growth_form != "fern")
table(species_growth_form$growth_form)
##################################################
# Species-level trait aggregation
##################################################

species_traits <- core_traits %>%
  dplyr::filter(
    TraitID %in% c(
      14,    # Leaf N
      47,    # LDMC
      3106,  # Height
      3117   # SLA
    )
  ) %>%
  dplyr::group_by(
    SpeciesName,
    TraitName
  ) %>%
  dplyr::summarise(
    trait_median = median(
      StdValue,
      na.rm = TRUE
    ),
    .groups = "drop"
  )
species_traits_wide <- species_traits %>%
  tidyr::pivot_wider(
    names_from = TraitName,
    values_from = trait_median
  )
final_trait_matrix <- species_traits_wide %>%
  dplyr::left_join(
    species_growth_form,
    by = "SpeciesName"
  )
##################################################
# Clean species names
##################################################

final_trait_matrix <- final_trait_matrix %>%
  dplyr::filter(
    !is.na(SpeciesName),
    SpeciesName != "",
    !grepl("^\\(", SpeciesName)
  )

final_trait_matrix <- final_trait_matrix %>%
  dplyr::filter(
    grepl(
      "^[A-Za-z\\*][A-Za-z\\.-]+\\s+[a-z]",
      SpeciesName
    )
  )

final_trait_matrix <- final_trait_matrix %>%
  dplyr::mutate(
    SpeciesName = iconv(
      SpeciesName,
      from = "",
      to = "UTF-8",
      sub = ""
    )
  )

final_trait_matrix <- final_trait_matrix %>%
  dplyr::filter(
    grepl(
      "^[A-Za-z\\*][A-Za-z\\.-]+\\s+[a-z]",
      SpeciesName
    )
  )

final_trait_matrix <- final_trait_matrix %>%
  dplyr::mutate(
    SpeciesName = gsub(
      "^\\*",
      "",
      SpeciesName
    )
  )
##################################################
# Rename traits
##################################################

final_trait_matrix <- final_trait_matrix %>%
  dplyr::rename(
    SLA = `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded`,
    
    Leaf_N = `Leaf nitrogen (N) content per leaf dry mass`,
    
    Height = `Plant height vegetative`,
    
    LDMC = `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`
  )
colnames(final_trait_matrix)
##################################################
# Complete LES dataset
##################################################

les_complete <- final_trait_matrix %>%
  dplyr::filter(
    !is.na(SLA),
    !is.na(Leaf_N),
    !is.na(LDMC),
    !is.na(Height),
    !is.na(growth_form)
  )

dim(les_complete)

table(les_complete$growth_form)
##################################################
# Log-transform LES traits
##################################################

les_log <- les_complete %>%
  dplyr::filter(
    SLA > 0,
    Leaf_N > 0,
    Height > 0,
    LDMC > 0
  ) %>%
  dplyr::mutate(
    log_SLA    = log10(SLA),
    log_Leaf_N = log10(Leaf_N),
    log_Height = log10(Height),
    log_LDMC   = log10(LDMC)
  )
dim(les_log)

summary(
  les_log %>%
    dplyr::select(
      log_SLA,
      log_Leaf_N,
      log_Height,
      log_LDMC
    )
)
# One species removed due to non-positive trait value(s),
# which cannot be log-transformed.
##################################################
# Principal Components Analysis
##################################################

les_pca <- prcomp(
  les_log %>%
    dplyr::select(
      log_SLA,
      log_Leaf_N,
      log_Height,
      log_LDMC
    ),
  scale. = TRUE
)

summary(les_pca)

les_pca$rotation
loadings_df <- data.frame(
  Trait = rownames(les_pca$rotation),
  PC1 = les_pca$rotation[,1],
  PC2 = les_pca$rotation[,2]
)
les_scores <- cbind(
  les_log,
  as.data.frame(les_pca$x)
)
summary(les_pca)
##################################################
# PCA scores
##################################################

les_scores <- cbind(
  les_log,
  as.data.frame(les_pca$x)
)
##################################################
# Growth-form overlap
# Bhattacharyya distances
##################################################

pc_space <- les_scores %>%
  dplyr::select(
    PC1,
    PC2,
    growth_form
  )

groups <- unique(pc_space$growth_form)

bhat <- matrix(
  NA,
  length(groups),
  length(groups),
  dimnames = list(groups, groups)
)

for(i in 1:length(groups)){
  for(j in i:length(groups)){
    
    g1 <- subset(
      pc_space,
      growth_form == groups[i],
      select = c(PC1, PC2)
    )
    
    g2 <- subset(
      pc_space,
      growth_form == groups[j],
      select = c(PC1, PC2)
    )
    
    m1 <- colMeans(g1)
    m2 <- colMeans(g2)
    
    S1 <- cov(g1)
    S2 <- cov(g2)
    
    S <- (S1 + S2)/2
    
    term1 <-
      0.125 *
      t(m1-m2) %*%
      solve(S) %*%
      (m1-m2)
    
    term2 <-
      0.5 *
      log(
        det(S) /
          sqrt(det(S1)*det(S2))
      )
    
    bhat[i,j] <- term1 + term2
    bhat[j,i] <- bhat[i,j]
  }
}

round(bhat, 3)

write.csv(
  round(bhat, 3),
  "data/processed/bhattacharyya_distances.csv"
)
##################################################
# Growth-form differences along LES axis
##################################################

anova_pc1 <- aov(
  PC1 ~ growth_form,
  data = les_scores
)

summary(anova_pc1)

TukeyHSD(anova_pc1)
eta_squared <- 
  summary(anova_pc1)[[1]][["Sum Sq"]][1] /
  sum(
    summary(anova_pc1)[[1]][["Sum Sq"]]
  )

eta_squared
summary(anova_pc1)
eta_squared
##################################################
# Species list for phylogeny
##################################################

species_list <- les_scores %>%
  dplyr::select(SpeciesName) %>%
  dplyr::distinct()

species_list <- species_list %>%
  tidyr::separate(
    SpeciesName,
    into = c("genus", "species"),
    sep = " ",
    remove = FALSE
  )
head(species_list)
nrow(species_list)
library(V.PhyloMaker2)

data(GBOTB.extended.WP)
data(nodes.info.1.WP)
phylo <- phylo.maker(
  sp.list = species_list,
  tree = GBOTB.extended.WP,
  nodes = nodes.info.1.WP
)
phy <- phylo$scenario.3
length(phy$tip.label)
les_scores <- les_scores %>%
  dplyr::mutate(
    phylo_name = gsub(" ", "_", SpeciesName)
  )

les_phylo <- les_scores %>%
  dplyr::filter(
    phylo_name %in% phy$tip.label
  )
dim(les_phylo)
##################################################
# Prune phylogeny
##################################################

phy_pruned <- drop.tip(
  phy,
  setdiff(
    phy$tip.label,
    les_phylo$phylo_name
  )
)

les_phylo <- les_phylo %>%
  dplyr::slice(
    match(
      phy_pruned$tip.label,
      phylo_name
    )
  )
all(
  phy_pruned$tip.label ==
    les_phylo$phylo_name
)
##################################################
# Phylogenetic signal
##################################################

library(phytools)

pc1_vector <- les_phylo$PC1
names(pc1_vector) <- les_phylo$phylo_name
blomberg_k <- phylosig(
  phy_pruned,
  pc1_vector,
  method = "K",
  test = TRUE
)

pagel_lambda <- phylosig(
  phy_pruned,
  pc1_vector,
  method = "lambda",
  test = TRUE
)
blomberg_k

pagel_lambda
##################################################
# Phylogenetic ANOVA
##################################################

growth_vector <- les_phylo$growth_form
names(growth_vector) <- les_phylo$phylo_name
set.seed(123)
phylo_anova <- phylANOVA(
  tree = phy_pruned,
  x = growth_vector,
  y = pc1_vector,
  nsim = 1000
)

phylo_anova
##################################################
# Mixed-effects model
##################################################

les_phylo <- les_phylo %>%
  tidyr::separate(
    phylo_name,
    into = c("genus", "species"),
    sep = "_",
    remove = FALSE
  )

model_genus <- lmer(
  PC1 ~ growth_form +
    (1 | genus),
  data = les_phylo
)

summary(model_genus)
library(performance)

r2(model_genus)
library(broom.mixed)

coef_df <- tidy(
  model_genus,
  effects = "fixed",
  conf.int = TRUE
)

coef_df
summary(model_genus)
coef_df
##################################################
# PGLS
##################################################

library(caper)

phy_pruned$node.label <- NULL

comp_data <- comparative.data(
  phy = phy_pruned,
  data = les_phylo,
  names.col = phylo_name,
  vcv = TRUE,
  warn.dropped = TRUE
)

pgls_model <- pgls(
  PC1 ~ growth_form,
  data = comp_data,
  lambda = "ML"
)

summary(pgls_model)
##################################################
# Balanced subsampling robustness
##################################################

robustness_results <- data.frame(
  iteration = integer(),
  F_value = numeric(),
  p_value = numeric(),
  eta_squared = numeric()
)

set.seed(123)

for(i in 1:100){
  
  temp_data <- les_phylo %>%
    dplyr::group_by(growth_form) %>%
    dplyr::sample_n(99) %>%
    dplyr::ungroup()
  
  temp_pca <- prcomp(
    temp_data %>%
      dplyr::select(
        log_SLA,
        log_Leaf_N,
        log_Height,
        log_LDMC
      ),
    scale. = TRUE
  )
  
  temp_scores <- cbind(
    temp_data,
    as.data.frame(temp_pca$x)
  )
  
  temp_aov <- aov(
    PC1 ~ growth_form,
    data = temp_scores
  )
  
  aov_table <- summary(temp_aov)[[1]]
  
  robustness_results <- rbind(
    robustness_results,
    data.frame(
      iteration = i,
      F_value = aov_table$`F value`[1],
      p_value = aov_table$`Pr(>F)`[1],
      eta_squared =
        aov_table$`Sum Sq`[1] /
        sum(aov_table$`Sum Sq`)
    )
  )
}
summary(robustness_results)
##################################################
# PERMANOVA
##################################################

trait_matrix <- les_phylo %>%
  dplyr::select(
    log_SLA,
    log_Leaf_N,
    log_Height,
    log_LDMC
  )

permanova_result <- adonis2(
  trait_matrix ~ growth_form,
  data = les_phylo,
  method = "euclidean",
  permutations = 999
)

permanova_result
##################################################
# Save analytical outputs
##################################################

saveRDS(
  les_pca,
  "data/processed/les_pca_model.rds"
)

saveRDS(
  phy_pruned,
  "data/processed/phylogeny_pruned.rds"
)

saveRDS(
  phylo_anova,
  "data/processed/phylogenetic_anova.rds"
)

saveRDS(
  pgls_model,
  "data/processed/pgls_model.rds"
)

saveRDS(
  model_genus,
  "data/processed/mixed_model.rds"
)

write.csv(
  robustness_results,
  "data/processed/robustness_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(permanova_result),
  "data/processed/permanova_results.csv"
)

writeLines(
  capture.output(sessionInfo()),
  "environment/sessionInfo.txt"
)
##################################################
# Sensitivity Analysis 1
# LES without Height
##################################################

les_noheight <- les_complete %>%
  dplyr::select(
    SpeciesName,
    SLA,
    Leaf_N,
    LDMC,
    growth_form
  )

les_noheight_log <- les_noheight %>%
  dplyr::mutate(
    log_SLA = log10(SLA),
    log_Leaf_N = log10(Leaf_N),
    log_LDMC = log10(LDMC)
  )
pca_noheight <- prcomp(
  les_noheight_log %>%
    dplyr::select(
      log_SLA,
      log_Leaf_N,
      log_LDMC
    ),
  center = TRUE,
  scale. = TRUE
)

summary(pca_noheight)
scores_noheight <- cbind(
  les_noheight,
  as.data.frame(pca_noheight$x)
)
anova_noheight <- aov(
  PC1 ~ growth_form,
  data = scores_noheight
)

summary(anova_noheight)
eta_noheight <-
  summary(anova_noheight)[[1]][["Sum Sq"]][1] /
  sum(summary(anova_noheight)[[1]][["Sum Sq"]])

eta_noheight
library(vegan)

trait_matrix_noheight <- les_noheight_log %>%
  dplyr::select(
    log_SLA,
    log_Leaf_N,
    log_LDMC
  )

permanova_noheight <- adonis2(
  trait_matrix_noheight ~ growth_form,
  data = les_noheight_log,
  method = "euclidean",
  permutations = 999
)

permanova_noheight
summary(pca_noheight)

summary(anova_noheight)

eta_noheight

permanova_noheight
set.seed(123)

pc1_variance <- numeric(1000)

for(i in 1:1000){
  
  tmp <- les_log %>%
    dplyr::sample_frac(
      replace = TRUE
    )
  
  pca_tmp <- prcomp(
    tmp %>%
      dplyr::select(
        log_SLA,
        log_Leaf_N,
        log_Height,
        log_LDMC
      ),
    center = TRUE,
    scale. = TRUE
  )
  
  pc1_variance[i] <-
    summary(pca_tmp)$importance[2,1]
}
summary(pc1_variance)

sd(pc1_variance)

quantile(
  pc1_variance,
  c(0.025, 0.975)
)
##################################################
# PERMDISP (matched to PERMANOVA dataset)
##################################################

library(vegan)

dist_mat_phylo <- dist(
  trait_matrix
)

disp_phylo <- betadisper(
  dist_mat_phylo,
  les_phylo$growth_form
)

anova(disp_phylo)

permutest(
  disp_phylo,
  permutations = 999
)

saveRDS(
  disp_phylo,
  "data/processed/permdisp_results.rds"
)
##################################################
# Additional Diagnostics and Robustness Analyses
##################################################

##################################################
# Growth-form assignment conflicts
##################################################

growth_form_conflicts <- growth_forms %>%
  dplyr::filter(!is.na(growth_form_standardized)) %>%
  dplyr::group_by(SpeciesName) %>%
  dplyr::summarise(
    n_forms = n_distinct(growth_form_standardized),
    forms = paste(
      sort(unique(growth_form_standardized)),
      collapse = ", "
    ),
    .groups = "drop"
  )

conflict_summary <- growth_form_conflicts %>%
  dplyr::count(
    forms,
    sort = TRUE
  )

conflict_summary

##################################################
# Filtering bias assessment
##################################################

before_counts <- species_growth_form %>%
  dplyr::count(growth_form) %>%
  dplyr::rename(
    n_before = n
  )

after_counts <- les_complete %>%
  dplyr::count(growth_form) %>%
  dplyr::rename(
    n_after = n
  )

filtering_bias <- before_counts %>%
  dplyr::left_join(
    after_counts,
    by = "growth_form"
  ) %>%
  dplyr::mutate(
    prop_before =
      n_before / sum(n_before),
    prop_after =
      n_after / sum(n_after),
    retention_rate =
      n_after / n_before
  )

filtering_bias

##################################################
# Chi-square test of growth-form representation
##################################################

chisq_filtering <- chisq.test(
  rbind(
    before_counts$n_before,
    after_counts$n_after
  )
)

chisq_filtering

##################################################
# Retention percentages by growth form
##################################################

filtering_bias_summary <- filtering_bias %>%
  dplyr::mutate(
    retention_pct =
      round(100 * retention_rate, 2),
    relative_to_mean =
      retention_rate /
      mean(retention_rate)
  )

filtering_bias_summary

##################################################
# Ambiguous species in final dataset
##################################################

ambiguous_final <- les_complete %>%
  dplyr::filter(
    SpeciesName %in%
      growth_form_conflicts$SpeciesName[
        growth_form_conflicts$n_forms > 1
      ]
  )

n_ambiguous <- nrow(ambiguous_final)

prop_ambiguous <-
  nrow(ambiguous_final) /
  nrow(les_complete)

n_ambiguous
prop_ambiguous

##################################################
# Sensitivity analysis:
# Unambiguous growth-form assignments only
##################################################

strict_species <- growth_form_conflicts %>%
  dplyr::filter(
    n_forms == 1
  ) %>%
  dplyr::pull(
    SpeciesName
  )

les_complete_strict <- les_complete %>%
  dplyr::filter(
    SpeciesName %in% strict_species
  )

nrow(les_complete_strict)

strict_log <- les_complete_strict %>%
  dplyr::filter(
    SLA > 0,
    Leaf_N > 0,
    Height > 0,
    LDMC > 0
  ) %>%
  dplyr::mutate(
    log_SLA = log10(SLA),
    log_Leaf_N = log10(Leaf_N),
    log_Height = log10(Height),
    log_LDMC = log10(LDMC)
  )

strict_pca <- prcomp(
  strict_log %>%
    dplyr::select(
      log_SLA,
      log_Leaf_N,
      log_Height,
      log_LDMC
    ),
  scale. = TRUE
)

strict_scores <- cbind(
  strict_log,
  as.data.frame(strict_pca$x)
)

strict_aov <- aov(
  PC1 ~ growth_form,
  data = strict_scores
)

summary(strict_aov)

strict_eta <-
  summary(strict_aov)[[1]][["Sum Sq"]][1] /
  sum(
    summary(strict_aov)[[1]][["Sum Sq"]]
  )

strict_eta

##################################################
# Save additional outputs
##################################################

write.csv(
  conflict_summary,
  "data/processed/growth_form_conflict_summary.csv",
  row.names = FALSE
)

write.csv(
  filtering_bias,
  "data/processed/filtering_bias_summary.csv",
  row.names = FALSE
)

write.csv(
  ambiguous_final,
  "data/processed/ambiguous_species_in_final_dataset.csv",
  row.names = FALSE
)
writeLines(
  capture.output(
    list(
      chi_square_filtering = chisq_filtering,
      ambiguous_species_n = n_ambiguous,
      ambiguous_species_prop = prop_ambiguous,
      strict_dataset_n = nrow(les_complete_strict),
      original_eta_squared = eta_squared,
      strict_eta_squared = strict_eta
    )
  ),
  "data/processed/additional_robustness_summary.txt"
)

########################################
#### TO GENERATE FIGURES ###############

##### Figure 1 #####
library(ggplot2)

cols <- c(
  climber   = "#F8766D",
  graminoid = "#B79F00",
  herb      = "#00BA38",
  shrub     = "#00B0F6",
  tree      = "#E76BF3"
)

legend_df <- data.frame(
  PC1 = 0,
  PC2 = 0,
  growth_form = names(cols)
)

figure1 <- ggplot(
  les_scores,
  aes(
    x = PC1,
    y = PC2,
    fill = growth_form
  )
) +
  
  stat_density_2d(
    aes(alpha = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    show.legend = FALSE
  ) +
  
  scale_alpha(
    range = c(0.05, 0.45),
    guide = "none"
  ) +
  
  geom_point(
    data = legend_df,
    aes(
      x = PC1,
      y = PC2,
      colour = growth_form
    ),
    alpha = 0,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = cols
  ) +
  
  scale_colour_manual(
    values = cols,
    name = "Growth form"
  ) +
  
  labs(
    title =
      "Growth-form occupancy across multidimensional LES space",
    x =
      paste0(
        "LES axis (PC1; ",
        round(
          100 *
            summary(les_pca)$importance[2,1],
          1
        ),
        "% variance)"
      ),
    y =
      paste0(
        "Secondary trait axis (PC2; ",
        round(
          100 *
            summary(les_pca)$importance[2,2],
          1
        ),
        "% variance)"
      )
  ) +
  
  theme_bw(base_size = 18) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 24
    ),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        alpha = 1,
        size = 5
      )
    )
  )

figure1
ggsave(
  filename = "Figure1_LES_space.pdf",
  plot = figure1,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "Figure1_LES_space.png",
  plot = figure1,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
##################################################
# Figure 2
# Distribution of LES positioning among growth forms
##################################################

library(ggplot2)
library(dplyr)
library(ggsignif)

##################################################
# Sample sizes from PCA dataset
##################################################

n_labels <- les_scores %>%
  dplyr::count(growth_form)

##################################################
# Color palette
##################################################

cols <- c(
  climber   = "#F8766D",
  graminoid = "#B79F00",
  herb      = "#00BA38",
  shrub     = "#00B0F6",
  tree      = "#E76BF3"
)

##################################################
# Figure
##################################################

figure2 <- ggplot(
  les_scores,
  aes(
    x = growth_form,
    y = PC1,
    fill = growth_form
  )
) +
  
  geom_violin(
    alpha = 0.70,
    colour = "black",
    linewidth = 0.8,
    trim = FALSE
  ) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.08,
    size = 1.2,
    colour = "black"
  ) +
  
  geom_boxplot(
    width = 0.12,
    fill = "white",
    alpha = 0.8,
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  
  geom_signif(
    comparisons = list(
      c("climber", "shrub"),
      c("climber", "tree"),
      c("herb", "shrub"),
      c("herb", "tree")
    ),
    annotations = c(
      "****",
      "****",
      "****",
      "****"
    ),
    y_position = c(
      5.2,
      7.8,
      7.0,
      6.2
    ),
    tip_length = 0.01,
    textsize = 6
  ) +
  
  geom_text(
    data = n_labels,
    aes(
      x = growth_form,
      y = -5.4,
      label = paste0("n = ", n)
    ),
    inherit.aes = FALSE,
    size = 5
  ) +
  
  scale_fill_manual(
    values = cols
  ) +
  
  coord_cartesian(
    ylim = c(-5.8, 8.5)
  ) +
  
  labs(
    title = "Growth forms occupy distinct positions along the LES axis",
    x = "Growth form",
    y = "LES axis (PC1)"
  ) +
  
  theme_classic(base_size = 20) +
  
  theme(
    legend.position = "none",
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 22
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      colour = "black"
    )
  )

figure2

##################################################
# Save Figure 2
##################################################

ggsave(
  filename = "Figure2_PC1_growthform_violin.png",
  plot = figure2,
  path = "figures",
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = "Figure2_PC1_growthform_violin.pdf",
  plot = figure2,
  path = "figures",
  width = 12,
  height = 8,
  units = "in",
  bg = "white",
  limitsize = FALSE
)
##################################################
# Figure 3
# Bhattacharyya distance heatmap
##################################################

library(ggplot2)
library(reshape2)

bhat_df <- melt(round(bhat, 3))

colnames(bhat_df) <- c(
  "GrowthForm1",
  "GrowthForm2",
  "Distance"
)

figure3 <- ggplot(
  bhat_df,
  aes(
    GrowthForm1,
    GrowthForm2,
    fill = Distance
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.3f", Distance)
    ),
    size = 6,
    fontface = "bold"
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "firebrick"
  ) +
  
  labs(
    title =
      "Pairwise growth-form separation in multidimensional LES space",
    x = "",
    y = "",
    fill = "Bhattacharyya\nDistance"
  ) +
  
  coord_fixed() +
  
  theme_classic(base_size = 18) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      face = "bold"
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    )
  )

figure3
##################################################
# Save Figure 3 
##################################################

ggsave(
  filename = "Figure3_Bhattacharyya_heatmap.pdf",
  plot = figure3,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = "Figure3_Bhattacharyya_heatmap.png",
  plot = figure3,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
##################################################
# Figure 4
# Robustness of growth-form differentiation
##################################################

library(ggplot2)
library(dplyr)
library(patchwork)

##################################################
# Panel A
# Sensitivity analyses
##################################################

robust_summary <- data.frame(
  Analysis = factor(
    c(
      "No height",
      "Full dataset",
      "Strict assignment"
    ),
    levels = c(
      "No height",
      "Full dataset",
      "Strict assignment"
    )
  ),
  Eta2 = c(
    eta_noheight,
    eta_squared,
    strict_eta
  )
)

panelA <- ggplot(
  robust_summary,
  aes(
    x = Analysis,
    y = Eta2
  )
) +
  
  geom_col(
    fill = "steelblue"
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.3f", Eta2)
    ),
    vjust = -0.4,
    size = 6
  ) +
  
  ylim(0, 0.60) +
  
  labs(
    title = "A. Sensitivity analyses",
    x = "",
    y = expression(eta^2)
  ) +
  
  theme_classic(base_size = 18)

##################################################
# Panel B
# Balanced subsampling robustness
##################################################

panelB <- ggplot(
  robustness_results,
  aes(x = eta_squared)
) +
  
  geom_histogram(
    bins = 15,
    fill = "tomato",
    colour = "black"
  ) +
  
  geom_vline(
    xintercept = eta_squared,
    linetype = 2,
    linewidth = 1
  ) +
  
  labs(
    title = "B. Balanced subsampling robustness",
    x = expression(eta^2),
    y = "Frequency"
  ) +
  
  theme_classic(base_size = 18)

##################################################
# Combine
##################################################

figure4 <- panelA + panelB

figure4

ggsave(
  filename = "Figure4_Robustness.pdf",
  plot = figure4,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = "Figure4_Robustness.png",
  plot = figure4,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)


##################################################
# Figure 5
# Phylogenetic persistence of growth-form effects
##################################################

library(ggplot2)
library(dplyr)
library(patchwork)

##################################################
# PGLS coefficients
##################################################

pgls_plot <- data.frame(
  term = c(
    "Graminoid",
    "Herb",
    "Shrub",
    "Tree"
  ),
  
  estimate = c(
    0.5479823,
    0.5096245,
    -0.5834247,
    -0.6386739
  ),
  
  se = c(
    0.2499252,
    0.1178218,
    0.1150259,
    0.1170897
  )
)

pgls_plot <- pgls_plot %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se
  )

##################################################
# Mixed model coefficients
##################################################

mixed_plot <- coef_df %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = c(
      "Graminoid",
      "Herb",
      "Shrub",
      "Tree"
    )
  )

##################################################
# Panel A
##################################################

panelA <- ggplot(
  pgls_plot,
  aes(
    x = estimate,
    y = term
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  
  geom_errorbar(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.15,
    linewidth = 1
  ) +
  
  geom_point(
    size = 4,
    colour = "firebrick"
  ) +
  
  labs(
    title = "A. PGLS",
    x = "Coefficient estimate",
    y = ""
  ) +
  
  theme_classic(base_size = 18)

##################################################
# Panel B
##################################################

panelB <- ggplot(
  mixed_plot,
  aes(
    x = estimate,
    y = term
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey40"
  ) +
  
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.15,
    linewidth = 1
  ) +
  
  geom_point(
    size = 4,
    colour = "steelblue"
  ) +
  
  labs(
    title = "B. Mixed-effects model",
    x = "Coefficient estimate",
    y = ""
  ) +
  
  theme_classic(base_size = 18)

##################################################
# Combine
##################################################

figure5 <- panelA + panelB

figure5
ggsave(
  filename = "Figure5_Phylogenetic_persistence.pdf",
  plot = figure5,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = "Figure5_Phylogenetic_persistence.png",
  plot = figure5,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
##################################################
# Supplementary Figure S1
# PERMDISP
##################################################

library(ggplot2)

disp_df <- data.frame(
  Distance = disp_phylo$distances,
  GrowthForm = disp_phylo$group
)

disp_df$GrowthForm <- factor(
  disp_df$GrowthForm,
  levels = c(
    "tree",
    "shrub",
    "herb",
    "graminoid",
    "climber"
  )
)

gf_cols <- c(
  tree      = "#E69F00",
  shrub     = "#56B4E9",
  herb      = "#009E73",
  graminoid = "#F0E442",
  climber   = "#CC79A7"
)

figureS1 <- ggplot(
  disp_df,
  aes(
    x = GrowthForm,
    y = Distance,
    fill = GrowthForm
  )
) +
  
  geom_violin(
    alpha = 0.8,
    colour = "black"
  ) +
  
  geom_boxplot(
    width = 0.15,
    fill = "white",
    outlier.shape = NA
  ) +
  
  scale_fill_manual(
    values = gf_cols
  ) +
  
  labs(
    x = "Growth form",
    y = "Distance to group centroid"
  ) +
  
  theme_classic(base_size = 18) +
  
  theme(
    legend.position = "none"
  )

figureS1
ggsave(
  filename = "FigureS1_PERMDISP.pdf",
  plot = figureS1,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = "FigureS1_PERMDISP.png",
  plot = figureS1,
  path = "figures",
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
##################################################
# Table 1
# Growth-form representation after filtering
##################################################

table1 <- filtering_bias_summary %>%
  dplyr::select(
    growth_form,
    n_before,
    n_after,
    retention_pct
  ) %>%
  dplyr::mutate(
    growth_form = tools::toTitleCase(growth_form),
    
    retention_pct = round(
      retention_pct,
      2
    )
  ) %>%
  dplyr::rename(
    `Growth form` = growth_form,
    `Before filtering` = n_before,
    `After filtering` = n_after,
    `Retention (%)` = retention_pct
  )

table1

write.csv(
  table1,
  "tables/Table1_GrowthForm_Filtering.csv",
  row.names = FALSE
)
##################################################
# Table 2
# Summary of phylogenetic analyses
##################################################

table2 <- data.frame(
  Analysis = c(
    "Blomberg's K",
    "Pagel's lambda",
    "Phylogenetic ANOVA",
    "PGLS",
    "PGLS",
    "Mixed model",
    "Mixed model"
  ),
  
  Statistic = c(
    "K",
    "Lambda",
    "F",
    "F",
    "Lambda (ML)",
    "Marginal R2",
    "Conditional R2"
  ),
  
  Value = c(
    0.012,
    0.918,
    524.25,
    58.54,
    0.870,
    0.306,
    0.680
  ),
  
  P_value = c(
    "0.001",
    "<0.001",
    "0.004",
    "<0.001",
    "-",
    "-",
    "-"
  )
)

table2
write.csv(
  table2,
  "tables/Table2_Phylogenetic_Analyses.csv",
  row.names = FALSE
)
