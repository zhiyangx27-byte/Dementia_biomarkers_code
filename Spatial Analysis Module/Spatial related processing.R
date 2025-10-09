pacman::p_load(
  dplyr, ggplot2, tidyr, purrr, stringr,  
  lme4, lmerTest,                          
  brms, rstan,                             
  spdep, sf, spatstat,                     
  patchwork                               
)

merfish_data <- read.csv("/Users/a1234/Desktop/merfish_for_r_analysis.csv", stringsAsFactors = FALSE)

cat("Original columns:\n")
print(names(merfish_data))

merfish_data <- merfish_data %>%
  rename(
    Donor_ID = 1,         
    Cognitive_Status = 2, 
    Age_at_Death = 3,     
    Sex = 4,              
    Subclass = 5,         
    Layer_annotation = 6, 
    Cell_ID = 7,          
    x = 8,                
    y = 9                
  ) %>%
  mutate(
    Donor_ID = as.factor(Donor_ID),
    Cognitive_Status = as.factor(Cognitive_Status),
    Sex = as.factor(Sex),
    Subclass = as.factor(Subclass),
    Layer_annotation = as.factor(Layer_annotation)
  ) %>%
  filter(Subclass %in% c("L2/3 IT", "L4 IT", "L5 IT", "Oligodendrocyte", "Vip"))

cat("\nDate dimension after processing:", dim(merfish_data), "\n")
cat("Number of donors:", length(unique(merfish_data$Donor_ID)), "\n")
cat("Distribution of subclasses of cells:\n")
print(table(merfish_data$Subclass))
cat("Distribution of organizational layers:\n")
print(table(merfish_data$Layer_annotation))

target_subclasses <- c("L2/3 IT", "L4 IT", "L5 IT", "Oligodendrocyte", "Vip")
all_layers <- levels(merfish_data$Layer_annotation)

donor_level_data <- merfish_data %>%
  group_by(Donor_ID, Subclass, Layer_annotation, Cognitive_Status, Age_at_Death, Sex) %>%
  summarise(
    total_cells = n(),  
    .groups = "drop"
  ) %>%
  group_by(Donor_ID, Subclass) %>%
  mutate(
    layer_ratio = total_cells / sum(total_cells) * 100 
  ) %>%
  ungroup()

print(donor_level_data)


corrected_data <- donor_level_data %>%
  mutate(
    Age_at_Death = case_when(
      Age_at_Death == "90+" ~ NA_character_,
      TRUE ~ Age_at_Death
    ),
    Age_at_Death = as.numeric(Age_at_Death)
  ) %>%
  group_by(Donor_ID, Subclass) %>%
  mutate(
    N_total = sum(total_cells),
    layer_prob = layer_ratio / 100, 
    layer_prob_corrected = case_when(
      N_total <= 2 ~ layer_prob, 
      TRUE ~ (layer_prob * (N_total - 2) + 1) / N_total
    )
  ) %>%
  ungroup() %>%
  filter(layer_prob_corrected > 0 & layer_prob_corrected < 1)

cat("\nThe Age_at_Death column has been corrected and converted to numeric type。\n")
print(unique(corrected_data$Age_at_Death))

cat("The data dimensions prepared by Beta GLMM:", dim(corrected_data), "\n")

corrected_data$Age_at_Death <- as.numeric(corrected_data$Age_at_Death)

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE) 


priors <- set_prior("normal(0, 2)", class = "b") # 'b' class for all population-level effects (fixed effects)

beta_model <- brms::brm(
  layer_prob_corrected ~ Cognitive_Status + Subclass + Layer_annotation + 
    Cognitive_Status:Subclass + 
    Cognitive_Status:Layer_annotation + 
    Subclass:Layer_annotation + 
    Age_at_Death + Sex + 
    (1 | Donor_ID),
  data = corrected_data,
  family = Beta(),
  prior = priors,        
  chains = 4, 
  iter = 4000, 
  warmup = 2000, 
  seed = 42,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

print(beta_model)

target_subclasses <- c("L2/3 IT", "L4 IT", "L5 IT", "Oligodendrocyte", "Vip")

run_moran_analysis <- function(subclass_name, all_merfish_data) {
  moran_results_list <- all_merfish_data %>%
    group_by(Donor_ID) %>%
    nest() %>%
    ungroup() %>%
    mutate(
      Cognitive_Status = map_chr(data, ~ as.character(unique(.$Cognitive_Status))), 
      
      moran_analysis = map(data, ~ {
        df_donor <- .
        coords_matrix <- as.matrix(df_donor[, c("x", "y")])
        
        if (nrow(df_donor) < 10) return(list(Moran_I = NA, P_Value = NA))
        
        nb_delaunay <- try(tri2nb(coords_matrix), silent = TRUE)
        if (inherits(nb_delaunay, "try-error")) return(list(Moran_I = NA, P_Value = NA))
        listw <- nb2listw(nb_delaunay, style = "B")
        
        is_target <- as.numeric(df_donor$Subclass == subclass_name)
        
        if (sum(is_target) < 5) return(list(Moran_I = NA, P_Value = NA))
        
        m_test <- try(moran.test(is_target, listw), silent = TRUE)
        
        if (inherits(m_test, "try-error")) {
          return(list(Moran_I = NA, P_Value = NA))
        } else {
          return(list(Moran_I = m_test$estimate[1], P_Value_Local = m_test$p.value))
        }
      })
    ) %>%
    unnest_wider(moran_analysis) %>%
    filter(!is.na(Moran_I))
  
  t_test_data <- moran_results_list %>%
    filter(Cognitive_Status %in% c("Dementia", "No dementia")) 
  
  if (nrow(t_test_data) > 3) { 
    t_result <- t.test(Moran_I ~ Cognitive_Status, data = t_test_data)
    
    return(data.frame(
      Subclass = subclass_name,
      N_donors = nrow(t_test_data),
      Dementia_Moran_I_Mean = mean(t_test_data$Moran_I[t_test_data$Cognitive_Status == "Dementia"]),
      NoDementia_Moran_I_Mean = mean(t_test_data$Moran_I[t_test_data$Cognitive_Status == "No dementia"]),
      T_Statistic = t_result$statistic,
      P_Value_Group = t_result$p.value
    ))
  } else {
    return(data.frame(
      Subclass = subclass_name,
      N_donors = nrow(t_test_data),
      Dementia_Moran_I_Mean = NA,
      NoDementia_Moran_I_Mean = NA,
      T_Statistic = NA,
      P_Value_Group = NA
    ))
  }
}

all_results <- map_dfr(target_subclasses, ~ run_moran_analysis(., merfish_data))

final_moran_comparison <- all_results %>%
  mutate(
    FDR_Group = p.adjust(P_Value_Group, method = "fdr")
  ) %>%
  arrange(FDR_Group)

cat("\n--- Final comparison of Moran's I differences between groups (FDR corrected) ---\n")
print(final_moran_comparison)
