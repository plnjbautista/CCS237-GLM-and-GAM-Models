# ==============================================================================
# GLM vs GAM Analysis - STAR98 Dataset (R Implementation)
# Author: Pauline Joy O. Bautista
# ==============================================================================

# Load libraries
library(tidyverse)
library(car)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(knitr)

# Verify they loaded
cat("\nPackages loaded successfully:\n")
cat("  tidyverse:", as.character(packageVersion("tidyverse")), "\n")
cat("  car:", as.character(packageVersion("car")), "\n")
cat("  mgcv:", as.character(packageVersion("mgcv")), "\n")
cat("  ggplot2:", as.character(packageVersion("ggplot2")), "\n")
cat("  gridExtra:", as.character(packageVersion("gridExtra")), "\n")
cat("  knitr:", as.character(packageVersion("knitr")), "\n")

# ==============================================================================
# 1. LOAD DATASET
# ==============================================================================

url <- "https://raw.githubusercontent.com/statsmodels/statsmodels/main/statsmodels/datasets/star98/star98.csv"
df <- read.csv(url)
names(df) <- toupper(names(df))

# Rename MATHTOT to NABOVE
df <- df %>% rename(NABOVE = MATHTOT)

cat("Dataset loaded successfully.\n")
cat("First few rows:\n")
print(head(df))

# Check for non-positive values in NABOVE
cat("\n### Checking NABOVE for Gamma compatibility:\n")
cat("Min NABOVE value:", min(df$NABOVE, na.rm = TRUE), "\n")
cat("Number of non-positive values:", sum(df$NABOVE <= 0, na.rm = TRUE), "\n")
cat("Number of NA values:", sum(is.na(df$NABOVE)), "\n")

# ==============================================================================
# 2. DESCRIPTIVE STATISTICS
# ==============================================================================

cat("\n### Descriptive Statistics:\n")
print(summary(df))

cat("\n### Missing Values:\n")
print(colSums(is.na(df)))

cat("\n### Data Types:\n")
str(df)

# ==============================================================================
# 3. VISUALIZATIONS - HISTOGRAMS
# ==============================================================================

key_variables <- c('NABOVE', 'NBELOW', 'LOWINC', 'PERASIAN', 'PERBLACK', 
                   'PERHISP', 'AVYRSEXP', 'AVSALK', 'PERSPENK')

# Filter only existing variables
key_variables <- key_variables[key_variables %in% names(df)]

cat("\n### Creating Histograms of Key Numerical Variables\n")

# Histograms
hist_data <- df %>% 
  select(all_of(key_variables)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

p_hist <- ggplot(hist_data, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Histograms of Key Numerical Variables",
       x = "Value", y = "Frequency")

print(p_hist)

# ==============================================================================
# 4. VISUALIZATIONS - BOXPLOTS
# ==============================================================================

cat("\n### Creating Box Plots of Key Numerical Variables\n")

p_box <- ggplot(hist_data, aes(y = value)) +
  geom_boxplot(fill = "orange", alpha = 0.7) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Box Plots of Key Numerical Variables",
       y = "Value")

print(p_box)

cat("Histograms and Box Plots for key numerical variables have been generated.\n")

# ==============================================================================
# 5. CORRELATION MATRIX AND HEATMAP
# ==============================================================================

cat("\n### Calculating Correlation Matrix\n")

# Identify independent variables (drop NABOVE and NBELOW if exists)
cols_to_drop <- c("NABOVE", "NBELOW")
cols_to_drop <- cols_to_drop[cols_to_drop %in% names(df)]
independent_vars <- df %>% select(-all_of(cols_to_drop))

# Calculate correlation matrix
correlation_matrix <- cor(independent_vars, use = "complete.obs")

cat("Correlation matrix calculated. First 5x5:\n")
print(round(correlation_matrix[1:min(5, nrow(correlation_matrix)), 
                               1:min(5, ncol(correlation_matrix))], 2))

# Create heatmap
cat("\nGenerating correlation heatmap...\n")

# Reshape for ggplot
corr_melted <- as.data.frame(as.table(correlation_matrix))
names(corr_melted) <- c("Var1", "Var2", "value")

p_heatmap <- ggplot(corr_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix of Independent Variables",
       x = "", y = "", fill = "Correlation")

print(p_heatmap)

# ==============================================================================
# 6. INITIAL VIF CALCULATION
# ==============================================================================

cat("\n### Calculating Initial VIF for Independent Variables\n")

# Initial VIF calculation
vif_formula <- as.formula(paste("NABOVE ~", paste(names(independent_vars), collapse = " + ")))
vif_model <- lm(vif_formula, data = df)
initial_vif <- vif(vif_model)

vif_data <- data.frame(
  variable = names(initial_vif),
  VIF = as.numeric(initial_vif)
) %>% arrange(desc(VIF))

cat("\n### Variance Inflation Factor (VIF) for Independent Variables:\n")
print(vif_data)

cat("\nCorrelation matrix and VIF values have been calculated and displayed.\n")

# ==============================================================================
# 7. ITERATIVE VIF REMOVAL (THRESHOLD = 5)
# ==============================================================================

cat("\n### Starting VIF analysis and variable removal...\n")

independent_vars_cleaned <- independent_vars
vif_threshold <- 5
high_vif_removed_vars <- c()

while (TRUE) {
  # Fit model with current variables
  current_formula <- as.formula(paste("NABOVE ~", 
                                      paste(names(independent_vars_cleaned), collapse = " + ")))
  current_model <- lm(current_formula, data = cbind(NABOVE = df$NABOVE, independent_vars_cleaned))
  
  # Calculate VIF
  current_vif <- vif(current_model)
  
  # Check if we need to remove any variables
  max_vif <- max(current_vif)
  
  if (max_vif > vif_threshold) {
    var_to_remove <- names(which.max(current_vif))
    independent_vars_cleaned <- independent_vars_cleaned %>% select(-all_of(var_to_remove))
    high_vif_removed_vars <- c(high_vif_removed_vars, var_to_remove)
    cat(sprintf("Removing '%s' with VIF: %.2f\n", var_to_remove, max_vif))
  } else {
    break
  }
}

cat("\n--- Final VIF Values ---\n")
final_vif <- vif(current_model)
final_vif_df <- data.frame(
  variable = names(final_vif),
  VIF = as.numeric(final_vif)
) %>% arrange(desc(VIF))

print(final_vif_df)

cat(sprintf("\nVariables removed due to high VIF: %s\n", 
            paste(high_vif_removed_vars, collapse = ", ")))
cat(sprintf("Final number of independent variables: %d\n", 
            ncol(independent_vars_cleaned)))
cat("Multicollinearity addressed. Final independent variables stored in 'independent_vars_cleaned'.\n")

# ==============================================================================
# 8. OUTLIER TREATMENT - IQR CAPPING
# ==============================================================================

cat("\n### Addressing Outliers using IQR Capping\n")

# Function to cap outliers
cap_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  
  x <- ifelse(x > upper_bound, upper_bound, x)
  x <- ifelse(x < lower_bound, lower_bound, x)
  
  return(x)
}

# Apply IQR capping to all columns
independent_vars_cleaned_capped <- independent_vars_cleaned %>%
  mutate(across(everything(), cap_outliers))

cat("Outliers addressed using IQR capping for all numerical independent variables.\n")
cat("\n### Descriptive Statistics after Outlier Capping:\n")
print(summary(independent_vars_cleaned_capped))

# ==============================================================================
# 9. PREPARE DATA FOR MODELING
# ==============================================================================

y <- df$NABOVE
X <- independent_vars_cleaned_capped

# Create adjusted y for Gamma models (must be positive)
y_gamma <- ifelse(y <= 0, 0.1, y)
if (any(y <= 0)) {
  cat(sprintf("\nNote: %d non-positive values in NABOVE adjusted to 0.1 for Gamma models\n", 
              sum(y <= 0)))
}

model_data <- cbind(NABOVE = y, X)
model_data_gamma <- cbind(NABOVE = y_gamma, X)

cat("\nData preparation complete.\n")
cat(sprintf("Final modeling data dimensions: %d rows, %d columns\n", 
            nrow(model_data), ncol(model_data)))

# ==============================================================================
# 10. FIT GLM MODELS
# ==============================================================================

cat("\n### Fitting GLM Models ###\n")

# --- Helper functions for metrics ---
r_squared <- function(y, yhat) {
  if (all(is.na(yhat)) || length(y) != length(yhat)) return(NA)
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

rmse <- function(y, yhat) {
  if (all(is.na(yhat)) || length(y) != length(yhat)) return(NA)
  sqrt(mean((y - yhat)^2))
}

mae <- function(y, yhat) {
  if (all(is.na(yhat)) || length(y) != length(yhat)) return(NA)
  mean(abs(y - yhat))
}

# Safe wrapper to catch errors in metric calculations
safe_metric <- function(metric_fun, y, yhat) {
  tryCatch(metric_fun(y, yhat), error = function(e) NA)
}

# --- Create GLM formula ---
glm_formula <- as.formula(paste("NABOVE ~", paste(names(X), collapse = " + ")))
cat(sprintf("\nGLM Formula: %s\n", deparse(glm_formula)))

# --- Storage for models and results ---
glm_models <- list()
glm_results <- list()

# --- 1. GLM Gaussian (Linear) ---
cat("\n--- Fitting GLM with Gaussian family (linear) ---\n")
tryCatch({
  glm_linear <- glm(glm_formula, data = model_data, family = gaussian())
  glm_models[['linear']] <- glm_linear
  
  predictions_linear <- predict(glm_linear, type = "response")
  
  glm_results[['linear']] <- list(
    AIC = AIC(glm_linear),
    BIC = BIC(glm_linear),
    'R-squared' = safe_metric(r_squared, y, predictions_linear),
    RMSE = safe_metric(rmse, y, predictions_linear),
    MAE = safe_metric(mae, y, predictions_linear)
  )
  
  cat("Linear GLM fitted successfully.\n")
}, error = function(e) {
  cat(sprintf("Error fitting Linear GLM: %s\n", e$message))
  glm_results[['linear']] <- list(Error = e$message)
})

# --- 2. GLM Poisson (Log link) ---
cat("\n--- Fitting GLM with Poisson family (log link) ---\n")
tryCatch({
  glm_poisson <- glm(glm_formula, data = model_data, family = poisson(link = "log"))
  glm_models[['log']] <- glm_poisson
  
  predictions_poisson <- predict(glm_poisson, type = "response")
  
  glm_results[['log']] <- list(
    AIC = AIC(glm_poisson),
    BIC = BIC(glm_poisson),
    'R-squared' = safe_metric(r_squared, y, predictions_poisson),
    RMSE = safe_metric(rmse, y, predictions_poisson),
    MAE = safe_metric(mae, y, predictions_poisson)
  )
  
  cat("Poisson GLM fitted successfully.\n")
}, error = function(e) {
  cat(sprintf("Error fitting Poisson GLM: %s\n", e$message))
  glm_results[['log']] <- list(Error = e$message)
})

# --- 3. GLM Gamma (log link) ---
cat("\n--- Fitting GLM with Gamma family (log link) ---\n")
tryCatch({
  glm_gamma <- glm(glm_formula, data = model_data_gamma, family = Gamma(link = "log"))
  glm_models[['gamma']] <- glm_gamma
  
  predictions_gamma <- predict(glm_gamma, type = "response")
  
  glm_results[['gamma']] <- list(
    AIC = AIC(glm_gamma),
    BIC = BIC(glm_gamma),
    'R-squared' = safe_metric(r_squared, y_gamma, predictions_gamma),
    RMSE = safe_metric(rmse, y_gamma, predictions_gamma),
    MAE = safe_metric(mae, y_gamma, predictions_gamma)
  )
  
  cat("Gamma GLM fitted successfully.\n")
  cat(sprintf("  AIC: %.4f, BIC: %.4f\n", 
              glm_results[['gamma']]$AIC,
              glm_results[['gamma']]$BIC))
}, error = function(e) {
  cat(sprintf("ERROR fitting Gamma GLM: %s\n", e$message))
  cat("This is likely due to convergence issues or data incompatibility with Gamma distribution.\n")
  glm_results[['gamma']] <- list(Error = e$message)
})

# --- Unified Summary of All GLMs ---
cat("\n--- Summary of All GLM Model Performance ---\n")
for (model_name in names(glm_results)) {
  cat(sprintf("\n%s GLM:\n", toupper(model_name)))
  for (metric_name in names(glm_results[[model_name]])) {
    value <- glm_results[[model_name]][[metric_name]]
    cat(sprintf("  %s: %s\n", metric_name, ifelse(is.na(value), "NA", sprintf("%.4f", value))))
  }
}

cat("\nAll GLM models fitted and performance metrics calculated.\n")

# ==============================================================================
# 11. FIT GAM MODELS
# ==============================================================================

cat("\n### Fitting GAM Models\n")

# Build GAM formula with smoothing splines
gam_terms <- c()
for (var_name in names(X)) {
  n_unique <- length(unique(X[[var_name]]))
  if (n_unique > 1) {
    k_val <- min(5, n_unique - 1)
    if (k_val >= 3) {
      gam_terms <- c(gam_terms, sprintf("s(%s, k=%d)", var_name, k_val))
    } else {
      gam_terms <- c(gam_terms, var_name)
    }
  } else {
    cat(sprintf("Skipping %s - only 1 unique value\n", var_name))
  }
}

gam_formula <- as.formula(paste("NABOVE ~", paste(gam_terms, collapse = " + ")))
cat("\nGAM Formula:\n")
print(gam_formula)

# Storage for GAM models and results
gam_models <- list()
gam_results <- list()

n_samples <- nrow(model_data)

# --- GAM Gaussian (Linear) ---
cat("\n--- Fitting GAM with Gaussian family (linear) ---\n")
tryCatch({
  gam_linear <- gam(gam_formula, data = model_data, family = gaussian())
  gam_models[['linear']] <- gam_linear
  
  predictions_linear <- predict(gam_linear, type = "response")
  
  # Calculate BIC manually
  edof_linear <- sum(gam_linear$edf)
  bic_linear <- -2 * logLik(gam_linear)[1] + log(n_samples) * edof_linear
  
  gam_results[['linear']] <- list(
    AIC = AIC(gam_linear),
    BIC = bic_linear,
    'R-squared' = r_squared(y, predictions_linear),
    RMSE = rmse(y, predictions_linear),
    MAE = mae(y, predictions_linear)
  )
  
  cat("Linear GAM fitted successfully.\n")
  cat(sprintf("  AIC: %.4f, BIC: %.4f, R²: %.4f\n", 
              gam_results[['linear']]$AIC,
              gam_results[['linear']]$BIC,
              gam_results[['linear']]$'R-squared'))
}, error = function(e) {
  cat(sprintf("Error fitting Linear GAM: %s\n", e$message))
  gam_results[['linear']] <- list(Error = e$message)
})

# --- GAM Poisson (Log) ---
cat("\n--- Fitting GAM with Poisson family (log link) ---\n")
tryCatch({
  if (any(y < 0)) {
    stop("Dependent variable 'y' contains negative values, not suitable for Poisson GAM.")
  }
  
  gam_poisson <- gam(gam_formula, data = model_data, family = poisson(link = "log"))
  gam_models[['log']] <- gam_poisson
  
  predictions_poisson <- predict(gam_poisson, type = "response")
  
  # Calculate BIC manually
  edof_poisson <- sum(gam_poisson$edf)
  bic_poisson <- -2 * logLik(gam_poisson)[1] + log(n_samples) * edof_poisson
  
  gam_results[['log']] <- list(
    AIC = AIC(gam_poisson),
    BIC = bic_poisson,
    'R-squared' = r_squared(y, predictions_poisson),
    RMSE = rmse(y, predictions_poisson),
    MAE = mae(y, predictions_poisson)
  )
  
  cat("Poisson GAM fitted successfully.\n")
  cat(sprintf("  AIC: %.4f, BIC: %.4f, R²: %.4f\n", 
              gam_results[['log']]$AIC,
              gam_results[['log']]$BIC,
              gam_results[['log']]$'R-squared'))
}, error = function(e) {
  cat(sprintf("Error fitting Poisson GAM: %s\n", e$message))
  gam_results[['log']] <- list(Error = e$message)
})

# --- GAM Gamma ---
cat("\n--- Fitting GAM with Gamma family (inverse link) ---\n")
tryCatch({
  gam_gamma <- gam(gam_formula, data = model_data_gamma, family = Gamma(link = "inverse"))
  gam_models[['gamma']] <- gam_gamma
  
  predictions_gamma <- predict(gam_gamma, type = "response")
  
  # Calculate BIC manually
  edof_gamma <- sum(gam_gamma$edf)
  bic_gamma <- -2 * logLik(gam_gamma)[1] + log(n_samples) * edof_gamma
  
  gam_results[['gamma']] <- list(
    AIC = AIC(gam_gamma),
    BIC = bic_gamma,
    'R-squared' = r_squared(y_gamma, predictions_gamma),
    RMSE = rmse(y_gamma, predictions_gamma),
    MAE = mae(y_gamma, predictions_gamma)
  )
  
  cat("Gamma GAM fitted successfully.\n")
  cat(sprintf("  AIC: %.4f, BIC: %.4f, R²: %.4f\n", 
              gam_results[['gamma']]$AIC,
              gam_results[['gamma']]$BIC,
              gam_results[['gamma']]$'R-squared'))
}, error = function(e) {
  cat(sprintf("Error fitting Gamma GAM: %s\n", e$message))
  gam_results[['gamma']] <- list(Error = e$message)
})

# Print GAM summary
cat("\n--- Summary of GAM Model Performance ---\n")
for (model_name in names(gam_results)) {
  cat(sprintf("\n%s GAM:\n", toupper(model_name)))
  for (metric_name in names(gam_results[[model_name]])) {
    value <- gam_results[[model_name]][[metric_name]]
    if (is.numeric(value)) {
      cat(sprintf("  %s: %.4f\n", metric_name, value))
    } else {
      cat(sprintf("  %s: %s\n", metric_name, value))
    }
  }
}

cat("\nGAM models fitted and performance metrics calculated.\n")

# ==============================================================================
# 12. GLM LINEAR MODEL INTERPRETATION
# ==============================================================================

cat("\n### GLM - Linear Model Interpretation\n")

if ("linear" %in% names(glm_models)) {
  glm_linear_model <- glm_models[['linear']]
  
  # Extract coefficients and p-values
  coef_summary <- summary(glm_linear_model)$coefficients
  coef_df <- data.frame(
    Predictor = rownames(coef_summary),
    Coefficient = coef_summary[, "Estimate"],
    P_value = coef_summary[, "Pr(>|t|)"]
  )
  
  # Sort by absolute coefficient value (excluding Intercept)
  coef_df_sorted <- coef_df %>%
    filter(Predictor != "(Intercept)") %>%
    arrange(desc(abs(Coefficient)))
  
  cat("\nTop Predictors by Absolute Coefficient Value (GLM - Linear):\n")
  print(coef_df_sorted)
  
  cat("\nStatistically Significant Predictors (p < 0.05) (GLM - Linear):\n")
  sig_predictors <- coef_df_sorted %>% filter(P_value < 0.05)
  if (nrow(sig_predictors) > 0) {
    print(sig_predictors)
  } else {
    cat("No statistically significant predictors at p < 0.05 level.\n")
  }
  
  cat("\n--- Notes for GLM - Linear ---\n")
  cat("* Coefficients: Positive coefficients indicate increase in predictor leads to increase in NABOVE,\n")
  cat("  while negative coefficients suggest the opposite, assuming other variables are constant.\n")
  cat("* P-values: Lower p-values (typically < 0.05) indicate statistical significance.\n")
}

# ==============================================================================
# 13. MODEL PERFORMANCE COMPARISON
# ==============================================================================

cat("\n### Model Performance Comparison (GLM vs. GAM)\n")

# Prepare comparison data frame
all_model_results <- data.frame()

for (model_type in names(glm_results)) {
  if (!("Error" %in% names(glm_results[[model_type]]))) {
    row_data <- data.frame(
      Model = sprintf("GLM - %s", toupper(model_type)),
      AIC = glm_results[[model_type]]$AIC,
      BIC = glm_results[[model_type]]$BIC,
      R_squared = glm_results[[model_type]]$`R-squared`,
      RMSE = glm_results[[model_type]]$RMSE,
      MAE = glm_results[[model_type]]$MAE
    )
    all_model_results <- rbind(all_model_results, row_data)
  }
}

for (model_type in names(gam_results)) {
  if (!("Error" %in% names(gam_results[[model_type]]))) {
    row_data <- data.frame(
      Model = sprintf("GAM - %s", toupper(model_type)),
      AIC = gam_results[[model_type]]$AIC,
      BIC = gam_results[[model_type]]$BIC,
      R_squared = gam_results[[model_type]]$`R-squared`,
      RMSE = gam_results[[model_type]]$RMSE,
      MAE = gam_results[[model_type]]$MAE
    )
    all_model_results <- rbind(all_model_results, row_data)
  }
}

# Format and display
comparison_df <- all_model_results
comparison_df[, 2:6] <- round(comparison_df[, 2:6], 4)

cat("\n╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║              Model Performance Comparison (GLM vs. GAM)                        ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n\n")

# Create a nicely formatted table
cat(sprintf("%-15s %12s %12s %12s %12s %12s\n", 
            "Model", "AIC", "BIC", "R-squared", "RMSE", "MAE"))
cat(strrep("─", 87), "\n")

for (i in 1:nrow(comparison_df)) {
  cat(sprintf("%-15s %12.4f %12.4f %12.4f %12.4f %12.4f\n",
              comparison_df$Model[i],
              comparison_df$AIC[i],
              comparison_df$BIC[i],
              comparison_df$R_squared[i],
              comparison_df$RMSE[i],
              comparison_df$MAE[i]))
}
cat(strrep("─", 87), "\n")

# ==============================================================================
# 14. IDENTIFY BEST PERFORMING MODELS
# ==============================================================================

cat("\n--- Best Performing Models ---\n")

best_aic <- comparison_df[which.min(comparison_df$AIC), ]
best_bic <- comparison_df[which.min(comparison_df$BIC), ]
best_r2 <- comparison_df[which.max(comparison_df$R_squared), ]
best_rmse <- comparison_df[which.min(comparison_df$RMSE), ]
best_mae <- comparison_df[which.min(comparison_df$MAE), ]

cat(sprintf("\nBased on AIC (lower is better): %s with AIC = %.4f\n", 
            best_aic$Model, best_aic$AIC))
cat(sprintf("Based on BIC (lower is better): %s with BIC = %.4f\n", 
            best_bic$Model, best_bic$BIC))
cat(sprintf("Based on R-squared (higher is better): %s with R-squared = %.4f\n", 
            best_r2$Model, best_r2$R_squared))
cat(sprintf("Based on RMSE (lower is better): %s with RMSE = %.4f\n", 
            best_rmse$Model, best_rmse$RMSE))
cat(sprintf("Based on MAE (lower is better): %s with MAE = %.4f\n", 
            best_mae$Model, best_mae$MAE))

cat("\n=== Analysis Complete ===\n")
