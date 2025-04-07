rm(list=ls())
# Load libraries
library(here)
library(dplyr)
library(flextable)
library(officer)
library(tidyr)

# Load data (assuming you have these files)
H1_adj_res <- readRDS(here("results/adjusted/H1_adj_res.RDS")) %>% mutate(hypothesis = "H1")
H2_adj_res <- readRDS(here("results/adjusted/H2_adj_res.RDS")) %>% mutate(hypothesis = "H2")
H3_adj_res <- readRDS(here("results/adjusted/H3_adj_res.RDS")) %>% mutate(hypothesis = "H3")

# Add FDR correction
H1_adj_res <- H1_adj_res %>% mutate(corrected.pval = p.adjust(pval, method = "fdr")) 
H2_adj_res <- H2_adj_res %>% mutate(corrected.pval = p.adjust(pval, method = "fdr"))
H3_adj_res <- H3_adj_res %>% mutate(corrected.pval = p.adjust(pval, method = "fdr"))

# Combine data
df <- bind_rows(H1_adj_res, H2_adj_res, H3_adj_res) %>%
  mutate(hypothesis = factor(hypothesis, levels = c("H1", "H2", "H3")))
rownames(df) <- NULL


# Create variable labels based on the existing script
variable_labels <- c(
  "logCRP" = "Ln CRP",
  "logAGP" = "LN AGP",
  "sumscore_t0_mom_Z" = "Sum score (Z-scored)",
  "ln_preg_cort" = "Ln Cortisol (ug/dL)",
  "ln_preg_estri" = "Ln Estriol (ng/mL)",
  "ln_L_conc_t1" = "Ln Lactulose Age 3 months (mmol/L)",
  "ln_M_conc_t1" = "Ln Mannitol Age 3 months (mmol/L)",
  "ln_mpo1" = "Ln MPO Age 3 months (ng/mL)",
  "ln_aat1" = "Ln AAT Age 3 months(mg/g)",
  "ln_neo1" = "Ln NEO Age 3 months (nmol/L)",
  "ln_L_conc_t2" = "Ln Lactulose Age 14 months (mmol/L)",
  "ln_M_conc_t2" = "Ln Mannitol Age 14 months (mmol/L)",
  "ln_mpo2" = "Ln MPO Age 14 months (ng/mL)",
  "ln_aat2" = "Ln AAT Age 14 months(mg/g)",
  "ln_neo2" = "Ln NEO Age 14 months (nmol/L)",
  "ln_reg2" = "Ln REG1β Age 14 months (µg/mL)",
  "ln_L_conc_t3" = "Ln Lactulose 28 months (mmol/L)",
  "ln_M_conc_t3" = "Ln Mannitol Age 28 months (mmol/L)",
  "ln_mpo3" = "Ln MPO Age 28 months (ng/mL)",
  "ln_aat3" = "Ln AAT Age 28 months(mg/g)",
  "ln_neo3" = "Ln NEO Age 28 months (nmol/L)"
)

# Add significance stars
df <- df %>%
  mutate(
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      pval < 0.10 ~ "†",
      TRUE ~ ""
    ),
    # Create combined estimate and CI column
    estimate_with_ci = sprintf("%.3f (%.3f, %.3f)%s", est, ci.lb, ci.ub, sig)
  )

# Function to create a formatted table for each hypothesis
create_hypothesis_table <- function(data, hyp) {
  data_subset <- data %>% 
    filter(hypothesis == hyp) %>%
    mutate(
      X = ifelse(X %in% names(variable_labels), variable_labels[X], X),
      Y = ifelse(Y %in% names(variable_labels), variable_labels[Y], Y)
    )
  
  # For better formatting, let's add a bit of preprocessing for the tables
  # Sort by outcome groups (t1, t2, t3 which represent 3 months, 14 months, 28 months)
  data_subset <- data_subset %>%
    mutate(
      time_period = case_when(
        grepl("t1|1$", Y) ~ 1,
        grepl("t2|2$", Y) ~ 2,
        grepl("t3|3$", Y) ~ 3,
        TRUE ~ 4
      )
    ) %>%
    arrange(time_period, Y) %>%
    select(-time_period)
  
  # Create flextable
  ft <- data_subset %>%
    select(X, Y, estimate_with_ci, pval, corrected.pval) %>%
    flextable() %>%
    set_header_labels(
      X = "Exposure",
      Y = "Outcome",
      estimate_with_ci = "Estimate (95% CI)",
      pval = "P-value",
      corrected.pval = "FDR-adjusted P-value"
    ) %>%
    colformat_double(j = c("pval", "corrected.pval"), digits = 3) %>%
    bold(part = "header") %>%
    bg(bg = "#f5f5f5", part = "header") %>%
    border_outer() %>%
    border_inner_h() %>%
    autofit()
  
  # Add bold formatting for significant results
  for (i in 1:nrow(data_subset)) {
    if (data_subset$pval[i] < 0.05) {
      ft <- ft %>% bold(i = i, j = "estimate_with_ci")
    }
  }
  
  return(ft)
}

# Create Word document
doc <- read_docx()

# Add title and subtitle
doc <- doc %>%
  body_add_par("Pregnancy Stress Markers and Environmental Enteric Dysfunction (EED) Biomarkers", style = "heading 1") %>%
  body_add_par("Analysis of Stress Hormones and EED Outcomes", style = "heading 2") %>%
  body_add_par("", style = "Normal")

# Create and add tables for each hypothesis
for (hyp in c("H1", "H2", "H3")) {
  if (sum(df$hypothesis == hyp) > 0) {
    # Add hypothesis title
    hyp_titles <- c(
      "H1" = "Table 1: Maternal Cortisol and Child EED Biomarkers",
      "H2" = "Table 2: Maternal Stress Markers and Child EED Biomarkers",
      "H3" = "Table 3: Maternal Estriol and Child EED Biomarkers"
    )
    
    doc <- doc %>%
      body_add_par(hyp_titles[hyp], style = "heading 3") %>%
      body_add_par("", style = "Normal")
    
    # Create table
    ft <- create_hypothesis_table(df, hyp)
    
    # Add table to document
    doc <- doc %>%
      body_add_flextable(ft) %>%
      body_add_par("", style = "Normal") %>%
      body_add_par("", style = "Normal")
  }
}

# Add footnote for significance levels
doc <- doc %>%
  body_add_par("Note: † p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001", style = "Normal") %>%
  body_add_par("FDR = False Discovery Rate", style = "Normal")

# Save the Word document
print(doc, target = here("Tables/pregnancy_eed_tables.docx"))

# Display one table as an example (for console viewing)
create_hypothesis_table(df, "H1")