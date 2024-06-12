library(dplyr)
library(tidyr)
library(flextable)
library(data.table)

mediqr1 <- function(vector){
  quantiles <- round(quantile(vector, na.rm = TRUE), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep = "")}

mediqr <- function(vector){
  # Calculate quantiles while ignoring NA values and round them
  quantiles <- round(quantile(vector, na.rm = TRUE), 2)
  
  # Convert quantiles back from log scale to normal scale using exp()
  quantiles <- exp(quantiles)
  
  # Format the median and interquartile range (Q1, Q3) as a string
  formatted_output <- paste(round(quantiles[3], 2), " (", round(quantiles[2], 2), ", ", round(quantiles[4], 2), ")", sep = "")
  
  return(formatted_output)
}

# Maternal exposures
maternal_col <- c(mediqr(d$ln_preg_estri), mediqr(d$ln_preg_cort), mediqr(d$agp), mediqr(d$crp), mediqr1(d$sumscore_t0_mom_Z))

# Create the data frame with check.names = FALSE to prevent automatic renaming
tblm <- data.frame("Maternal Biomarker At Enrollment" = c("Estriol (ng/mL)", "Cortisol (ug/dL)", 
                                                         "Alpha(1)-acid glycoprotein (g/L)", "C-reactive protein (mg/L)", 
                                                         "Cytokine sum score"),
                  "Median (25th, 75th percentile)" = maternal_col,
                  check.names = FALSE)

# Create flextable
tblmflex <- flextable(tblm)
tblmflex <- set_header_labels(tblmflex, values = list(C1 = "Maternal Biomarker At Enrollment", C2 = "Median (25th, 75th percentile)"))
tblmflex <- hline_top(tblmflex, part = "header", border = fp_border(color = "black", width = 1))
tblmflex <- hline_bottom(tblmflex, part = "all", border = fp_border(color = "black", width = 1))
tblmflex <- autofit(tblmflex, part = "all")
tblmflex <- align(tblmflex, j = c(1, 2), align = "left", part = "all")
tblmflex <- align(tblmflex, j = 2, align = "center", part = "all")
tblmflex <- fit_to_width(tblmflex, max_width = 8)

# Child outcomes 

child_col_t1 <- c("Median (25th, 75th percentile)", mediqr(d$ln_L_conc_t1), mediqr(d$ln_M_conc_t1), mediqr(d$ln_mpo1), mediqr(d$ln_aat1), mediqr(d$ln_neo1), "N/A")
child_col_t2 <- c("Median (25th, 75th percentile)",mediqr(d$ln_L_conc_t2), mediqr(d$ln_M_conc_t2), mediqr(d$ln_mpo2), mediqr(d$ln_aat2), mediqr(d$ln_neo2), mediqr(d$ln_reg2))
child_col_t3 <- c("Median (25th, 75th percentile)",mediqr(d$ln_L_conc_t3), mediqr(d$ln_M_conc_t3), mediqr(d$ln_mpo3), mediqr(d$ln_aat3), mediqr(d$ln_neo3), "N/A")

# Create the data frame with check.names = FALSE to prevent automatic renaming
tblc <- data.frame("Child Biomarker" = c("Overview","Lactulose (mmol/L)", "Mannitol (mmol/L)", "MPO (ng/mL)",
                                         "A1AT (mg/g)", "NEO (nmol/L)", "REG1β (µg/mL)"),
                   "Age 3 months" = child_col_t1,
                   "Age 14 months" = child_col_t2,
                   "Age 28 months" = child_col_t3,
                   check.names = FALSE)

# Create flextable
tblcflex <- flextable(tblc)
tblcflex <- hline_top(tblcflex, part = "header", border = fp_border(color = "black", width = 1))
tblcflex <- hline_bottom(tblcflex, part = "all", border = fp_border(color = "black", width = 1))
tblcflex <- autofit(tblcflex, part = "all")
tblcflex <- align(tblcflex, j = c(1, 2, 3), align = "left", part = "all")
tblcflex <- align(tblcflex, j = c(2, 3, 4), align = "center", part = "all")
tblcflex <- fit_to_width(tblcflex, max_width = 8)

# Save


write.csv(tblm, here("tblm.csv"))
write.csv(tblc, here("tblc.csv"))

save_as_docx("Maternal Biomarkers" = tblmflex,
             "Child Biomarkers" = tblcflex,
             path=here("Maternal+Child-Biomarkers.docx"),
             pr_section = sect_properties)

