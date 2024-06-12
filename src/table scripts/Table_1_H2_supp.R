rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

data_path <- here("bangladesh-cleaned-master-data.RDS")
d <- readRDS(data_path)
H2 <- readRDS(here('H2_res.RDS'))
H2_adj <- readRDS(here('H2_adj_res.RDS'))

### Enrollment charcteristics
nperc <- function(vector, value){
  n <- sum(vector == value, na.rm = TRUE)
  total_non_missing <- sum(!is.na(vector))
  if (total_non_missing == 0) {
    return("No valid data")
  }
  perc <- round(n / total_non_missing * 100)
  paste(n, " (", perc, "%)", sep = "")
}

mediqr <- function(vector){
  quantiles <- round(quantile(vector, na.rm=T), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

n_med_col <- c(nperc(d$sex, "female"),
               mediqr(d$laz_t1), mediqr(d$waz_t1), mediqr(d$whz_t1), mediqr(d$hcz_t1),
               mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
               mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3), 
               nperc(d$diar7d_t1, 1), nperc(d$diar7d_t2, 1), nperc(d$diar7d_t3, 1), 
               nperc(d$ari7d_t2, 1), nperc(d$ari7d_t3, 1),
               mediqr(d$momage), mediqr(d$gest_age_weeks), mediqr(d$momheight), 
               mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
               nperc(d$life_viol_any_t3, 1))

tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","Mother","","","","","","",""),
                   "C2" = c("","Anthropometry (3 months)","","","",
                            "Anthropometry (14 months)","","","",
                            "Anthropometry (28 months)","","","",
                            "Diarrhea (3 months)", "Diarrhea (14 months)", "Diarrhea (28 months)",
                            "Acute respiratory illness (14 months)", "Acute respiratory illness (28 months)", "","",
                            "Anthropometry at enrollment", "Education", "Depression (14 months)", "Depression (28 months)", "Perceived stress (28 months)", 
                            "Intimate partner violence"),
                   "C3" = c("Female",
                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", 
                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", 
                            "Age (years)", "Gestational age (weeks)","Height (cm)", "Schooling completed (years)",
                            "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
                   "C4" = n_med_col)


tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
tbl1flex <- set_header_labels(tbl1flex,
                              values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
tbl1flex <- autofit(tbl1flex, part = "all")
tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
tbl1flex <- fit_to_width(tbl1flex, max_width=8)
names(tbl1)<- c("","","","n (%) or median (IQR)")


### Table H2
BH_H2 <- H2 %>% 
  mutate(BH.Pval = p.adjust(Pval, method="BH")) %>%
  as.data.frame()

BH_H2_adj <- H2_adj %>% 
  mutate(BH.Pval = p.adjust(Pval, method="BH")) %>%
  as.data.frame()

saveRDS(BH_H2, file = here("BH_H2.RDS"))
saveRDS(BH_H2_adj, file = here("BH_H2_adj.RDS"))

exposure <- c("agp", "crp", "sumscore_t0_mom_Z")
outcome <- c("ln_L_conc_t1", "ln_M_conc_t1", "ln_mpo1", "ln_aat1", "ln_neo1", "ln_L_conc_t2", "ln_M_conc_t2", "ln_mpo2", "ln_aat2", "ln_neo2", "ln_reg2", "ln_L_conc_t3", "ln_M_conc_t3","ln_mpo3", "ln_aat3", "ln_neo3")
expo_var <- c("Ln Alpha(1)-acid glycoprotein (g/L)", "Ln C-reactive protein (mg/L)", "Maternal Sumscore of 13 plasma cytokines")
out_var <- c("Ln Lactulose 3 months (mmol/L)", "Ln Mannitol 3 months (mmol/L)", "Ln MPO Age 3 months (ng/mL)", "Ln A1AT Age 3 months(mg/g)", "Ln NEO Age 3 months (nmol/L)", "Ln Lactulose 14 months (mmol/L)", "Ln Mannitol 14 months (mmol/L)", "Ln MPO Age 14 months (ng/mL)", "Ln A1AT Age 14 months(mg/g)", "Ln NEO Age 14 months (nmol/L)", "Ln REG1β Age 14 months (µg/mL)", "Ln Lactulose 28 months (mmol/L)", "Ln Mannitol 28 months (mmol/L)", "Ln MPO Age 28 months (ng/mL)", "Ln A1AT Age 28 months(mg/g)", "Ln NEO Age 28 months (nmol/L)")

tbl2 <- growth_tbl("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, BH_H2, BH_H2_adj, T)
tbl2flex <- growth_tbl_flex("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, BH_H2, BH_H2_adj, T, 1.3, .7)
tbl2supp <- growth_tbl("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, BH_H2, BH_H2_adj)
tbl2flexsupp <- growth_tbl_flex("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, BH_H2, BH_H2_adj, F, 1.3, .7)



#### SAVE TABLES ####

write.csv(tbl1, here("tbl1.csv"))
write.csv(tbl2, here("tbl2.csv"))


save_as_docx("Table 1" = tbl1flex, 
             "Table 2" = tbl2flex,
             path=here("tbl-H2.docx"),
             pr_section = sect_properties)

write.csv(tbl2supp, here("tbl2supp.csv"))


save_as_docx("Table S2" = tbl2flexsupp,
             path=here("tbl-H2-supp.docx"),
             pr_section = sect_properties)

