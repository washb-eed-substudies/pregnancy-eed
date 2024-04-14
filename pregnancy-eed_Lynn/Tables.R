rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

# load enrollment characteristics and results
data_path <- here("bangladesh-cleaned-master-data.RDS")
d <- readRDS(data_path)
H2 <- readRDS(here('H2_res.RDS'))
H2adj <- readRDS(here('H2_adj_res.RDS'))

#### Table 2 ####

exposure <- c("agp", "crp", "sumscore_t0_mom_Z")
outcome <- c("ln_mpo1","ln_mpo2", "ln_mpo3", "ln_aat1", "ln_aat2", "ln_aat3", "ln_neo1", "ln_neo2", "ln_neo3", "ln_reg2", "ln_L_conc_t1", "ln_L_conc_t2", "ln_L_conc_t3", "ln_M_conc_t1", "ln_M_conc_t2", "ln_M_conc_t3")
expo_var <- c("Maternal Alpha(1)-acid glycoprotein", "Maternal C-reactive protein", "Maternal Sumscore of 13 plasma cytokines")
out_var <- c("Child Myeloperoxidase Month 3", "Child  Myeloperoxidase Year 1", "Child Myeloperoxidase Year 2", "Child Alpha-1 Antitrypsin Month 3", "Child Alpha-1 Antitrypsin Year 1", "Child Alpha-1 Antitrypsin Year 2", "Child Neopetrin Month 3", "Child Neopetrin Year 1", "Child Neopetrin Year 2", "Child REG1B Year 1", "Child Lactulose Month 3", "Child Lactulose Year 1", "Child Lactulose Year 2", "Child Mannitol Month 3", "Child Mannitol Year 1", "Child Mannitol Year 2")

tbl2 <- growth_tbl("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl2flex <- growth_tbl_flex("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, H2, H2adj, T, 1.3, .7)
tbl1supp <- growth_tbl("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl1flexsupp <- growth_tbl_flex("Serum maternal inflammation biomarkers", expo_var, out_var, exposure, outcome, H2, H2adj, F, 1.3, .7)



# #### Supplementary Tables ####
# #### Table S1 ####
# 
# exposure <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "t2_iso_pca")
# outcome <- c("waz_t2", "whz_t2", "hcz_t2", "waz_t3", "whz_t3", "hcz_t3",
#              "wei_velocity_t2_t3", "hc_velocity_t2_t3",
#              "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")
# expo_var <- c("IPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI", "Combined urinary oxidative stress biomarkers")
# out_var <- c("WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
#              "WAZ Year 2", "WLZ Year 2", "HCZ Year 2",
#              "Weight velocity (kg/month) Year 1 to Year 2",
#              "Head circumference velocity (cm/month) Year 1 to Year 2",
#              "Change in child WAZ from Year 1 to Year 2",
#              "Change in WLZ from Year 1 to Year 2",
#              "Change in HCZ from Year 1 to Year 2")
# 
# tbls1 <- growth_tbl("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)
# tbls1flex <- growth_tbl_flex("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)
# 
# 
# #### Table S2 ####
# 
# exposure <- c("t3_cort_z01", "t3_cort_z03", "t3_cort_slope", "t3_residual_cort", 
#               "t3_saa_z01", "t3_saa_z02", "t3_saa_slope", "t3_residual_saa")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Pre-stressor cortisol", "Post-stressor cortisol", "Pre to post-stress change in slope of cortisol", "Cortisol residualized gain score", 
#               "Pre-stressor sAA", "Post-stressor sAA", "Pre to post-stress change in slope of sAA", "sAA residualized gain score")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")
# 
# tbls2 <- growth_tbl("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)
# tbls2flex <- growth_tbl_flex("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)
# 
# 
# #### Table S3 ####
# 
# exposure <- c("t3_map", "t3_hr_mean")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Mean arterial pressure", "Mean resting heart rate")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")
# 
# tbls3 <- growth_tbl("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)
# tbls3flex <- growth_tbl_flex("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)
# 
# 
# 
# #### Table S4 ####
# 
# exposure <- c("t3_gcr_mean", "t3_gcr_cpg12")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Entire promoter region (39 assayed CpG sites)", "NGFI-A transcription factor binding site (CpG site #12)")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")
# 
# tbls4 <- growth_tbl("Methylation site", expo_var, out_var, exposure, outcome, H4, H4adj)
# tbls4flex <- growth_tbl_flex("Methylation site", expo_var, out_var, exposure, outcome, H4, H4adj)


#### SAVE TABLES ####

write.csv(tbl1, file=here("tables/main/stress-growth-table1.csv"))
write.csv(tbl2, here('tables/main/stress-growth-table2.csv'))
write.csv(tbl3, here('tables/main/stress-growth-table3.csv'))
write.csv(tbl4, here('tables/main/stress-growth-table4.csv'))
write.csv(tbl5, here('tables/main/stress-growth-table5.csv'))
write.csv(tbl6, here('tables/main/stress-growth-table6.csv'))
write.csv(tbl7, here('tables/main/stress-growth-table7.csv'))
write.csv(tbl8, here('tables/main/stress-growth-table8.csv'))

save_as_docx("Table 1" = tbl1flex, 
             "Table 2" = tbl2flex,
             "Table 3" = tbl3flex,
             "Table 4" = tbl4flex, 
             "Table 5" = tbl5flex, 
             "Table 6" = tbl6flex, 
             "Table 7" = tbl7flex, 
             "Table 8" = tbl8flex, 
             path='C:/Users/Sophia/Documents/WASH/WASH Stress and Growth/stress-growth main v10.docx',
             pr_section = sect_properties)

write.csv(tbl1supp, here('tables/supplementary/stress-growth-tables1.csv'))
write.csv(tbl2supp, here('tables/supplementary/stress-growth-tables2.csv'))
write.csv(tbl3supp, here('tables/supplementary/stress-growth-tables3.csv'))
write.csv(tbl4supp, here('tables/supplementary/stress-growth-tables4.csv'))
write.csv(tbl5supp, here('tables/supplementary/stress-growth-tables5.csv'))
write.csv(tbl6supp, here('tables/supplementary/stress-growth-tables6.csv'))
write.csv(tbl7supp, here('tables/supplementary/stress-growth-tables7.csv'))


save_as_docx("Table S1" = tbl1flexsupp, "Table S2" = tbl2flexsupp, "Table S3" = tbl3flexsupp, 
             "Table S4" = tbl4flexsupp, "Table S5" = tbl5flexsupp, "Table S6" = tbl6flexsupp, "Table S7" = tbl7flexsupp, 
             path='C:/Users/Sophia/Documents/WASH/WASH Stress and Growth/stress-growth supplementary v10.docx',
             pr_section = sect_properties)
