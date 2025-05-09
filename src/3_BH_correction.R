
rm(list=ls())

# Load libraries and configuration
source(here::here("src/0-config.R"))
source(here::here("functions/cowboy_glm.R"))

H1_res <- readRDS(here("results/unadjusted/H1_unadj_res.RDS")) %>% mutate(H=1)
H2_res <- readRDS(here("results/unadjusted/H2_unadj_res.RDS")) %>% mutate(H=1)
H3_res <- readRDS(here("results/unadjusted/H3_unadj_res.RDS")) %>% mutate(H=1)
full_res <- bind_rows(H1_res, H2_res, H3_res)


H1_adj_res <- readRDS(here("results/adjusted/H1_adj_res.RDS")) %>% mutate(H=1)
H2_adj_res <- readRDS(here("results/adjusted/H2_adj_res.RDS")) %>% mutate(H=1)
H3_adj_res <- readRDS(here("results/adjusted/H3_adj_res.RDS")) %>% mutate(H=1)
full_adj_res <- bind_rows(H1_adj_res, H2_adj_res, H3_adj_res)


full_res <- full_res %>% group_by(H) %>% 
  mutate(BH.Pval=p.adjust(pval , method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_adj_res <- full_adj_res %>% group_by(H) %>% 
  mutate(BH.Pval=p.adjust(pval , method="BH")) %>%
  ungroup() %>%
  as.data.frame()


saveRDS(full_res, here("results/unadjusted/clean_res.RDS"))
saveRDS(full_adj_res, here("results/adjusted/clean_adj_res.RDS"))

full_adj_res %>% arrange(pval )

