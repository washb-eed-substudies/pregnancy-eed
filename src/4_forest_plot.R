
rm(list=ls())

# Load libraries and configuration
source(here::here("src/0-config.R"))
source(here::here("functions/cowboy_glm.R"))

res <- readRDS(here("results/adjusted/clean_adj_res.RDS"))

ggplot(res, aes(x=Y,y=est)) + 
  geom_point() +
  geom_linerange(aes(ymin=ci.lb, ymax=ci.ub)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_flip() +
  facet_wrap(~X, scales="free")
