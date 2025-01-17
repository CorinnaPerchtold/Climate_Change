library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(reshape2)
library(viridis)
library(viridisLite)
library(LaplacesDemon)
library(ggplot2)
library(tidyverse)

load("01_result_pred_n.binom_interaction.R")
load("01_prediction_stack.R")
load("01_elev_data.R")

index_early<-inla.stack.index(stack_early_pred, "pred_stack")$data
index_late<-inla.stack.index(stack_late_pred, "pred_stack")$data

########################### negative binomial
summary_lp_n.binom_early<-result4_n.binom_early_pred$summary.linear.predictor[index_early,"mean"]
summary_lp_n.binom_late<- result4_n.binom_late_pred$summary.linear.predictor[index_late, "mean"]

#apply link function, mu=exp(eta)
dp$expectd_dry_spell_early<-exp(summary_lp_n.binom_early)
dp$expectd_dry_spell_late<-exp(summary_lp_n.binom_late)
dp$difference_dry_spell<-exp(summary_lp_n.binom_late)-exp(summary_lp_n.binom_early)


dpm_n.binom_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("expectd_dry_spell_early"))
dpm_n.binom_late<- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("expectd_dry_spell_late"))
dpm_n.binom_diff<- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("difference_dry_spell"))


dpm_n.binom_diff<-dpm_n.binom_diff %>% unite("Months", variable, Month, sep="_", remove=F) %>% arrange(Month)

#plot
ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_n.binom_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~~reorder(Month,Month)) +
  scale_fill_viridis("Difference in\nlength of a dry spell\nin days", option="turbo",limits=c(-6,6), oob=scales::squish) +  #without nugget effect
  theme_bw()

