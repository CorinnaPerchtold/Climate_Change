library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(reshape2)
library(viridis)
library(viridisLite)
library(LaplacesDemon)
library(ggplot2)

load("01_result_pred_mean.R")
load("01_prediction_stack.R")
load("01_elev_data.R")

index_early<-inla.stack.index(stack_early_pred, "pred_stack")$data
index_late<-inla.stack.index(stack_late_pred, "pred_stack")$data

summary_lp_gamma_early<-result4_early_pred$summary.linear.predictor[index_early,"mean"]
summary_lp_gamma_late<- result4_late_pred$summary.linear.predictor[index_late, "mean"]



#apply link function mu=exp(eta)
dp$expected_mean_early<-exp(summary_lp_gamma_early)
dp$expected_mean_late<-exp(summary_lp_gamma_late)
dp$difference_mean_precip<-exp(summary_lp_gamma_late)-exp(summary_lp_gamma_early)

  
dpm_gamma_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("expected_mean_early"))
dpm_gamma_late <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("expected_mean_late"))
dpm_gamma_diff<-melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("difference_mean_precip"))


dpm_gamma_diff<-dpm_gamma_diff %>% unite("Months", variable, Month, sep="_", remove=F) %>% arrange(Month)

ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_gamma_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~~reorder(Months,Month)) +
  scale_fill_viridis("Inferred difference\nin monthly mean\nprecipitation in mm",option = "turbo", limits=c(-6,4), oob=scales::squish) +  #without nugget effect
  theme_bw()

#pdf(file="mean_diff.pdf", width=10, heigth=5)
#print(diff_mean)
#dev.off()


save(Austria,dpm_gamma_diff, file="01_plots_mean_stat_interaction.R")
