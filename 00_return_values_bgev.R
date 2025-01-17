library(ggplot2)
library(inlabru)
library(sf)
library(RColorBrewer)
library(INLA)
library(foreach)
library(Matrix)
library(foreach)
library(viridis)
library(viridisLite)

source("00_functions_return_values.R")
load("01_elev_data.R")
load("01_prediction_stack.R")

index_data_early_pred<-inla.stack.index(stack_max_early_pred, "pred_stack_max" )$data
index_data_late_pred<- inla.stack.index(stack_max_late_pred, "pred_stack_max")$data


############## return values with bgev ###################

#20 years return value
return_level_early_bgev<-return_level_bgev(240, result4_max_early.pred$summary.linear.predictor$mean,result4_max_early.pred$summary.hyperpar[1,1],result4_max_early.pred$summary.hyperpar[2,1])
return_level_late_bgev<-return_level_bgev(240, result4_max_late.pred$summary.linear.predictor$mean,result4_max_late.pred$summary.hyperpar[1,1],result4_max_late.pred$summary.hyperpar[2,1])

return_level_early_bgev<-return_level_early_bgev[index_data_early_pred]
return_level_late_bgev<-return_level_late_bgev[index_data_late_pred]


dp$early<-return_level_early_bgev
dp$late<-return_level_late_bgev
dp$difference_return_values<-return_level_late_bgev-return_level_early_bgev

dp$var_early<-rep(factor("return_level_early"), length(return_level_early_bgev))
dp$var_late<-rep(factor("return_level_late"), length(return_level_early_bgev))
dp$variable<-rep(factor("difference_return_values"), length(return_level_early_bgev))

dp<-dp %>% unite("Months", variable, Month, sep="_", remove=F) %>% arrange(Month)


ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dp, aes(x = Longitude, y = Latitude, fill = difference_return_values)) +
  labs(x = "", y = "") +
  facet_wrap(~~reorder(Months,Month)) +
  scale_fill_viridis("Difference in \nreturn values\nin mm", option="turbo",limits=c(-20,5), oob=scales::squish) 
