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
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)


source("00_functions_return_values.R")
load("01_prediction_stack_30y_long_lat.R")
load("01_elev_data_30y_long_lat.R")

index_data_early_pred<-inla.stack.index(stack_max_early_pred, "pred_stack_max" )$data
index_data_late_pred<- inla.stack.index(stack_max_late_pred, "pred_stack_max")$data


############## return values with bgev ###################

#20 years return value
return_level_early_bgev<-return_level_bgev(240, result4_max_early.pred_wt_prior$summary.linear.predictor$mean,result4_max_early.pred_wt_prior$summary.hyperpar[1,1],result4_max_early.pred_wt_prior$summary.hyperpar[2,1])
return_level_late_bgev<-return_level_bgev(240, result4_max_late.pred_wt_prior$summary.linear.predictor$mean,result4_max_late.pred_wt_prior$summary.hyperpar[1,1],result4_max_late.pred_wt_prior$summary.hyperpar[2,1])

return_level_early_bgev<-return_level_early_bgev[index_data_early_pred]
return_level_late_bgev<-return_level_late_bgev[index_data_late_pred]


dp$early<-return_level_early_bgev
dp$late<-return_level_late_bgev
dp$difference_return_values<-(return_level_late_bgev/return_level_early_bgev)-1

dpm_bgev_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("early"))
dpm_bgev_late <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("late"))
dpm_bgev_diff<-melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("difference_return_values"))


months_level<-c("January", "February", "March",
                "April", "May", "June", "July",
                "August", "September", "October", "November", "December")

dpm_bgev_early$Month<-factor(month.name[dpm_bgev_early$Month], levels=months_level)
dpm_bgev_late$Month<-factor(month.name[dpm_bgev_late$Month], levels=months_level)
dpm_bgev_diff$Month<-factor(month.name[dpm_bgev_diff$Month], levels=months_level)

dpm_bgev_diff$value<-dpm_bgev_diff$value*100

setwd("/home/guests/corinnap/Documents/Modelling_climate_change/Data_Austria")

districts<-read_sf("gadm41_AUT_2.shp",)


setwd("~/Documents/Prediction_precip/30_year_model")
png(file="diff_return_values.png",width=900)

ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_bgev_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~~Month) +
  scale_fill_gradientn(colours = c("#78410A", "#D8B365", "#F6E8C3", "white","#C7EAE5", "#5AB4AC", "#01665E"),
                       values = rescale(c(-10, -5, -2, 0, 2, 5, 10)),
                limits = c(-10,10),
              oob = scales::squish, name="Relative change \nin 20-year \nreturn values (%)")+
  #scale_fill_viridis("Difference in \nreturn values\nin mm", option="turbo",limits=c(40,80), oob=scales::squish)+
  geom_sf(data=districts, fill=NA, linewidth=0.1, color="black") +geom_sf(data=Austria, fill=NA, linewidth=0.2, color="black") +
  theme(axis.text.x = element_text(size=15, face="bold"), axis.text.y = element_text(size=15, face="bold"), 
        text = element_text(size=15, face="bold") ,axis.title.y = element_text(size=15,face="bold") ,
        axis.title.x.bottom = element_text(size=15,face="bold"), 
        axis.title.x = element_text(size=15,face="bold"), strip.text.x = element_text(size=15),   #facet size
        strip.text.y = element_text(size=15))+
  scale_x_continuous(labels = function(x) round(x,3))+
  scale_y_continuous(labels=c("46.5","47.0","47.5", "48.0", "48.5","49.0"))


dev.off()


###############bgev values
dp$early_max_value<-result4_max_early.pred_wt_prior$summary.fitted.values[index_data_early_pred,"mean"]
dp$late_max_value<-result4_max_late.pred_wt_prior$summary.fitted.values[index_data_late_pred,"mean"]
dp$difference_max_precip<-(dp$late_max_value/dp$early_max_value)-1

dpm_max_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("early_max_value"))
dpm_max_late <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("late_max_value"))
dpm_max_diff<-melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("difference_max_precip"))


months_level<-c("January", "February", "March",
                "April", "May", "June", "July",
                "August", "September", "October", "November", "December")

dpm_max_diff$Month<-factor(month.name[dpm_max_diff$Month], levels=months_level)
dpm_max_diff$value<-dpm_max_diff$value*100


setwd("~/Documents/Prediction_precip/30_year_model")
png(file="diff_daily_max_values.png",width=900)

ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_max_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~~Month) +
  scale_fill_gradientn(colours = c("#78410A", "#D8B365", "#F6E8C3", "white","#C7EAE5", "#5AB4AC", "#01665E"),
                       values = rescale(c(-20,-10,-5,0,5,10,20)),#rescale(c(-15, -7, -3, 0, 3, 7, 10)),
                       limits = c(-30, 30),#c(-20,20)
                       oob = scales::squish, name="Relative change in\nposterior daily \nmax. precipitation (%)")+
  #scale_fill_viridis("Difference in \nreturn values\nin mm", option="turbo",limits=c(40,80), oob=scales::squish)+
  geom_sf(data=districts, fill=NA, linewidth=0.1, color="black") +geom_sf(data=Austria, fill=NA, linewidth=0.2, color="black") +
  theme(axis.text.x = element_text(size=15, face="bold"), axis.text.y = element_text(size=15, face="bold"), 
        text = element_text(size=15, face="bold") ,axis.title.y = element_text(size=15,face="bold") ,
        axis.title.x.bottom = element_text(size=15,face="bold"), 
        axis.title.x = element_text(size=15,face="bold"), strip.text.x = element_text(size=15),   #facet size
        strip.text.y = element_text(size=15))+
  scale_x_continuous(labels = function(x) round(x,3))+
  scale_y_continuous(labels=c("46.5","47.0","47.5", "48.0", "48.5","49.0"))


dev.off()
