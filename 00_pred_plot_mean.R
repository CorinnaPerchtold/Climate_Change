library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(reshape2)
library(viridis)
library(viridisLite)
library(LaplacesDemon)
library(ggplot2)
library(tidyr)
library(dplyr)
library(sf)

load("01_prediction_stack_30y_long_lat.R")
load("01_elev_data_30y_long_lat.R")

index_early<-inla.stack.index(stack_early_pred, "pred_stack")$data
index_late<-inla.stack.index(stack_late_pred, "pred_stack")$data


summary_lp_gamma_early<-result4_early_pred_30y_wt_prior$summary.linear.predictor[index_early,"mean"]
summary_lp_gamma_late<-result4_late_pred_30y_wt_prior$summary.linear.predictor[index_late,"mean"]

#apply link function mu=exp(eta)
dp$expected_mean_early<-exp(summary_lp_gamma_early)
dp$expected_mean_late<-exp(summary_lp_gamma_late)

#relative difference and exp(0)=1 would mean no difference, so scale with -1 accordingly
dp$difference_mean_precip<-exp(summary_lp_gamma_late)/exp(summary_lp_gamma_early)-1

dpm_gamma_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("expected_mean_early"))
dpm_gamma_late <- melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("expected_mean_late"))
dpm_gamma_diff<-melt(dp,  id.vars = c("Longitude", "Latitude", "Month"), measure.vars = c("difference_mean_precip"))


dpm_gamma_diff<-dpm_gamma_diff %>% unite("Months", variable, Month, sep="_", remove=F) %>% arrange(Month)


months_level<-c("January", "February", "March",
                "April", "May", "June", "July",
                "August", "September", "October", "November", "December")

dpm_gamma_diff$Month<-factor(month.name[dpm_gamma_diff$Month], levels=months_level)

#in %
dpm_gamma_diff$value<-dpm_gamma_diff$value*100


library(scales)
setwd("/home/guests/corinnap/Documents/Modelling_climate_change/Data_Austria")

districts<-read_sf("gadm41_AUT_2.shp",)

setwd("~/Documents/Prediction_precip/30_year_model")

png(file="diff_mean.png",width=900)

#test_early<-melt(test_early23,  id.vars = c("Longitude.cent", "Latitude.cent", "Month"), measure.vars = c("expected_mean_early"))
#diff_mean<-
  ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_gamma_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~~Month) +
    scale_fill_gradientn(colours = c("#78410A", "#D8B365", "#F6E8C3", "white","#C7EAE5", "#5AB4AC", "#01665E"),
                         values = rescale(c(-40, -20, -10, 0, 10, 20, 40)),
                         limits = c(-50, 50),
                         oob = scales::squish, name="Relative change in\nposterior mean\nprecipitation (%)")+
  #scale_fill_viridis("Difference\nin posterior of monthly mean\nprecipitation in %",option = "turbo",limits=c(-60,50), oob=scales::squish) +
  geom_sf(data=districts, fill=NA, linewidth=0.1, color="black") +geom_sf(data=Austria, fill=NA, linewidth=0.2, color="black")+ #add borders of states
    theme(axis.text.x = element_text(size=15, face="bold"), axis.text.y = element_text(size=15, face="bold"), 
          text = element_text(size=15, face="bold") ,axis.title.y = element_text(size=15,face="bold") ,
          axis.title.x.bottom = element_text(size=15,face="bold"), 
          axis.title.x = element_text(size=15,face="bold"), strip.text.x = element_text(size=15),   #facet size
          strip.text.y = element_text(size=15))+
    scale_x_continuous(labels = function(x) round(x,3))+
    scale_y_continuous(labels=c("46.5","47.0","47.5", "48.0", "48.5","49.0"))
  
  
  dev.off()
  
  png(file="random_effect_mean_early.png", width=600)
  autoplot(result4_early_pred_30y_wt_prior)[3]
  dev.off()
  
  png(file="random_effect_mean_latey.png", width=600)
  autoplot(result4_late_pred_30y_wt_prior)[3]
  dev.off()
  
  