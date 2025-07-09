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
library(dplyr)
library(sf)
library(scales)

load("01_prediction_stack_30y_long_lat.R")
load("01_elev_data_30y_long_lat.R")

index_early<-inla.stack.index(stack_early_pred, "pred_stack")$data
index_late<-inla.stack.index(stack_late_pred, "pred_stack")$data

########################### negative binomial
summary_lp_n.binom_early<-result4_n.binom_early_pred_30y_wt_prior$summary.linear.predictor[index_early,"mean"]
summary_lp_n.binom_late<- result4_n.binom_late_pred_30y_wt_prior$summary.linear.predictor[index_late, "mean"]

#apply link function, mu=exp(eta)
dp$expectd_dry_spell_early<-exp(summary_lp_n.binom_early)
dp$expectd_dry_spell_late<-exp(summary_lp_n.binom_late)
dp$difference_dry_spell<-exp(summary_lp_n.binom_late)/exp(summary_lp_n.binom_early)-1 #to subtract 100 % so to say


dpm_n.binom_early <- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("expectd_dry_spell_early"))
dpm_n.binom_late<- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("expectd_dry_spell_late"))
dpm_n.binom_diff<- melt(dp,  id.vars = c("Longitude", "Latitude", "Time","Month"), measure.vars = c("difference_dry_spell"))


dpm_n.binom_diff<-dpm_n.binom_diff %>% unite("Months", variable, Month, sep="_", remove=F) %>% arrange(Month)


months_level<-c("January", "February", "March",
           "April", "May", "June", "July",
              "August", "September", "October", "November", "December")

dpm_n.binom_diff$Month<-factor(month.name[dpm_n.binom_diff$Month], levels=months_level)


#in %
dpm_n.binom_diff$value<-dpm_n.binom_diff$value*100
setwd("/home/guests/corinnap/Documents/Modelling_climate_change/Data_Austria")

districts<-read_sf("gadm41_AUT_2.shp",)


setwd("~/Documents/Prediction_precip/30_year_model")

png(file="diff_dry_spell.png",width=900)


#plot
ggplot(Austria) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpm_n.binom_diff, aes(x = Longitude, y = Latitude, fill = value)) +
  labs(x = "", y = "") +
  facet_wrap(~Month)+
  #scale_fill_viridis("Relative difference in\nposterior maximum \nlength of a dry \nspell in %", option="turbo",limits=c(0.7,1.2), oob=scales::squish) +  #without nugget effect
  #scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0)+
  scale_fill_gradientn(colours = c("#2166AC", "#67A9CF", "#D1E5F0", "white","#FDDBC7", "#EF8A62", "#B2182B"),
                       values = rescale(c(-40, -20, -10, 0, 10, 20, 40)),
                       limits = c(-30, 30),
                       oob = scales::squish, name="Relative change\nin posterior max.\nlength of a\ndry spell (%)")+
  geom_sf(data=districts, fill=NA, linewidth=0.1, color="black") +geom_sf(data=Austria, fill=NA, linewidth=0.2, color="black") +
  theme(axis.text.x = element_text(size=15, face="bold"), axis.text.y = element_text(size=15, face="bold"), 
        text = element_text(size=15, face="bold") ,axis.title.y = element_text(size=15,face="bold") ,
        axis.title.x.bottom = element_text(size=15,face="bold"), 
        axis.title.x = element_text(size=15,face="bold"), strip.text.x = element_text(size=15),   #facet size
        strip.text.y = element_text(size=15))+
  scale_x_continuous(labels = function(x) round(x,3))+
  scale_y_continuous(labels=c("46.5","47.0","47.5", "48.0", "48.5","49.0"))

dev.off()

setwd("~/Documents/Prediction_precip/final pics")
#ggsave(file="dry_spell_diff.png", height = 5, width=10)

dev.off()

png(file="random_effect_n.binom_early.png", width=600)
autoplot(result4_n.binom_early_pred_30y_wt_prior)[3]
dev.off()

png(file="random_effect_n.binom_later.png", width=600)
autoplot(result4_n.binom_late_pred_30y_wt_prior)[3]
dev.off()

