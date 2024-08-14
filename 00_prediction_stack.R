library(dplyr)
library(tidyverse)
library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(raster)

load("01_rain_data.R")
load("01_elev_data.R")
load("01_stack.R")


#elevation map
Aut.elev.aggregate<-rasterToPoints(Aut.elev.aggregate)
Aut.elev.aggregate<-as.data.frame(Aut.elev.aggregate)
colnames(Aut.elev.aggregate)<-c("Longitude", "Latitude", "Elevation")
Aut.elev.aggregate<- Aut.elev.aggregate %>% mutate( Latitude.cent=Latitude-mean(temp_all$Latitude),
                                                    Longitude.cent=Longitude-mean(temp_all$Longitude),
                                                    Elevation=Elevation/1000)

#replicate data frame
dp<-sapply(Aut.elev.aggregate,rep.int,times=max(temp_early$Time))
dp<-as.data.frame(dp)

#add time and month to prediction data frame
dp$Time<-rep(1:max(temp_early$Time), each=length(Aut.elev.aggregate$Longitude))
dp$Month<-rep(1:max(temp_early$Month), each=length(Aut.elev.aggregate$Longitude))

#set the amount of days per month
dp$Days_per_month[dp$Month==1| dp$Month ==3 | dp$Month==5 |dp$Month==7|dp$Month==8|dp$Month==10 |dp$Month==12]<-31
dp$Days_per_month[dp$Month==2]<-28
dp$Days_per_month[dp$Month==4|dp$Month==6 |dp$Month==9| dp$Month==11]<-30

dp$Days_per_month<-as.integer(dp$Days_per_month)

#preparation for binomial distribution
Ntrials_early<-c(temp_early$Days_per_month, dp$Days_per_month)
Ntrials_late<-c(temp_late$Days_per_month, dp$Days_per_month)

#preparation for bgev 
spread.x_early<- c(temp_early$Elevation, dp$Elevation)
spread.x_late<- c(temp_late$Elevation, dp$Elevation)

n_early_pred<- length(c(temp_early$Elevation, dp$Elevation))
n_late_pred<- length(c(temp_late$Elevation, dp$Elevation))

null.matrix_early_pred<- matrix(nrow=n_early_pred, ncol=0)
null.matrix_late_pred<- matrix(nrow=n_late_pred, ncol=0)

tail.x_early<-null.matrix_early_pred
tail.x_late<-null.matrix_late_pred

#prediction points
pred.points <- cbind(dp$Longitude,dp$Latitude)

#observation matrix
A.pred<-inla.spde.make.A(mesh,loc=pred.points, group=dp$Month)

#response is set to NA at prediction points
obs.data_pred<- data.frame(matrix(data=NA, ncol=5, nrow=length(dp$Longitude)))
colnames(obs.data_pred)<-c("Rain_no_days", "Rain_no_length", "Rain_max", "Rain_mean", "Rain_sum")

covars.data_pred<-data.frame(matrix(data=NA, ncol=3, nrow=length(dp$Longitude)))
covars.data_pred[,1]<-dp$Elevation
covars.data_pred[,2]<-dp$Latitude.cent
covars.data_pred[,3]<-dp$Longitude.cent
colnames(covars.data_pred)<-c("Elevation", "Latitude.cent", "Longitude.cent")

#stack for mean and minimum precipitation
stack.pred<-inla.stack(data=obs.data_pred,
                       A=list(A.pred,1),
                       effects=list(c(index_mean, list(Intercept=1)),
                                    list(covars.data_pred)),
                       tag="pred_stack")

#stack for maximum precipitation
stack.pred.max<-inla.stack(data=obs.data_pred,
                           A=list(A.pred,1),
                           effects=list(c(index_max, list(Intercept=1)),
                                        list(covars.data_pred)),
                           tag="pred_stack_max")


#combine stacks for mean and min
stack_early_pred<-inla.stack(stack_early, stack.pred)
stack_late_pred<-inla.stack(stack_late, stack.pred)

#combine stacks for max
stack_max_early_pred<-inla.stack(stack_max_early, stack.pred.max)
stack_max_late_pred<-inla.stack(stack_max_late, stack.pred.max)


save(dp, A.pred,
     Ntrials_early, Ntrials_late,
     spread.x_early, spread.x_late, tail.x_early,tail.x_late, 
     stack_early_pred, stack_late_pred, 
     stack_max_early_pred, stack_max_late_pred, 
     file="01_prediction_stack.R")

