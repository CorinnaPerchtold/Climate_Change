library(dplyr)
library(tidyverse)
library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(raster)

load("01_rain_data_30y.R")
load("01_elev_data_30y.R")
load("01_stack_30y.R")


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
covars.data_pred[,2]<-dp$Slope
covars.data_pred[,3]<-dp$Aspect
colnames(covars.data_pred)<-c("Elevation", "Slope", "Aspect")

#stack for mean and minimum precipitation
stack.pred<-inla.stack(data=obs.data_pred,
                       A=list(A.pred,1),
                       effects=list(c(index,index_mean, list(Intercept=1)),
                                    list(covars.data_pred)),
                       tag="pred_stack")

#stack for maximum precipitation
stack.pred.max<-inla.stack(data=obs.data_pred,
                           A=list(A.pred,1),
                           effects=list(c(index,index_max, list(Intercept=1)),
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
     file="01_prediction_stack_30y_long_lat.R")
