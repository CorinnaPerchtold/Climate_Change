library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(raster)

load("01_elev_data_30y_long_lat.R")
load("01_stack_30y_long_lat.R")
load("01_prediction_stack_30y_long_lat.R")

CV<-function(data_raw_early,data_raw_late,   data_inferred_early,data_inferred_late){
  early<-data.frame(Log_score=-round(mean(log(data_inferred_early$cv), na.rm=T),3),
                    Mean_square_error=round(mean((data_raw_early-data_inferred_early$mean)^2),4),
                    RMSE=round(sqrt(mean((data_raw_early-data_inferred_early$mean)^2)),4),
                    Mean_absolute_error=round(mean(abs(data_raw_early-data_inferred_early$mean)),4),
                    Period="Early")
  
  late<-data.frame(Log_score=-round(mean(log(data_inferred_late$cv), na.rm=T),3),
                   Mean_square_error=round(mean((data_raw_late-data_inferred_late$mean)^2),4),
                   RMSE=round(sqrt(mean((data_raw_late-data_inferred_late$mean)^2)),4),
                   Mean_absolute_error=round(mean(abs(data_raw_late-data_inferred_late$mean)),4),
                   Period="Late")
  
  rbind(early,late)
}



####################### mean CV #####################
formula4_mean<-Rain_mean~ -1+Intercept+scale(Elevation):Month+Slope+scale(Aspect)+
  f(spatial_field, model=spde.nonstat)+
  f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group,
    control.group=list(model="ar1", hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))


result4_early_30y<-inla(formula4_mean, family="gamma",
                        data=inla.stack.data(stack_early, spde=spde.nonstat),
                        control.family = list(link='log'), 
                        control.predictor = list(A=inla.stack.A(stack_early),compute=TRUE,link=1),
                        control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                               openmp.strategy="pardiso.parallel"),
                        verbose=TRUE)

result4_late_30y<-inla(formula4_mean, family="gamma",
                       data=inla.stack.data(stack_late, spde=spde.nonstat),
                       control.family = list(link='log'),
                       control.predictor = list(A=inla.stack.A(stack_late), compute=TRUE,link=1),
                       control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                              openmp.strategy="pardiso.parallel"),
                       verbose=TRUE)#, inla.mode = "experimental")


inla.stack.index(stack_early,"mean_early")$data->index_early
result4_early_30y$summary.fitted$mean[index_early]->result4_early_30y_data

inla.stack.index(stack_late,"mean_late")$data->index_late
result4_late_30y$summary.fitted$mean[index_late]->result4_late_30y_data

#leave one out cross validation
loocv_early<-inla.group.cv(result=result4_early_30y, num.level.sets=-1)
loocv_late<-inla.group.cv(result=result4_late_30y, num.level.sets=-1)

#leave group out cross validation
lgocv_early_3<-inla.group.cv(result=result4_early_30y, num.level.sets=3, strategy="posterior")
lgocv_late_3<-inla.group.cv(result=result4_late_30y, num.level.sets=3, strategy="posterior")

LOO_mean<- CV(result4_early_30y_data, result4_late_30y_data,loocv_early, loocv_late )
LGO_mean<-CV(result4_early_30y_data, result4_late_30y_data,lgocv_early_3, lgocv_late_3 )

LOO_mean$CV<-"Single-leave-out"
LGO_mean$CV<-"Leave-3 groups-out"

final_CV_mean<-rbind(LOO_mean, LGO_mean)


####################### n.binom CV #####################
E_early<-rep(1, length(temp_early$Days_per_month))
E_late<-rep(1, length(temp_late$Days_per_month))

formula4_dry_spell<-Rain_no_length ~ -1 + Intercept + scale(Elevation):Month + Slope + 
  scale(Aspect) + f(spatial_field, model = spde.nonstat) + 
  f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                                                                                                     hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))


result4_n.binom_early_30y<-inla(formula4_dry_spell,  family="nbinomial",E=E_early,
                                data=inla.stack.data(stack_early, spde=spde.nonstat),
                                control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                control.predictor = list(A=inla.stack.A(stack_early), compute=TRUE,link=1),
                                control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                       openmp.strategy="pardiso"),
                                verbose=TRUE)

result4_n.binom_late_30y<-inla(formula4_dry_spell,  family="nbinomial", E=E_late,
                               data=inla.stack.data(stack_late, spde=spde.nonstat),
                               control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                               control.predictor = list(A=inla.stack.A(stack_late), compute=TRUE,link=1),
                               control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                      openmp.strategy="pardiso"),
                               verbose=TRUE)

result4_n.binom_early_30y$summary.fitted$mean[index_early]->result4_n.binom_early_30y_data
result4_n.binom_late_30y$summary.fitted$mean[index_late]->result4_n.binom_late_30y_data


#leave one out cross validation
loocv_early_n.binom<-inla.group.cv(result=result4_n.binom_early_30y, num.level.sets=-1)
loocv_late_n.binom<-inla.group.cv(result=result4_n.binom_late_30y, num.level.sets=-1)

lgocv_early_3_n.binom<-inla.group.cv(result=result4_n.binom_early_30y, num.level.sets=3, strategy="posterior")
lgocv_late_3_n.binom<-inla.group.cv(result=result4_n.binom_late_30y, num.level.sets=3, strategy="posterior")


LOO_n.binom<-CV(result4_early_30y_data, result4_late_30y_data,loocv_early_n.binom, loocv_late_n.binom)
LGO_n.binom<- CV(result4_early_30y_data, result4_late_30y_data,lgocv_early_3_n.binom, lgocv_late_3_n.binom )
LOO_n.binom$CV<-"Single-leave-out"
LGO_n.binom$CV<-"Leave-3 groups-out"

final_CV_n.binom<-rbind(LOO_n.binom, LGO_n.binom)

####################### bgev CV #####################
formula4_max_early<-inla.mdata(Rain_max, spread.x_early, tail.x_early) ~ -1+Intercept +scale(Elevation):Month+Slope+scale(Aspect)+
  f(spatial_field, model=spde.nonstat)+ f(spatio_temporal_field, model=spde.nonstat, group = spatio_temporal_field.group, 
                                          control.group = list(model="ar1",hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))

result4_max_early_30y<-inla(formula4_max_early, family="bgev",
                            data=inla.stack.data(stack_max_early,spde=spde.nonstat),
                            control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                            control.predictor = list(A=inla.stack.A(stack_max_early),compute=TRUE,link=1),
                            control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                   openmp.strategy="pardiso.parallel"),
                            verbose=TRUE, 
                            control.inla = list( cmin=1e-6))


formula4_max_late<- inla.mdata(Rain_max, spread.x_late, tail.x_late) ~ -1 + Intercept + 
  scale(Elevation):Month + Slope + scale(Aspect) + f(spatial_field, 
                                                     model = spde.nonstat) + f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                                                                                                                                                                                hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))

result4_max_late_30y<-inla(formula4_max_late, family="bgev", 
                           data=inla.stack.data(stack_max_late,spde=spde.nonstat),
                           control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                           control.predictor = list(A=inla.stack.A(stack_max_late),compute=TRUE,link=1),
                           control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                  openmp.strategy="pardiso.parallel"),
                           verbose=TRUE,
                           control.inla =list(cmin=1e-6))

result4_max_early_30y$summary.fitted$mean[index_early]->result4_max_early_30y_data
result4_max_late_30y$summary.fitted$mean[index_late]->result4_max_late_30y_data


#leave one out cross validation
loocv_early_max<-inla.group.cv(result=result4_max_early_30y, num.level.sets=-1)
loocv_late_max<-inla.group.cv(result=result4_max_late_30y, num.level.sets=-1)

lgocv_early_3_max<-inla.group.cv(result=result4_max_early_30y, num.level.sets=3, strategy="posterior")
lgocv_late_3_max<-inla.group.cv(result=result4_max_late_30y, num.level.sets=3, strategy="posterior")


LOO_max<- CV(result4_early_30y_data, result4_late_30y_data,loocv_early_max, loocv_late_max )
LGO_max<-CV(result4_early_30y_data, result4_late_30y_data,lgocv_early_3_max, lgocv_late_3_max )

LOO_max$CV<-"Single-leave-out"
LGO_max$CV<-"Leave-3 groups-out"

final_CV_max<-rbind(LOO_max, LGO_max)

save(final_CV_mean, final_CV_n.binom, final_CV_max, file="01_Crossvalidation_30y.R")



############## check how well inference is in high mountain areas mean ###########
A_early<-inla.spde.make.A(mesh, loc=as.matrix(temp_early[,c('Longitude','Latitude')]), 
                          group = temp_early$Month)

A_late<-inla.spde.make.A(mesh, loc=as.matrix(temp_late[,c('Longitude','Latitude')]), 
                         group = temp_late$Month)


cv_results_early<-list()
cv_results_late<-list()

temp_early->cv_data_early
temp_late->cv_data_late

#remove Rain values for elevation above 2000m
cv_data_early$Rain_mean[cv_data_early$Elevation>=2]<-NA
cv_data_late$Rain_mean[cv_data_late$Elevation>=2]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))


stack.train_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                              effects=list(c(index,index_mean, list(Intercept=1)),
                                           list(covars.data_early)),
                              tag="train_early")

stack.train_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                             effects=list(c(index,index_mean, list(Intercept=1)),
                                          list(covars.data_late)),
                             tag="train_late")


cv_fit_elevation_above2_early<-inla(formula4_mean, family="gamma",
                                    data=inla.stack.data(stack.train_early, spde=spde.nonstat),
                                    control.family = list(link='log'), 
                                    control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                    control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                           openmp.strategy="pardiso.parallel"),
                                    verbose=TRUE)

cv_fit_elevation_above2_late<-inla(formula4_mean, family="gamma",
                                   data=inla.stack.data(stack.train_late, spde=spde.nonstat),
                                   control.family = list(link='log'), 
                                   control.predictor = list(A=inla.stack.A(stack.train_late),compute=TRUE,link=1),
                                   control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                          openmp.strategy="pardiso.parallel"),
                                   verbose=TRUE)

which(temp_early$Elevation>=2)->index_elev_early
which(temp_late$Elevation>=2)->index_elev_late

cv_results_early[[1]]<-list(fold=1, observed=temp_early$Rain_mean[index_elev_early], predicted_mean=cv_fit_elevation_above2_early$summary.fitted.values[index_elev_early,"mean"], predicted_sd=cv_fit_elevation_above2_early$summary.fitted.values[index_elev_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_above2_early$cpo$cpo), na.rm=T),3))
cv_results_late[[1]]<-list(fold=1, observed=temp_late$Rain_mean[index_elev_late], predicted_mean=cv_fit_elevation_above2_late$summary.fitted.values[index_elev_late,"mean"],  predicted_sd=cv_fit_elevation_above2_late$summary.fitted.values[index_elev_late,"sd"],log_score=-round(mean(log(cv_fit_elevation_above2_late$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation between 1000m and 2000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_mean[cv_data_early$Elevation<2 & cv_data_early$Elevation>=1]<-NA
cv_data_late$Rain_mean[cv_data_late$Elevation<2 & cv_data_late$Elevation>=1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))


stack.train1_2_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                 effects=list(c(index,index_mean, list(Intercept=1)),
                                              list(covars.data_early)),
                                 tag="train_early")

stack.train1_2_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_late)),
                                tag="train_late")

cv_fit_elevation_between1_2_early<-inla(formula4_mean, family="gamma",
                                        data=inla.stack.data(stack.train1_2_early, spde=spde.nonstat),
                                        control.family = list(link='log'), 
                                        control.predictor = list(A=inla.stack.A(stack.train1_2_early),compute=TRUE,link=1),
                                        control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                               openmp.strategy="pardiso.parallel"),
                                        verbose=TRUE)

cv_fit_elevation_between1_2_late<-inla(formula4_mean, family="gamma",
                                       data=inla.stack.data(stack.train1_2_late, spde=spde.nonstat),
                                       control.family = list(link='log'), 
                                       control.predictor = list(A=inla.stack.A(stack.train1_2_late),compute=TRUE,link=1),
                                       control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                              openmp.strategy="pardiso.parallel"),
                                       verbose=TRUE)

which(temp_early$Elevation<2 & temp_early$Elevation >=1)->index_elev_1_2_early
which(temp_late$Elevation<2 & temp_late$Elevation >=1)->index_elev_1_2_late

cv_results_early[[2]]<-list(fold=2, observed=temp_early$Rain_mean[index_elev_1_2_early], predicted_mean=cv_fit_elevation_between1_2_early$summary.fitted.values[index_elev_1_2_early,"mean"], predicted_sd=cv_fit_elevation_between1_2_early$summary.fitted.values[index_elev_1_2_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_early$cpo$cpo), na.rm=T),3))
cv_results_late[[2]]<-list(fold=2, observed=temp_late$Rain_mean[index_elev_1_2_late], predicted_mean=cv_fit_elevation_between1_2_late$summary.fitted.values[index_elev_1_2_late,"mean"], predicted_sd=cv_fit_elevation_between1_2_late$summary.fitted.values[index_elev_1_2_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_late$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation below 1000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_mean[cv_data_early$Elevation<1]<-NA
cv_data_late$Rain_mean[cv_data_late$Elevation<1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))

stack.train_1_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_early)),
                                tag="train_early")

stack.train_1_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                               effects=list(c(index,index_mean, list(Intercept=1)),
                                            list(covars.data_late)),
                               tag="train_late")

cv_fit_elevation_below_1_early<-inla(formula4_mean, family="gamma",
                                     data=inla.stack.data(stack.train_1_early, spde=spde.nonstat),
                                     control.family = list(link='log'), 
                                     control.predictor = list(A=inla.stack.A(stack.train_1_early),compute=TRUE,link=1),
                                     control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                            openmp.strategy="pardiso.parallel"),
                                     verbose=TRUE)


cv_fit_elevation_below_1_late<-inla(formula4_mean, family="gamma",
                                    data=inla.stack.data(stack.train_1_late, spde=spde.nonstat),
                                    control.family = list(link='log'), 
                                    control.predictor = list(A=inla.stack.A(stack.train_1_late),compute=TRUE,link=1),
                                    control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                           openmp.strategy="pardiso.parallel"),
                                    verbose=TRUE)

which(temp_early$Elevation<1 )->index_elev_1_early
which(temp_late$Elevation<1 )->index_elev_1_late



cv_results_early[[3]]<-list(fold=3, observed=temp_early$Rain_mean[index_elev_1_early], predicted_mean=cv_fit_elevation_below_1_early$summary.fitted.values[index_elev_1_early,"mean"], predicted_sd=cv_fit_elevation_below_1_early$summary.fitted.values[index_elev_1_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_below_1_early$cpo$cpo), na.rm=T),3))
cv_results_late[[3]]<-list(fold=3, observed=temp_late$Rain_mean[index_elev_1_late], predicted_mean=cv_fit_elevation_below_1_late$summary.fitted.values[index_elev_1_late,"mean"], predicted_sd=cv_fit_elevation_below_1_late$summary.fitted.values[index_elev_1_late,"sd"],log_score=-round(mean(log(cv_fit_elevation_below_1_late$cpo$cpo), na.rm=T),3))


cv_df_early<-do.call(rbind, lapply(cv_results_early, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean, predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_early<-cv_df_early %>% mutate( residual=observed-predicted_mean, elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_early<-cv_df_early %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_Score=mean(log_score),MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))

cv_df_late<-do.call(rbind, lapply(cv_results_late, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean,predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_late<-cv_df_late %>% mutate(residual=observed-predicted_mean,elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_late<-cv_df_late %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_Score=mean(log_score),MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))


############## check how well inference is in high mountain areas dry spell ###########
cv_results_early_n.binom<-list()
cv_results_late_n.binom<-list()

temp_early->cv_data_early_n.binom
temp_late->cv_data_late_n.binom

#remove Rain values for elevation above 2000m
cv_data_early_n.binom$Rain_no_length[cv_data_early_n.binom$Elevation>=2]<-NA
cv_data_late_n.binom$Rain_no_length[cv_data_late_n.binom$Elevation>=2]<-NA

obs.data_early <- dplyr::select(cv_data_early_n.binom, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late_n.binom, starts_with('Rain'))


stack.train_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                              effects=list(c(index,index_mean, list(Intercept=1)),
                                           list(covars.data_early)),
                              tag="train_early")

stack.train_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                             effects=list(c(index,index_mean, list(Intercept=1)),
                                          list(covars.data_late)),
                             tag="train_late")

formula4_dry_spell<-Rain_no_length ~ -1 + Intercept + scale(Elevation):Month + Slope + 
  scale(Aspect) + f(spatial_field, model = spde.nonstat) + 
  f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                                                                                                     hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))


cv_fit_elevation_above2_early_n.binom<-inla(formula4_dry_spell, family="nbinomial",
                                            data=inla.stack.data(stack.train_early, spde=spde.nonstat),
                                            control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                            control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                            control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                   openmp.strategy="pardiso.parallel"),
                                            verbose=TRUE)

cv_fit_elevation_above2_late_n.binom<-inla(formula4_dry_spell, family="nbinomial",
                                           data=inla.stack.data(stack.train_late, spde=spde.nonstat),
                                           control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                           control.predictor = list(A=inla.stack.A(stack.train_late),compute=TRUE,link=1),
                                           control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                  openmp.strategy="pardiso.parallel"),
                                           verbose=TRUE)

which(temp_early$Elevation>=2)->index_elev_early
which(temp_late$Elevation>=2)->index_elev_late

cv_results_early_n.binom[[1]]<-list(fold=1, observed=temp_early$Rain_no_length[index_elev_early], predicted_mean=cv_fit_elevation_above2_early_n.binom$summary.fitted.values[index_elev_early,"mean"], predicted_sd=cv_fit_elevation_above2_early_n.binom$summary.fitted.values[index_elev_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_above2_early_n.binom$cpo$cpo), na.rm=T),3))
cv_results_late_n.binom[[1]]<-list(fold=1, observed=temp_late$Rain_no_length[index_elev_late], predicted_mean=cv_fit_elevation_above2_late_n.binom$summary.fitted.values[index_elev_late,"mean"], predicted_sd=cv_fit_elevation_above2_late_n.binom$summary.fitted.values[index_elev_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_above2_late_n.binom$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation between 1000m and 2000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_no_length[cv_data_early$Elevation<2 & cv_data_early$Elevation>=1]<-NA
cv_data_late$Rain_no_length[cv_data_late$Elevation<2 & cv_data_late$Elevation>=1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))

stack.train1_2_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                 effects=list(c(index,index_mean, list(Intercept=1)),
                                              list(covars.data_early)),
                                 tag="train_early")

stack.train1_2_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_late)),
                                tag="train_late")

cv_fit_elevation_between1_2_early_n.binom<-inla(formula4_dry_spell, family="nbinomial",
                                                data=inla.stack.data(stack.train_early, spde=spde.nonstat),
                                                control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                                control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                                control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                       openmp.strategy="pardiso.parallel"),
                                                verbose=TRUE)

cv_fit_elevation_between1_2_late_n.binom<-inla(formula4_dry_spell, family="nbinomial",
                                               data=inla.stack.data(stack.train1_2_late, spde=spde.nonstat),
                                               control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                               control.predictor = list(A=inla.stack.A(stack.train1_2_late),compute=TRUE,link=1),
                                               control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                      openmp.strategy="pardiso.parallel"),
                                               verbose=TRUE)

which(temp_early$Elevation<2 & temp_early$Elevation >=1)->index_elev_1_2_early
which(temp_late$Elevation<2 & temp_late$Elevation >=1)->index_elev_1_2_late

cv_results_early_n.binom[[2]]<-list(fold=2, observed=temp_early$Rain_no_length[index_elev_1_2_early], predicted_mean=cv_fit_elevation_between1_2_early_n.binom$summary.fitted.values[index_elev_1_2_early,"mean"], predicted_sd=cv_fit_elevation_between1_2_early_n.binom$summary.fitted.values[index_elev_1_2_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_early_n.binom$cpo$cpo), na.rm=T),3))
cv_results_late_n.binom[[2]]<-list(fold=2, observed=temp_late$Rain_no_length[index_elev_1_2_late], predicted_mean=cv_fit_elevation_between1_2_late_n.binom$summary.fitted.values[index_elev_1_2_late,"mean"], predicted_sd=cv_fit_elevation_between1_2_late_n.binom$summary.fitted.values[index_elev_1_2_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_late_n.binom$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation below 1000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_no_length[cv_data_early$Elevation<1]<-NA
cv_data_late$Rain_no_length[cv_data_late$Elevation<1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))

stack.train_1_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_early)),
                                tag="train_early")

stack.train_1_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                               effects=list(c(index,index_mean, list(Intercept=1)),
                                            list(covars.data_late)),
                               tag="train_late")

cv_fit_elevation_below_1_early<-inla(formula4_dry_spell, family="nbinomial",
                                     data=inla.stack.data(stack.train_1_early, spde=spde.nonstat),
                                     control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                     control.predictor = list(A=inla.stack.A(stack.train_1_early),compute=TRUE,link=1),
                                     control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                            openmp.strategy="pardiso.parallel"),
                                     verbose=TRUE)


cv_fit_elevation_below_1_late<-inla(formula4_dry_spell, family="nbinomial",
                                    data=inla.stack.data(stack.train_1_late, spde=spde.nonstat),
                                    control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                    control.predictor = list(A=inla.stack.A(stack.train_1_late),compute=TRUE,link=1),
                                    control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                           openmp.strategy="pardiso.parallel"),
                                    verbose=TRUE)
which(temp_early$Elevation<1 )->index_elev_1_early
which(temp_late$Elevation<1 )->index_elev_1_late


cv_results_early_n.binom[[3]]<-list(fold=3, observed=temp_early$Rain_no_length[index_elev_1_early], predicted_mean=cv_fit_elevation_below_1_early$summary.fitted.values[index_elev_1_early,"mean"], predicted_sd=cv_fit_elevation_below_1_early$summary.fitted.values[index_elev_1_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_below_1_early$cpo$cpo), na.rm=T),3))
cv_results_late_n.binom[[3]]<-list(fold=3, observed=temp_late$Rain_no_length[index_elev_1_late], predicted_mean=cv_fit_elevation_below_1_late$summary.fitted.values[index_elev_1_late,"mean"], predicted_sd=cv_fit_elevation_below_1_late$summary.fitted.values[index_elev_1_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_below_1_late$cpo$cpo), na.rm=T),3))


cv_df_early_n.binom<-do.call(rbind, lapply(cv_results_early_n.binom, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean, predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_early_n.binom<-cv_df_early_n.binom %>% mutate(residual=observed-predicted_mean, elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_early_n.binom<-cv_df_early_n.binom %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_score=mean(log_score), MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))

cv_df_late_n.binom<-do.call(rbind, lapply(cv_results_late_n.binom, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean, predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_late_n.binom<-cv_df_late_n.binom %>% mutate(residual=observed-predicted_mean,elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_late_n.binom<-cv_df_late_n.binom %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_score=mean(log_score), MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))

############## check how well inference is in high mountain areas max ###########
map.tail = function(x, interval, inverse = FALSE) {
  if (!inverse) {
    return (interval[1] + (interval[2] - interval[1]) * exp(x)/(1.0 + exp(x)))
  } else {
    return (log((x-interval[1])/(interval[2]-x)))
  }
}

n_early<- length(temp_early$Rain_max)
null.matrix_early<- matrix(nrow=n_early, ncol=0)

n_late<- length(temp_late$Rain_max)
null.matrix_late<- matrix(nrow=n_late, ncol=0)

#models for spread and tail parameter
spread.x_early<-temp_early$Elevation
tail.x_early<-null.matrix_early

spread.x_late<-temp_late$Elevation
tail.x_late<-null.matrix_late

#default prior dist for spread is Gamma with shape and rate parameters=3
hyper.spread = list(initial = 1,
                    fixed=FALSE,
                    prior = "loggamma",
                    param = c(4,1))# c(3, 3)) #instead 3,3

#initial value for tail
tail<- 0.1
tail.interval= c(0, 0.5) #low=0 and high=0.5

tail.intern<- map.tail(tail, tail.interval, inverse=TRUE)

#default prior dist for tail is PC prior with parameter lambda=7, low=0 and high=0.5
hyper.tail <-  list(initial = if (tail == 0.0) -Inf else tail.intern,
                    prior = "pc.gevtail",
                    param = c(7, tail.interval),  #lambda=7
                    fixed= if (tail == 0.0) TRUE else FALSE)

#default hyperparameter specification
hyper.bgev<- list(spread=hyper.spread,
                  tail=hyper.tail)

cv_results_early_max<-list()
cv_results_late_max<-list()

temp_early->cv_data_early_max
temp_late->cv_data_late_max

#remove Rain values for elevation above 2000m
cv_data_early_max$Rain_max[cv_data_early_max$Elevation>=2]<-NA
cv_data_late_max$Rain_max[cv_data_late_max$Elevation>=2]<-NA

obs.data_early <- dplyr::select(cv_data_early_max, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late_max, starts_with('Rain'))


stack.train_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                              effects=list(c(index,index_mean, list(Intercept=1)),
                                           list(covars.data_early)),
                              tag="train_early")

stack.train_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                             effects=list(c(index,index_mean, list(Intercept=1)),
                                          list(covars.data_late)),
                             tag="train_late")

formula4_max_early<- inla.mdata(Rain_max, spread.x_early, tail.x_early) ~ -1 + Intercept + 
  scale(Elevation):Month + Slope + scale(Aspect) + f(spatial_field, 
    model = spde.nonstat) + f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                hyper=list(rho=list(prior="pc.cor1", param=c(0.7,0.8)))))

formula4_max_late<- inla.mdata(Rain_max, spread.x_late, tail.x_late) ~ -1 + Intercept + 
  scale(Elevation):Month + Slope + scale(Aspect) + f(spatial_field, 
                                                     model = spde.nonstat) + f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                                                                                                                                                                                hyper=list(rho=list(prior="pc.cor1", param=c(0.7,0.8)))))


cv_fit_elevation_above2_early_max<-inla(formula4_max_early, family="bgev",
                                             data=inla.stack.data(stack.train_early,spde=spde.nonstat),
                                             control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                             control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                             control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                                    openmp.strategy="pardiso.parallel"),
                                             verbose=TRUE,#inla.mode = "experimental", 
                                             control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                               cmin=0))

cv_fit_elevation_above2_late_max<-inla(formula4_max_late, family="bgev",
                                           data=inla.stack.data(stack.train_late, spde=spde.nonstat),
                                           control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                           control.predictor = list(A=inla.stack.A(stack.train_late),compute=TRUE,link=1),
                                           control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                  openmp.strategy="pardiso.parallel"),
                                           verbose=TRUE,control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                             cmin=0))


which(temp_early$Elevation>=2)->index_elev_early
which(temp_late$Elevation>=2)->index_elev_late

cv_results_early_max[[1]]<-list(fold=1, observed=temp_early$Rain_max[index_elev_early], predicted_mean=cv_fit_elevation_above2_early_max$summary.fitted.values[index_elev_early,"mean"], predicted_sd=cv_fit_elevation_above2_early_max$summary.fitted.values[index_elev_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_above2_early_max$cpo$cpo), na.rm=T),3))
cv_results_late_max[[1]]<-list(fold=1, observed=temp_late$Rain_max[index_elev_late], predicted_mean=cv_fit_elevation_above2_late_max$summary.fitted.values[index_elev_late,"mean"], predicted_sd=cv_fit_elevation_above2_late_max$summary.fitted.values[index_elev_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_above2_late_max$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation between 1000m and 2000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_max[cv_data_early$Elevation<2 & cv_data_early$Elevation>=1]<-NA
cv_data_late$Rain_max[cv_data_late$Elevation<2 & cv_data_late$Elevation>=1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))

stack.train1_2_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                 effects=list(c(index,index_mean, list(Intercept=1)),
                                              list(covars.data_early)),
                                 tag="train_early")

stack.train1_2_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_late)),
                                tag="train_late")

cv_fit_elevation_between1_2_early_max<-inla(formula4_max_early, family="bgev",
                                                data=inla.stack.data(stack.train_early,spde=spde.nonstat),
                                                control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                                control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                                control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                                       openmp.strategy="pardiso.parallel"),
                                                verbose=TRUE,#inla.mode = "experimental", 
                                                control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                                  cmin=0))

cv_fit_elevation_between1_2_late_max<-inla(formula4_max_late, family="bgev",
                                               data=inla.stack.data(stack.train_late, spde=spde.nonstat),
                                               control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                               control.predictor = list(A=inla.stack.A(stack.train_late),compute=TRUE,link=1),
                                               control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                      openmp.strategy="pardiso.parallel"),
                                               verbose=TRUE,control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                                 cmin=0))

which(temp_early$Elevation<2 & temp_early$Elevation >=1)->index_elev_1_2_early
which(temp_late$Elevation<2 & temp_late$Elevation >=1)->index_elev_1_2_late

cv_results_early_max[[2]]<-list(fold=2, observed=temp_early$Rain_no_length[index_elev_1_2_early], predicted_mean=cv_fit_elevation_between1_2_early_max$summary.fitted.values[index_elev_1_2_early,"mean"], predicted_sd=cv_fit_elevation_between1_2_early_max$summary.fitted.values[index_elev_1_2_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_early_max$cpo$cpo), na.rm=T),3))
cv_results_late_max[[2]]<-list(fold=2, observed=temp_late$Rain_no_length[index_elev_1_2_late], predicted_mean=cv_fit_elevation_between1_2_late_max$summary.fitted.values[index_elev_1_2_late,"mean"], predicted_sd=cv_fit_elevation_between1_2_late_max$summary.fitted.values[index_elev_1_2_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_between1_2_late_max$cpo$cpo), na.rm=T),3))

#remove Rain values for elevation below 1000m
temp_early->cv_data_early
temp_late->cv_data_late

cv_data_early$Rain_max[cv_data_early$Elevation<1]<-NA
cv_data_late$Rain_max[cv_data_late$Elevation<1]<-NA

obs.data_early <- dplyr::select(cv_data_early, starts_with('Rain'))
obs.data_late <- dplyr::select(cv_data_late, starts_with('Rain'))

stack.train_1_early<-inla.stack(data=obs.data_early, A=list(A_early,1),
                                effects=list(c(index,index_mean, list(Intercept=1)),
                                             list(covars.data_early)),
                                tag="train_early")

stack.train_1_late<-inla.stack(data=obs.data_late, A=list(A_late,1),
                               effects=list(c(index,index_mean, list(Intercept=1)),
                                            list(covars.data_late)),
                               tag="train_late")

cv_fit_elevation_below_1_early_max<-inla(formula4_max_early, family="bgev",
                                     data=inla.stack.data(stack.train_early,spde=spde.nonstat),
                                     control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                     control.predictor = list(A=inla.stack.A(stack.train_early),compute=TRUE,link=1),
                                     control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                            openmp.strategy="pardiso.parallel"),
                                     verbose=TRUE,#inla.mode = "experimental", 
                                     control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                       cmin=0))


cv_fit_elevation_below_1_late_max<-inla(formula4_max_late, family="bgev",
                                    data=inla.stack.data(stack.train_late, spde=spde.nonstat),
                                    control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                    control.predictor = list(A=inla.stack.A(stack.train_late),compute=TRUE,link=1),
                                    control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                           openmp.strategy="pardiso.parallel"),
                                    verbose=TRUE,control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                      cmin=0))

which(temp_early$Elevation<1 )->index_elev_1_early
which(temp_late$Elevation<1 )->index_elev_1_late


cv_results_early_max[[3]]<-list(fold=3, observed=temp_early$Rain_max[index_elev_1_early], predicted_mean=cv_fit_elevation_below_1_early_max$summary.fitted.values[index_elev_1_early,"mean"], predicted_sd=cv_fit_elevation_below_1_early_max$summary.fitted.values[index_elev_1_early,"sd"], log_score=-round(mean(log(cv_fit_elevation_below_1_early_max$cpo$cpo), na.rm=T),3))
cv_results_late_max[[3]]<-list(fold=3, observed=temp_late$Rain_max[index_elev_1_late], predicted_mean=cv_fit_elevation_below_1_late_max$summary.fitted.values[index_elev_1_late,"mean"], predicted_sd=cv_fit_elevation_below_1_late_max$summary.fitted.values[index_elev_1_late,"sd"], log_score=-round(mean(log(cv_fit_elevation_below_1_late_max$cpo$cpo), na.rm=T),3))


cv_df_early_max<-do.call(rbind, lapply(cv_results_early_max, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean, predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_early_max<-cv_df_early_max %>% mutate(residual=observed-predicted_mean, elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_early_max<-cv_df_early_max %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_score=mean(log_score), MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))

cv_df_late_max<-do.call(rbind, lapply(cv_results_late_max, function(x) {data.frame(test_group=x$fold, observed=x$observed, predicted_mean=x$predicted_mean, predicted_sd=x$predicted_sd, log_score=x$log_score)}))
cv_df_late_max<-cv_df_late_max %>% mutate(residual=observed-predicted_mean,elevation_level=factor(test_group, levels=c("1","2","3"), labels=c("below 1000m" ,"between 1000 and 2000m" , "above 2000m")))
elevation_errors_late_max<-cv_df_late_max %>% group_by(elevation_level) %>% summarise(SD=mean(predicted_sd), Mean=mean(predicted_mean),Log_score=mean(log_score), MAE=mean(abs(residual), na.rm=T), RMSE=sqrt(mean(residual^2, na.rm=T)))



save(cv_df_early_max, cv_df_late_max, cv_df_early, cv_df_late, cv_df_early_n.binom, cv_df_late_n.binomx, file="01_CV_elevation_level_30y.R")
