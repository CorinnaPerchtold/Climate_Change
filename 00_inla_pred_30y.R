library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(raster)

setwd("~/Documents/Prediction_precip/30_year_model")

load("01_elev_data_30y_long_lat.R")

#for stationary model
load("01_stack_30y_long_lat.R")
load("01_prediction_stack_30y_long_lat.R")

####### results for mean rain ############

formula4_mean<-Rain_mean~ -1+Intercept+scale(Elevation):Month+Slope+scale(Aspect)+
  #f(spatial_field, model=spde.nonstat)+
  f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
                                                                                                     hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))

result4_early_pred_30y_wt.prior<-inla(formula4_mean, family="gamma",
                                                       data=inla.stack.data(stack_early_pred, spde=spde.nonstat),
                                                       control.family = list(link='log'), 
                                                       control.predictor = list(A=inla.stack.A(stack_early_pred),compute=TRUE,link=1),
                                                       control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                      openmp.strategy="pardiso.parallel"),verbose=TRUE)#, inla.mode = "experimental")

result4_late_pred_30y_wt_prior<-inla(formula4_mean, family="gamma",
                                      data=inla.stack.data(stack_late_pred, spde=spde.nonstat),
                                      control.family = list(link='log'),
                                      control.predictor = list(A=inla.stack.A(stack_late_pred), compute=TRUE,link=1),
                                      control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                      openmp.strategy="pardiso.parallel"),
                                      verbose=TRUE)

####### results for max rain ############


formula4_max_early<- inla.mdata(Rain_max, spread.x_early, tail.x_early) ~ -1 + Intercept + 
  scale(Elevation):Month + Slope + scale(Aspect) + f(spatial_field, 
     model = spde.nonstat) + f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
         hyper=list(rho=list(prior="pc.cor1", param=c(0.7,0.8)))))

formula4_max_late<- inla.mdata(Rain_max, spread.x_late, tail.x_late) ~ -1 + Intercept + 
  scale(Elevation):Month + Slope + scale(Aspect) + f(spatial_field, 
     model = spde.nonstat) + f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",    
       hyper=list(rho=list(prior="pc.cor1", param=c(0.7,0.8)))))

result4_max_early.pred_wt_prior<-inla(formula4_max_early, family="bgev",
                                      data=inla.stack.data(stack_max_early_pred,spde=spde.nonstat),
                                      control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                      control.predictor = list(A=inla.stack.A(stack_max_early_pred),compute=TRUE,link=1),
                                      control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                             openmp.strategy="pardiso.parallel"),
                                      verbose=TRUE,#inla.mode = "experimental", 
                                      control.inla = list( #int.strategy="eb", strategy="simplified.laplace", 
                                        cmin=0))#1e-6


result4_max_late.pred_wt_prior<-inla(formula4_max_late, family="bgev", 
                                     data=inla.stack.data(stack_max_late_pred,spde=spde.nonstat),
                                     control.family = list(hyper=hyper.bgev, control.bgev=control.bgev),
                                     control.predictor = list(A=inla.stack.A(stack_max_late_pred),compute=TRUE,link=1),
                                     control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE,config=TRUE, return.marginals.predictor=FALSE,
                                                            openmp.strategy="pardiso.parallel"),
                                     verbose=TRUE,
                                     #inla.mode = "experimental")#,
                                     control.inla =list(cmin=0))


####### results for dry spell ############

E_early<-rep(1, length(Ntrials_early))
E_late<-rep(1, length(Ntrials_late))



formula4_dry_spell<-Rain_no_length ~ -1 + Intercept + scale(Elevation):Month + Slope + 
  scale(Aspect) + f(spatial_field, model = spde.nonstat) + 
  f(spatio_temporal_field, model=spde.nonstat, group=spatio_temporal_field.group, control.group=list(model="ar1",
             hyper=list(rho=list(prior="pc.cor1", param=c(0.9,0.9)))))

## stationary ##
result4_n.binom_early_pred_30y_wt_prior<-inla(formula4_dry_spell,  family="nbinomial",E=E_early,
                                              data=inla.stack.data(stack_early_pred, spde=spde.nonstat),
                                              control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                              control.predictor = list(A=inla.stack.A(stack_early_pred), compute=TRUE,link=1),
                                              control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                     openmp.strategy="pardiso"),
                                              verbose=TRUE)

result4_n.binom_late_pred_30y_wt_prior<-inla(formula4_dry_spell,  family="nbinomial", E=E_late,
                                             data=inla.stack.data(stack_late_pred, spde=spde.nonstat),
                                             control.family=list(link='log', variant=1, hyper=list(theta=list(prior="loggamma", param=c(1,0.01)))),
                                             control.predictor = list(A=inla.stack.A(stack_late_pred), compute=TRUE,link=1),
                                             control.compute = list(dic=TRUE,cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=FALSE,
                                                                    openmp.strategy="pardiso"),
                                             verbose=TRUE)#, inla.mode = "experimental")

