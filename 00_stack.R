library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(dplyr)
library(tidyverse)

load("01_elev_data.R")
load("01_rain_data.R")

#################################################
#in this file:  set up observation matrix, non-stationary spde, stack and formulas
#grouping: according to Month
#bgev: arrangements for bgev distribution are included
#also check (inla.doc("bgev"))
##################################################

#observation matrices for early and late time period
A_early<-inla.spde.make.A(mesh, loc=as.matrix(temp_early[,c('Longitude','Latitude')]), 
                          group = temp_early$Month)

A_late<-inla.spde.make.A(mesh, loc=as.matrix(temp_late[,c('Longitude','Latitude')]), 
                         group = temp_late$Month)

#define non-stationary spde
spde.nonstat<-inla.spde2.matern(mesh, B.tau=cbind(0, 0, 0, 1, values_on_mesh), #first col offset
                                B.kappa=cbind(0, 1, values_on_mesh, 0, 0),
                                theta.prior.mean=rep(0,4),#c(-4,0,4,0),
                                theta.prior.prec=rep(1,4))

#index for mean precipitation
index_mean<- inla.spde.make.index("field", n.spde=spde.nonstat$n.spde, 
                                  n.group = dim(A_early)[2]/spde.nonstat$n.spde)

#oindex for maximum precipitation
index_max<- inla.spde.make.index("field_max", n.spde=spde.nonstat$n.spde,
                                 n.group = dim(A_early)[2]/spde.nonstat$n.spde)

#separate response variable from fixed effects
obs.data_early <- dplyr::select(temp_early, starts_with('Rain'))
covars.data_early <-dplyr::select(temp_early, -starts_with('Rain'))

obs.data_late <- dplyr::select(temp_late, starts_with('Rain'))
covars.data_late <-dplyr::select(temp_late, -starts_with('Rain'))

#stack for mean precipitation in early and late time period 
stack_early<-inla.stack(data=obs.data_early,
                        A=list(A_early,1),
                        effects=list(c(index_mean, list(Intercept=1)),
                                     list(covars.data_early)),
                        tag="mean_early")

stack_late<-inla.stack(data=obs.data_late,
                       A=list(A_late,1),
                       effects=list(c(index_mean, list(Intercept=1)),
                                    list(covars.data_late)),
                       tag="mean_late")

#stack for maximum precipitation in early and late time period 
stack_max_early<-inla.stack(data=obs.data_early,
                            A=list(A_early,1),
                            effects=list(c(index_max, list(Intercept=1)),
                                         list(covars.data_early)),
                            tag="max_early")

stack_max_late<-inla.stack(data=obs.data_late,
                           A=list(A_late,1),
                           effects=list(c(index_max, list(Intercept=1)),
                                        list(covars.data_late)),
                           tag="max_late")


#formulas for mean precipitation, for days without rain, for the lenght of dry spells
formula4_mean<-Rain_mean ~ -1+Intercept+Latitude.cent+Longitude.cent+Elevation+
  f(field,model=spde.nonstat, group = field.group, 
    control.group = list(model="ar1"))

formula4_min<- Rain_no_days~ -1 +Intercept + Latitude.cent+Longitude.cent+Elevation+
  f(field,model=spde.nonstat, group = field.group, 
    control.group = list(model="ar1"))

formula4_dry_spell<-Rain_no_length~ -1 +Intercept + Latitude.cent+Longitude.cent+Elevation+
  f(field,model=spde.nonstat, group = field.group, 
    control.group = list(model="ar1"))

#functions for bgev distribution
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
                    param = c(3, 3))

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

control.bgev<- list(q.location = 0.5,
                    q.spread = 0.8,   #choose =0.8 the higher the more numerically stable
                    # quantile levels for the mixing part
                    q.mix= c(0.05, 0.20),
                    # the Beta(s1, s2) mixing distribution parameters.
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)

formula4_max_early<-inla.mdata(Rain_max, spread.x_early, tail.x_early) ~ -1+Intercept +Latitude.cent+Longitude.cent+Elevation+
  f(field_max, model=spde.nonstat, group = field_max.group, 
    control.group = list(model="ar1"))

formula4_max_late<-inla.mdata(Rain_max, spread.x_late, tail.x_late) ~ -1+Intercept +Latitude.cent+Longitude.cent+Elevation+ 
  f(field_max, model=spde.nonstat, group = field_max.group, 
    control.group = list(model="ar1"))

save(spde.nonstat, 
     stack_early, stack_late, stack_max_early, stack_max_late,
     formula4_mean, formula4_min, formula4_max_early, formula4_max_late, formula4_dry_spell,
     obs.data_early, obs.data_late, covars.data_early, covars.data_late,
     index_mean, index_max,
     hyper.bgev,control.bgev, 
     spread.x_early,tail.x_early,tail.x_late,  spread.x_late, 
      file="01_stack.R")

