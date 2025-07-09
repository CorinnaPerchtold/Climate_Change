library(sf)
library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(raster)
library(inlabru)
library(ggplot2)
library(viridis)
library(viridisLite)
library(terra)
setwd("~/Documents/Prediction_precip/30_year_model")

load("01_rain_data_30y.R")

setwd("/home/guests/corinnap/Documents/Modelling_climate_change/Data_Austria")
#in this file: load in Austria shapefile, create mesh, 
#and compute altitude values on mesh


#### try to keep degrees as longitude latitude
Aut.elev<-raster("AUT_msk_alt.grd")

#res(Aut.elev)   # 0.008333333 0.008333333 degrees

temp_early_sf<-st_as_sf(temp_early, coords=c("Longitude","Latitude"), crs=4326)
temp_late_sf<-st_as_sf(temp_late, coords=c("Longitude","Latitude"), crs=4326)

st_transform(temp_early_sf, crs=32633)->projected_early
st_transform(temp_late_sf, crs=32633)->projected_late

Aut.elev.aggregate<-aggregate(Aut.elev, fact=3)  
#res(Aut.elev.aggregate)   #with fact=3 0.025 0.025 in degrees

dem_proj<-terra::project(rast(Aut.elev.aggregate), "EPSG:32633")

as.data.frame(dem_proj,xy=T)->dem_proj_df
dem_proj<-rast(dem_proj_df, type="xyz")

slope_raster<-terrain(dem_proj, v="slope", unit="degrees", neighbors=8)

#aspect is the compass direction the slope is facing. 0°=N, 90°=E, 180°=S, 270°=W
aspect_raster<-terrain(dem_proj, v="aspect", unit="degrees", neighbors=8)


temp_early_vec<-vect(projected_early)
temp_late_vec<-vect(projected_late)

Slope_early<-terra::extract(slope_raster, temp_early_vec)[2]
Slope_late<-terra::extract(slope_raster, temp_late_vec)[2]

Aspect_early<-terra::extract(aspect_raster, temp_early_vec)[2]
Aspect_late<-terra::extract(aspect_raster, temp_late_vec)[2]

temp_early$Slope<-Slope_early$slope
temp_late$Slope<-Slope_late$slope

temp_early$Aspect<-Aspect_early$aspect
temp_late$Aspect<-Aspect_late$aspect

temp_early<-temp_early[,c(1:13,16:17)]
temp_late<-temp_late[,c(1:13,16:17)]


elev_masked<-mask(dem_proj, slope_raster)   #takes elev values to dimensions of slope_raster
slope_raster<-mask(slope_raster, elev_masked)  #takes slope values to dimensions of elev_ raster
aspect_raster<-mask(aspect_raster, elev_masked)

crs(elev_masked)<-"EPSG:32633"
project(elev_masked, "EPSG:4326")->elev_masked_long_lat

crs(slope_raster)<-"EPSG:32633"
crs(aspect_raster)<-"EPSG:32633"

project(slope_raster, "EPSG:4326")->slope_raster_long_lat
project(aspect_raster, "EPSG:4326")->aspect_raster_long_lat

#################create prediction data frame ###################

prediction_df<-as.data.frame(elev_masked_long_lat, xy=T)

prediction_df$Slope<-as.data.frame(slope_raster_long_lat)$slope
prediction_df$Aspect<-as.data.frame(aspect_raster_long_lat)$aspect

colnames(prediction_df)<-c("Longitude","Latitude", "Elevation","Slope","Aspect")

temp_early<-temp_early %>% mutate(Time=(Year-min(Year))*12+Month)
temp_early<- temp_early %>% arrange(Time, Station) 

temp_late<-temp_late %>% mutate(Time=(Year-min(Year))*12+Month)
temp_late<- temp_late %>% arrange(Time, Station) 

dp<-sapply(prediction_df,rep.int,times=max(temp_early$Time))
dp<-as.data.frame(dp)

#add time and month to prediction data frame
dp$Time<-rep(1:max(temp_early$Time), each=length(prediction_df$Longitude))
dp$Month<-rep(1:max(temp_early$Month), each=length(prediction_df$Longitude))

#set the amount of days per month
dp$Days_per_month[dp$Month==1| dp$Month ==3 | dp$Month==5 |dp$Month==7|dp$Month==8|dp$Month==10 |dp$Month==12]<-31
dp$Days_per_month[dp$Month==2]<-28
dp$Days_per_month[dp$Month==4|dp$Month==6 |dp$Month==9| dp$Month==11]<-30

dp$Days_per_month<-as.integer(dp$Days_per_month)



#################create mesnh ###################

shapefile<-read_sf("gadm41_AUT_1.shp",)
Austria<-st_geometry(shapefile)

#unify Austria to drop inner state boarders
m<-as(st_union(Austria), "Spatial")


#define boundary of Upper Austria
boundary<-inla.sp2segment(m)

coords.stations<-temp_all[c('Longitude','Latitude')]

#define mesh
mesh<-fm_mesh_2d(loc.domain=coords.stations, boundary = boundary,
                 max.edge = c(0.1,0.35),cutoff=0.1, offset=c(1,0.5))

coords.mesh <- mesh$loc[,1:2]

Aut.elev<-rasterToPoints(Aut.elev)
Aut.elev<-as.data.frame(Aut.elev)
colnames(Aut.elev)<-c("Longitude", "Latitude", "Elevation")

coords.elev<-cbind(Aut.elev$Longitude, Aut.elev$Latitude)

#create projector matrix for elev on mesh points
A.elev<-inla.spde.make.A(mesh, loc=coords.elev)

#compute weighted averages of elevation, with weights equal to
#the "hat" basis functions on the mesh.
values_on_mesh <- (t(A.elev) %*% Aut.elev$Elevation) / colSums(A.elev)

#extrapolate missing values outside boundary
dat <- SpatialPointsDataFrame(coords.mesh,
                              data = data.frame(values = as.vector(values_on_mesh)))

values_on_mesh <- inlabru::bru_fill_missing(
  data = dat[!is.na(dat$values),],
  where = as.data.frame(coords.mesh),
  as.vector(values_on_mesh),
  layer = "values"
)


setwd("/home/guests/corinnap/Documents/Prediction_precip/30_year_model")
save(temp_early, temp_late, Austria,dp, m, mesh, values_on_mesh,Aut.elev,Aut.elev.aggregate,boundary, coords.mesh, file = "01_elev_data_30y_long_lat.R")


