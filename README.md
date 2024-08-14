The Austrian Central Institute for Meteorology and Geodynamics, called GeoSphere (formerly named ZAMG) offers a data hub, https://data.hub.geosphere.at/, for the weather data,
that they have collected for any region in Austria and period in the past. This study is based on daily precipitation data (in mm) of all monitoring stations
throughout Austria in the time periods 1973-1982 and 2013-2022. In particular, the pre-processed data used for our setup can be loaded from the file 01_rain_data.R. 

In order to model the spatial dependencies in the data more realistically, we intended to take the mountainous landscape of Austria into account. Therefore, we downloaded an elevation
map of Austria from this website: https://gadm.org/. We used this information, to build e.g. the mesh. This information can be loaded from file 01_elev_data.R.

In 00_stack.R you can find the implementation of the observation matrices, the non-stationary spde construction, stacks and formulas but also everything you need to implement the
extreme value distribuion "bgev" (blended generalised extreme value). 

The file 01_prediction_stack.R contains the prediction grid based on the elevation map and the prediction stack itself. Then you call the inla() function with the respective 
formula, distribution (mean precipitation -> gamma distribution, maximum precipitation -> bgev distribution, maximum length of a dry spell -> negative binomial distribution, 
days without precipitation -> binomial distribution). 

The return value functions for the bgev distribution was taken from https://github.com/siliusmv/inlaBGEV
