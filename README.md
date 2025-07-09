
Climate Change in Austria: Precipitation and Dry Spells over 60 years

This repository contains the code used for creating all the results in the paper "Climate Change in Austria: Precipitation and Dry Spells over 60 years", which performs spatio-temporal Bayesian hierarchical modelling of monthly precipitation patterns in Austria. We propose a generalised additive model to investigate changes in precipitation patterns between two climate normal periods (1961-1990 and 1991-2020) over the past 60 years in Austria. Specifically, we analyse three scenarios: monthly mean, monthly maximum daily precipitation sums, and the monthly maximum length of a dry spell.
The respective data-generating processes were identified as the gamma, blended generalised extreme value, and negative binomial distribution.

The Austrian Central Institute for Meteorology and Geodynamics, called GeoSphere (formerly named ZAMG) offers a data hub, https://data.hub.geosphere.at/, for the weather data, that they have collected for any region in Austria and period in the past. This study is based on daily precipitation data (in mm) of all monitoring stations throughout Austria. In particular, the pre-processed data used for our setup can be loaded from the file 01_rain_data_30y.R. 

In order to model the spatial dependencies in the data more realistically, we intended to take the mountainous landscape of Austria into account. Therefore, we downloaded an elevation map of Austria from this website: https://gadm.org/. We used this information, to build e.g. the mesh, the topographic covariates slope and aspec and the prediction data frame. This information can be compiled from file 00_elev_data_30y.R.

In 00_stack.R you can find the implementation of the observation matrices, the non-stationary spde construction, stacks and formulas but also everything you need to implement the extreme value distribuion "bgev" (blended generalised extreme value). 

The file 00_prediction_stack.R contains the prediction stack itself. Then you open 00_inla_pred_30y.R  and run the inla() function with the respective formula, distribution and link function.

The return value functions for the bgev distribution were taken from https://github.com/siliusmv/inlaBGEV and can be found in 00_functions_return_values.R

The plot files apply the respective link function to the model results and then we compute difference maps and show the plot. For the bgev distribution this file is called 00_return_values_bgev.R.

We do Leave-One-Out, Leave-Group-Out and Leave-Elevation range-Out crossvalidation in 00_Crossvalidation.R
