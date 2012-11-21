#' @title This is the data to be included into the package
#' @name genus
#' @docType data
#' @author CoForChange project 
#' @references \url{http://www.coforchange.eu/fr/}
#' @keywords data
#' @note The use of this dataset for publication must make reference to the CoForChange project. 
#' @description genus gives the abundance of 27 common tree genera in the tropical moist forest 
#' of the Congo-Basin and 40 geo-referenced environmental variables on one thousand 8 by 8 km plots
#'  (observations). Each plot's data was obtained by aggregating the data measured on a variable 
#'  number of previously sampled 0.5 ha sub-plots. Geo-referenced environmental variables were 
#'  used to describe the physical factors as well as vegetation characteristics. 
#'  14 physical factors were used pertaining the description of topography, geology and rainfall 
#'  of each plot. Vegetation is characterized through 16-days enhanced vegetation index (EVI) data. 
#' @param gen1 to gen27 are the abundance of the 27 common genera on each plot
#' @param altitude is the above-sea level in meters
#' @param pluvio_yr is the mean annual rainfall
#' @param forest forest type classified into seven classes
#' @param pluvio_1 to pluvio_12 are the monthly rainfalls
#' @param geology is a 5-level factor giving the geological substrate
#' @param evi_1 to evi_23 are the means over 10 years of the enhanced vegetation indexes 
#' for each successive 16 days period using  MODIS images
#' @param lon and lat are the longitude/latitude of the plot centers
#' @param surface is the surface of the sampled area contained in each 8 by 8 kms aggregated plot
NULL
