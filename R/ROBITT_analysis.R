#### ROBITT assessment for UK data provided to the BRC by national recording schemes

#devtools::install_github("robboyd/occAssess")

library(occAssess)
library(soaR)
library(BRCmap)
library(rgdal)
library(colorRamps)
library(gridExtra)
library(ggplot2)

## load required data 

# UK shapefile from BRCmap

data(UK)

shp <- UK[UK$REGION == "Great Britain", ]

shp2 <- UK[UK$REGION == "Ireland", ]

# fortify shapefile for use with ggplot2

mapGB <- ggplot2::fortify(shp)

mapIr <- ggplot2::fortify(shp2)

# species occurrence data 

dat <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_dataflow/Research Datasets/Soldierflies/2022/Research dataset/Soldierflies_2022.csv")

names(dat)

## process species occurrence data

# first remove data not resolved to one day

dat <- dat[- which(dat$startdate != dat$enddate)]

# then remove duplicates (in terms of species name, date and monad)

dat <- dat[- which(duplicated(dat[, c("recommended_name", "startdate", "monad")]))]

# drop columns that are not needed for analysis 

dat <- dat[, c("recommended_name", "monad", "startdate")]

# extract coordinates from grid references (needed by occAssess)

coords <- BRCmap::gr_let2num(gridref = dat$monad,
                             centre = TRUE,
                             return_projection = TRUE)

dat <- cbind(dat, coords)

# check if there are any coordinates on the OSNI projection

table(coords$PROJECTION)

# if yes then reproject these onto OSGB

if ("OSNI" %in% coords$PROJECTION) {
  
  GBCRS <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
  
  NICRS <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1 +x_0=200000 +y_0=250000 +ellps=airy +towgs84=482.5,-130.6,564.6,-1.042,-0.214,-0.631,8.15 +units=m +no_defs")
  
  datNI <- dat[which(dat$PROJECTION == "OSNI"), ]
  
  datGB <- dat[which(dat$PROJECTION == "OSGB"), ]

  NIcoords <- datNI[, c("EASTING", "NORTHING")]

  coordinates(NIcoords) <- c("EASTING", "NORTHING")

  proj4string(NIcoords) <- NICRS

  NIcoords <- spTransform(NIcoords, GBCRS)

  datNI[,c("EASTING", "NORTHING")] <- data.frame(NIcoords)
  
  dat <- rbind(datGB, datNI)
  
}

# remove more columns that aren't needed

dat <- dat[, c("recommended_name", "startdate", "EASTING", "NORTHING")]

head(dat)

# create a new column for year (needed by occAssess). Note we'll keep date as it will allow
# us to look specifically at repeat visits later 

dat$year <- substr(dat$startdate, 1, 4)

# create identifier and sptialUncertainty fields (again, needed by occAssess)

dat$identifier <- "soldierflies"

dat$spatialUncertainty <- 1000

#### Now for the occAssess analysis 

# set time periods in which data will be assessed

periods <- as.list(1970:2020)

# assess the geographic representativeness of the data

# map of the density of records

spatCov <- assessSpatialCov(periods = list(1970:2020),
                            dat = dat,
                            species = "recommended_name", 
                            year = "year",
                            identifier = "identifier",
                            x = "EASTING", 
                            y = "NORTHING",
                            spatialUncertainty = "spatialUncertainty",
                            res = 1000,
                            output = "density",
                            logCount = FALSE)

# plot map showing density of records

myCol <- rgb(255, 255, 255, max = 255, alpha = 0, names = "blue50")

spatCov$soldierflies +
  geom_polygon(data = mapGB, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  geom_polygon(data = mapIr, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F)

# map showing the number of years in which each grid cell has been sampled

spatCov2 <- assessSpatialCov(periods = periods,
                            dat = dat,
                            species = "recommended_name", 
                            year = "year",
                            identifier = "identifier",
                            x = "EASTING", 
                            y = "NORTHING",
                            spatialUncertainty = "spatialUncertainty",
                            res = 1000,
                            output = "nPeriods",
                            logCount = FALSE)

spatCov2$soldierflies + 
  geom_polygon(data = mapGB, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  geom_polygon(data = mapIr, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.4, na.value = myCol) 

assessRecordNumber(periods = periods,
                   dat = dat,
                   species = "recommended_name", 
                   year = "year",
                   identifier = "identifier",
                   x = "EASTING", 
                   y = "NORTHING",
                   spatialUncertainty = "spatialUncertainty")

assessSpeciesNumber(periods = periods,
                   dat = dat,
                   species = "recommended_name", 
                   year = "year",
                   identifier = "identifier",
                   x = "EASTING", 
                   y = "NORTHING",
                   spatialUncertainty = "spatialUncertainty")


