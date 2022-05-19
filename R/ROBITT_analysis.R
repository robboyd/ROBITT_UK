#### ROBITT assessment for UK data provided to the BRC by national recording schemes

#devtools::install_github("robboyd/occAssess")

library(occAssess)
library(BRCmap)
library(rgdal)
library(colorRamps)
library(gridExtra)
library(ggplot2)
library(rasterVis)

## load required data 

# UK shapefile from BRCmap

BRCmap::data(UK)

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

which(duplicated(dat[, c("year", "monad")]) &
        !duplicated(dat$startdate))

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

dat$identifier <- "all_data"

dat$spatialUncertainty <- 1000

head(dat)

## now create a second dataset with just the repeat visits 

repeats <- dat[which(duplicated(dat[, c("EASTING", "NORTHING", "year")]) &
                       !duplicated(dat[, c("EASTING", "NORTHING", "startdate")])), ]

repeats$identifier <- "repeat_visits" # set identifier to distinguish from the rest

# append to dat for analysis with occAssess

dat <- rbind(dat, repeats)

#### Now for the occAssess analysis 

# set time periods into which the data will be split and assessed

periods <- as.list(1970:2020)

# assess the geographic representativeness of the data

# map grid cells sampled at some point 

spatCov <- assessSpatialCov(periods = periods,
                            dat = dat,
                            species = "recommended_name", 
                            year = "year",
                            identifier = "identifier",
                            x = "EASTING", 
                            y = "NORTHING",
                            spatialUncertainty = "spatialUncertainty",
                            res = 1000,
                            output = "overlap",
                            minPeriods = 1,
                            returnRaster = TRUE)

myCol <- rgb(255, 255, 255, max = 255, alpha = 0, names = "blue50")

gplot(spatCov$rasters) +
  geom_tile(aes(fill = value)) +
  facet_wrap(~variable) +
  geom_polygon(data = mapGB, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  geom_polygon(data = mapIr, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  theme_linedraw() +
  theme(axis.text.x=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank()) +
  labs(fill = "Proportion
       of years
       sampled") +
  labs(x = "",
       y = "") +
  scale_fill_continuous(na.value = myCol) +
  guides(fill = "none")


# map of the density of records

spatCov <- assessSpatialCov(periods = periods,
                            dat = dat,
                            species = "recommended_name", 
                            year = "year",
                            identifier = "identifier",
                            x = "EASTING", 
                            y = "NORTHING",
                            spatialUncertainty = "spatialUncertainty",
                            res = 1000,
                            output = "nPeriods",
                            logCount = FALSE,
                            returnRaster = TRUE)


# plot map showing density of records

mean(raster::getValues(spatCov$rasters$all_data), na.rm=T) / 51

mean(raster::getValues(spatCov$rasters$repeat_visits), na.rm=T) / 51

myCol <- rgb(255, 255, 255, max = 255, alpha = 0, names = "blue50")

rasterVis::gplot(spatCov$rasters) +
  ggplot2::geom_tile(ggplot2::aes(fill = value / 51)) +
  ggplot2::facet_wrap(~variable) +
  ggplot2::geom_polygon(data = mapGB, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  ggplot2::geom_polygon(data = mapIr, ggplot2::aes(x = long, 
                                          y = lat, group = group), colour = "black", 
               fill = myCol, inherit.aes = F) +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.4, na.value = myCol) + 
  ggplot2::theme_linedraw() +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank()) +
  ggplot2::labs(fill = "Proportion
of years
sampled") +
  ggplot2::labs(x = "",
       y = "") 

# nearest neighbour index indicating whether the data are randomly distributed

mask <- raster::raster("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/SDMs/Data/SDMOutputs_Jan_Feb_2021/Bryophytes/Bry_986_LPT_1.asc")

NNI <- assessSpatialBias(dat = dat,
                         periods = periods, 
                         nSamps = 20,
                         degrade = TRUE,
                         mask = mask,
                         species = "recommended_name", 
                         year = "year",
                         identifier = "identifier",
                         x = "EASTING", 
                         y = "NORTHING",
                         spatialUncertainty = "spatialUncertainty",) 


ggplot(data = NNI$data, aes(x = as.numeric(Period) + 1969, y = mean,
                            group = identifier, fill = identifier,
                            colour = identifier)) + 
  geom_line() + 
  theme_linedraw() +
  xlab("Year") +
  ylab("NNI") +
  geom_hline(yintercept = 1, colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3) +
  labs(group = "",
       fill = "",
       colour = "")


lc <- raster::raster("W:/PYWELL_SHARED/Pywell Projects/BRC/Rob Boyd/TSDA/SDMs/Data/targetLandCover2007/selectedCovs/built_up.asc")

dat$built_up <- raster::extract(lc, 
                                cbind(dat$EASTING, dat$NORTHING))

ggplot(data=dat[dat$year>=1970, ], aes(x=year, y=built_up, fill=identifier)) +
  geom_boxplot() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "",
       y = "Proportion of cell classed as built up") +
  labs(fill = "")
