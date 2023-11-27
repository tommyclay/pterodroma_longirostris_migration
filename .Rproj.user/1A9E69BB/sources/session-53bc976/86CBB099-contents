library(GeoLight)
library(GeoLocTools)
library(probGLS)
library(TwGeos)
library(SGAT)
library(MASS)
library(ggplot2)


#### LOAD IN DEPLOYMENT AND RETRIEVAL INFORMATION ####

meta <- read.csv("./Data_inputs/calibration_timings.csv")
# deployment and retrieval
meta$deploy_datetime <- as.POSIXct(paste(as.Date(strptime(meta$deploy_date, format = "%d/%m/%Y")), meta$deploy_time_UTC, sep = " "), tz = "UTC")
meta$retrieve_datetime <- as.POSIXct(paste(as.Date(strptime(meta$retrieve_date, format = "%d/%m/%Y")), meta$retrieve_time_UTC, sep = " "), tz = "UTC")


#### MANUALLY ITERATE THROUGH AND ANNOTATE AND CORRET TWILIGHTS ####


threshold <- 2 # threshold in light
offset <- 21 # offset from UTC

# CHANGE i HERE MANUALLY for each individual, from 1 to 10.
i <- 1

wd <- "./Data_inputs/Raw light/"
files <- dir(wd)

# loading in light file
raw <- readMTlux(paste0(wd, files[i]), skip = 20)
names(raw) <- c("Date", "Light")
raw$Date <- as.POSIXct(raw$Date, tz = "GMT")

# crop by metadata spreadsheet
meta_sub <- subset(meta, Filename == gsub(".lux", "", files[i]),)
raw <- subset(raw, Date >= meta_sub$deploy_datetime & Date < meta_sub$retrieve_datetime ,)
nrow(raw)

# manually check light curves
twl <- preprocessLight(raw, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 20)         # max. light value

head(twl)
twl <- subset(twl, Deleted == FALSE,)
out_wd <- "./Data_outputs/Full twilights/"
write.csv(twl, paste0(out_wd, gsub(".lux", "", files[i]),  ".csv"), row.names = F)
