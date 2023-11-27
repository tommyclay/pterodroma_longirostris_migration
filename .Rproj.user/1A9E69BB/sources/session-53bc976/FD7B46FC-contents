library(GeoLight)
library(GeoLocTools)
library(probGLS)
library(TwGeos)
library(SGAT)
library(lubridate)
library(MASS)

# This code gets raw light files, crops by start and end calibration timings, identifies twilights manually and uses twilights to get zenith angles. 


#### LOAD IN CALIBRATION TIMINGS ####


cal <- read.csv("./Data_inputs/calibration_timings.csv")
# start of deployment calibration
cal$cal1_st_datetime <- as.POSIXct(paste(as.Date(strptime(cal$cal1_st_date, format = "%d/%m/%Y")), cal$cal1_st_time_UTC, sep = " "), tz = "UTC")
cal$cal1_end_datetime <- as.POSIXct(paste(as.Date(strptime(cal$cal1_end_date, format = "%d/%m/%Y")), cal$cal1_end_time_UTC, sep = " "), tz = "UTC")
# end of deployment calibration
cal$cal2_st_datetime <- as.POSIXct(paste(as.Date(strptime(cal$cal2_st_date, format = "%d/%m/%Y")), cal$cal2_st_time_UTC, sep = " "), tz = "UTC")
cal$cal2_end_datetime <- as.POSIXct(paste(as.Date(strptime(cal$cal2_end_date, format = "%d/%m/%Y")), cal$cal2_end_time_UTC, sep = " "), tz = "UTC")


#### 1. START CALIBRATION ####

# load in light data, crop calibration, annotate twilights and paste out - a manual process

threshold <- 2 # threshold in lux
offset <- 21 # offset from UTC

# CHANGE i HERE MANUALLY for each individual, from 1 to 10.
i <- 1

wd <- "./Data_inputs/Raw light/"
files <- dir(wd)

# loading in light file
raw <- readMTlux(paste0(wd, files[i]), skip = 20)
names(raw) <- c("Date", "Light")
raw$Date <- force_tz(raw$Date, tzone = "UTC")

# crop by metadata spreadsheet
cal_sub <- subset(cal, Filename == gsub(".lux", "", files[i]),)
raw <- subset(raw, Date >= cal_sub$cal1_st_datetime & Date < cal_sub$cal1_end_datetime ,)

# manually check light curves and assign start and end twilights
twl <- preprocessLight(raw, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 20)         # max. light value

out_wd <- "./Data_outputs/Calibration twilights/"
write.csv(twl, paste0(out_wd, gsub(".lux", "", files[i]),  " cal1 twl.csv"), row.names = F)

# LOAD IN GROUNDTRUTH PERIODS TO GET ZENITHS #

lon.calib <- -80.73
lat.calib <- -33.75

# getting zenith and alpha values
zen <- vector(length= length(files))
zen0 <- vector(length=  length(files))
alpha_mean <- vector(length= length(files))
alpha_sd <- vector(length= length(files))

for (i in 1:length(files)) {
  #  load in calibration period 
  cal_twl <- read.csv(paste0(out_wd, gsub(".lux", "", files[i]),  " cal1 twl.csv"))
  cal_twl$Twilight <- as.POSIXct(cal_twl$Twilight,tz = "UTC") 
  # define calibration period
  calib <- thresholdCalibration(cal_twl$Twilight, cal_twl$Rise, lon.calib, lat.calib, method = "gamma")
  zen[i]  <- calib[1]
  zen0[i] <- calib[2]
  alpha_mean[i] <- calib[3]
  alpha_sd[i] <- calib[4]
}   

range(zen)
# 92.28572 93.97971
range(zen0)
# 93.61578 98.72228

# PASTE OUT #

out_df <- data.frame(File = files, ID = gsub(".lux", "", files),
                     zen = zen, zen0 = zen0, alpha_mean = alpha_mean, alpha_sd = alpha_sd)
out_name <- "./Data_outputs/cal1_zenith.csv"
write.csv(out_df, out_name)



#### 2. START CALIBRATION ####

# load in light data, crop calibration, annotate twilights and paste out - a manual process

threshold <- 2 # threshold in lux
offset <- 21 # offset from UTC

# CHANGE i HERE MANUALLY for each individual, from 1 to 10.
i <- 1

wd <- "./Data_inputs/Raw light/"
files <- dir(wd)

# loading in light file
raw <- readMTlux(paste0(wd, files[i]), skip = 20)
names(raw) <- c("Date", "Light")
raw$Date <- force_tz(raw$Date, tzone = "UTC")

# crop by metadata spreadsheet
cal_sub <- subset(cal, Filename == gsub(".lux", "", files[i]),)
raw <- subset(raw, Date >= cal_sub$cal2_st_datetime & Date < cal_sub$cal2_end_datetime ,)

# manually check light curves and assign start and end twilights
twl <- preprocessLight(raw, 
                       threshold = threshold,
                       offset = offset, 
                       lmax = 20)         # max. light value

out_wd <- "./Data_outputs/Calibration twilights/"
write.csv(twl, paste0(out_wd, gsub(".lux", "", files[i]),  " cal2 twl.csv"), row.names = F)

# LOAD IN GROUNDTRUTH PERIODS TO GET ZENITHS #

lon.calib <- -80.73
lat.calib <- -33.75

# getting zenith and alpha values
zen <- vector(length= length(files))
zen0 <- vector(length=  length(files))
alpha_mean <- vector(length= length(files))
alpha_sd <- vector(length= length(files))

for (i in 1:length(files)) {
  #  load in calibration period 
  cal_twl <- read.csv(paste0(out_wd, gsub(".lux", "", files[i]),  " cal2 twl.csv"))
  cal_twl$Twilight <- as.POSIXct(cal_twl$Twilight,tz = "UTC") 
  # define calibration period
  calib <- thresholdCalibration(cal_twl$Twilight, cal_twl$Rise, lon.calib, lat.calib, method = "gamma")
  zen[i]  <- calib[1]
  zen0[i] <- calib[2]
  alpha_mean[i] <- calib[3]
  alpha_sd[i] <- calib[4]
}   

range(zen)
# 94.7362 95.6623
range(zen0)
# 95.10534 104.68340

# PASTE OUT #

out_df <- data.frame(File = files, ID = gsub(".lux", "", files),
                     zen = zen, zen0 = zen0, alpha_mean = alpha_mean, alpha_sd = alpha_sd)
out_name <- "./Data_outputs/cal2_zenith.csv"
write.csv(out_df, out_name)



#### 3. COMBINE TWO CALIBRATION PERIODS AND GET ZENITH ####


lon.calib <- -80.73
lat.calib <- -33.75

zen <- vector(length= length(files))
zen0 <- vector(length=  length(files))
alpha_mean <- vector(length= length(files))
alpha_sd <- vector(length= length(files))

for (i in 1:length(files)) {
  #  load in both calibration periods
  cal_twl1 <-  read.csv(paste0(out_wd, gsub(".lux", "", files[i]),  " cal1 twl.csv"))
  cal_twl1$Twilight <- as.POSIXct(cal_twl1$Twilight,tz = "GMT") 
  cal_twl2 <-  read.csv(paste0(out_wd, gsub(".lux", "", files[i]),  " cal2 twl.csv"))
  cal_twl2$Twilight <- as.POSIXct(cal_twl2$Twilight,tz = "GMT") 
  cal_twl_both <- rbind(cal_twl1,cal_twl2)
  # define calibration period
  calib <- thresholdCalibration(cal_twl_both$Twilight, cal_twl_both$Rise, lon.calib, lat.calib, method = "gamma")
  zen[i]  <- calib[1]
  zen0[i] <- calib[2]
  alpha_mean[i] <- calib[3]
  alpha_sd[i] <- calib[4]
}   

range(zen)
# 93.04323 94.50927
range(zen0)
# 95.10534 104.68340



##### 3C. PASTE OUT ######


out_df <- data.frame(File = files, ID = gsub(".lux", "", files),
                     zen = zen, zen0 = zen0, alpha_mean = alpha_mean, alpha_sd = alpha_sd)
out_name <- "./Data_outputs/cal_both_zenith.csv"
write.csv(out_df, out_name)

