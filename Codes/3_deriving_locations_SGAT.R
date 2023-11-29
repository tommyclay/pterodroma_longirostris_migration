library(GeoLight)
library(GeoLocTools)
library(probGLS)
library(TwGeos)
library(SGAT)
library(MASS)
library(ggplot2)
setupGeolocation()
library(cowplot)


#### 1. LOAD IN TWILIGHTS AND CALCULATE GAPS TO LATER FIX THESE PERIODS ON THE COLONY ####

wd_twl <- "./Data_outputs/Full twilights/"
twilights <- dir(wd_twl)

# for each individual, load in twilight data and calculate gaps (based on 24 hours or more in the dark) to later define as on the colony

twl <- list()
gaps <- list()
for (i in 1:length(twilights)) {
  twl[[i]] <- read.csv(paste0(wd_twl, twilights[i]))
  twl[[i]]$Twilight <- as.POSIXct(twl[[i]]$Twilight, tz = "UTC")
  twl[[i]]$t_diff <- c(as.numeric(difftime(twl[[i]]$Twilight3[2:nrow(twl[[i]])], 
                                      twl[[i]]$Twilight3[1:(nrow(twl[[i]])-1)],
                                      units = "hours")), NA)
  # splitting into dark periods of > 24 hours
  twl[[i]]$Seq <- seq(1, nrow(twl[[i]]), 1)
  gaps[[i]] <- subset(twl[[i]], t_diff >24 & Rise == FALSE,)
}


#### 2. LOAD IN ZENITHS ####

zen <- read.csv("./Data_outputs/cal_both_zenith.csv")
zen_1 <- read.csv("./Data_outputs/cal1_zenith.csv")
zen_2 <- read.csv("./Data_outputs/cal2_zenith.csv")

# for GLS8 D02 B67, the zenith0 value is oddly high (104), so using the mean of all other individuals.
zen$zen0[zen$ID == "GLS8 D02 B67"] <- mean(zen[!(zen$ID ==  "GLS8 D02 B67"),]$zen0)
zen$alpha_mean[zen$ID == "GLS8 D02 B67"] <- mean(zen[!(zen$ID ==  "GLS8 D02 B67"),]$alpha_mean)
zen$alpha_sd[zen$ID == "GLS8 D02 B67"] <- mean(zen[!(zen$ID ==  "GLS8 D02 B67"),]$alpha_sd)
zen$zen[zen$ID == "GLS8 D02 B67"] <- mean(zen[!(zen$ID ==  "GLS8 D02 B67"),]$zen)

range(zen$zen)
# 93.20690 94.50927
range(zen$zen0)
# 95.10534 97.79259


#### 3. PREPARING FUNCTIONS FOR RUNNING SGAT ####

# CREATING LAND/SEA MASK
earthseaMask <- function(xlim, ylim, n = 2, pacific=FALSE) {
  
  if (pacific) { wrld_simpl <- nowrapRecenter(wrld_simpl, avoidGEOS = TRUE)}
  
  # create empty raster with desired resolution
  r = raster(nrows = n * diff(ylim), ncols = n * diff(xlim), xmn = xlim[1],
             xmx = xlim[2], ymn = ylim[1], ymx = ylim[2], crs = proj4string(wrld_simpl))
  
  # create a raster for the stationary period, in this case by giving land a value of 1 and sea NA
  mask = cover(rasterize(elide(wrld_simpl, shift = c(-360, 0)), r, 1, silent = TRUE),
               rasterize(wrld_simpl, r, 1, silent = TRUE), 
               rasterize(elide(wrld_simpl,shift = c(360, 0)), r, 1, silent = TRUE))
  
  xbin = seq(xmin(mask),xmax(mask),length=ncol(mask)+1)
  ybin = seq(ymin(mask),ymax(mask),length=nrow(mask)+1)
  
  function(p) mask[cbind(.bincode(p[,2],ybin),.bincode(p[,1],xbin))]
}

# CREATING MOVEMENT MODEL
beta  <- c(2, 0.1) # from Franklin et al. 2022
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")

# DOWNLOAD WORLD MAP
mapworld <- map_data("world", wrap = c(0, 360), ylim = c(-75,75))



#### 4. RUNNING SGAT AND OUTPUTTING TRACK 10 TIMES PER INDIVIDUAL ####


n_iter <- 10

# SPECIFYING OUTPUT DIRECTORIES
track_dir <- "./Data_outputs/SGAT full tracks/"
plot_mcmc <- "./Plots/MCMC plots/"
map_whole_m <- "./Plots/Whole extent month map/"

# setting lon and lat of colony
lon.calib <- -80.73
lat.calib <- -33.75

# RUNNING FOR EACH INDIVIDUAL 
for (i in 1:length(twilights)) {
  print(paste(gsub(".csv", "", twilights[i])))
  
  # RUNNING INITIAL PATH #
  path <- thresholdPath(twl[[i]]$Twilight, twl[[i]]$Rise, zenith = zen$zen[i], tol=0.01)
  x0 <- path$x  
  z0 <- trackMidpts(x0) # getting midpoints of a path
  
  # DEFINE KNOWN LOCATIONS #
  
  # fixing start and end locations as colony
  fixedx <- rep(F, nrow(x0))
  fixedx[1] <- T # first location estimate
  fixedx[nrow(x0)] <- T # last location estimate
  
  # assigning gap periods at colony
  for (a in 1:nrow(gaps[[i]])) {
    fixedx[gaps[[i]]$Seq[a]] <- T 
    fixedx[gaps[[i]]$Seq[a]+1] <- T 
  }
  # assigning fixed values as colony locations
  x0[fixedx, 1] <- lon.calib
  x0[fixedx, 2] <- lat.calib
  z0 <- trackMidpts(x0) # we need to update the z0 locations
  
  # SETTING LAND MASK PRIOR
  
  xlim <- range(x0[,1]+c(-5,5))
  ylim <- range(x0[,2]+c(-5,5))
  mask <- earthseaMask(xlim, ylim, n = 1, pacific = TRUE)
  log.prior <- function(p) {
    f <- mask(p)
    ifelse(f | is.na(f), -100, -100) # giving these very small probability of occurences of being on land
  }
  
  
  # RUN MODEL 10 TIMES 
  
  for (a in 1:n_iter) {
    
    print(paste(a, n_iter, sep = "..."))
    
    # RUNNING ESTELLE MODEL #
    
    model <- thresholdModel(twilight = twl[[i]]$Twilight,
                            rise = twl[[i]]$Rise,
                            twilight.model = "ModifiedGamma",
                            alpha = c(zen$alpha_mean[i], zen$alpha_sd[i]),
                            beta = beta,
                            logp.x = log.prior, logp.z = log.prior, 
                            x0 = x0,
                            z0 = z0,
                            zenith = zen$zen0[i],
                            fixedx = fixedx)
    # We also need to define the error distribution around each location. We set that using a multivariate normal distribution. Then we can fit the model:
    proposal.x <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0))
    proposal.z <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0))
    
    fit <- estelleMetropolis(model, proposal.x, proposal.z, iters =2000, thin = 20)
    
    # TUNING THE PROPOSALS #
    
    x0 <- chainLast(fit$x)
    z0 <- chainLast(fit$z)
    
    model <- thresholdModel(twilight = twl[[i]]$Twilight,
                            rise = twl[[i]]$Rise,
                            twilight.model = "Gamma",
                            alpha = c(zen$alpha_mean[i], zen$alpha_sd[i]),
                            beta = beta,
                            logp.x = log.prior, logp.z = log.prior, 
                            x0 = x0,
                            z0 = z0,
                            zenith = zen$zen0[i],
                            fixedx = fixedx)
    
    x.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl[[i]]))
    z.proposal <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl[[i]]) - 1)
    #A number of short runs are conducted to tune the proposals. At the end of each run, new proposal distributions are defined based on the dispersion observed in the previous run.
    
    for (k in 1:3) {
      fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                               z0 = chainLast(fit$z), iters = 300, thin = 20)
      
      x.proposal <- mvnorm(chainCov(fit$x), s = 0.2)
      z.proposal <- mvnorm(chainCov(fit$z), s = 0.2)
    }
    #The samples drawn through this process need to be examined to ensure the chain mixes adequately
    # convert before plotting
    dat1 <- as.data.frame(fit$x[[1]][!fixedx, 1, ])
    dat1_2 <- data.frame(x=unlist(dat1))
    names(dat1_2) <- "Lon"
    dat1_2$x <- rep(seq(1,300, 1), each = nrow(dat1))
    dat1_2$Seq <- rep(seq(1,nrow(dat1), 1), 300)
    dat2 <- as.data.frame(fit$x[[1]][!fixedx, 2, ])
    dat2_2 <- data.frame(x=unlist(dat2))
    names(dat2_2) <- "Lat"
    dat2_2$x <- rep(seq(1,300, 1), each = nrow(dat2))
    dat2_2$Seq <- rep(seq(1,nrow(dat2), 1), 300)
    p1 <- ggplot(dat1_2, aes(x, Lon)) + geom_line(aes(group = Seq), size = 0.05, col = "dodgerblue") + xlab("")
    p2 <- ggplot(dat2_2, aes(x, Lat)) + geom_line(aes(group = Seq), size = 0.05, col = "firebrick") + xlab("")
    both <- plot_grid(p1, p2, nrow =2)
    print(both)
    # paste out
    out_name <- paste0(plot_mcmc, gsub(" ", "_", gsub(".csv", "", twilights[i])),"_", "iter_", a,  ".png")
    ggsave(out_name, width = 2000, height = 1500, units = "px", dpi = 300)
    dev.off()
    
    # FINAL RUN #
    
    x.proposal <- mvnorm(chainCov(fit$x), s = 0.25)
    z.proposal <- mvnorm(chainCov(fit$z), s = 0.25)
    
    fit <- estelleMetropolis(model, x.proposal, z.proposal, x0 = chainLast(fit$x), 
                             z0 = chainLast(fit$z), iters = 3000, thin = 20, chains = 4)
    
    # SUMMARIZE  AND PASTE OUT #
    
    sm <- locationSummary(fit$z, time=fit$model$time)
    sm$ID <- gsub(".csv", "", twilights[i])
    sm$Iter <- a
    
    out_name <- paste0(track_dir,  gsub(" ", "_", gsub(".csv", "", twilights[i])),"_", "iter_", a, ".csv")
    write.csv(sm, out_name)
    
    
    # CREATE GENERAL EXTENT PLOT # 
    sm$Lon2 <- sm$Lon.mean
    sm$Lon2[sm$Lon2<0] <- sm$Lon2[sm$Lon2<0] +360
    sm$Month <- as.factor(as.character(substr(sm$Time1, 6, 7)))
    p1 <- ggplot(sm[sm$Lat.mean<90,], aes(Lon2, Lat.mean)) + 
      geom_polygon(data = mapworld, aes(x = long, y = lat, group = group), col = "grey50", fill = "grey70") +
      geom_point(aes(colour = Month)) + geom_path() + 
      theme_bw() + coord_fixed(xlim = c(120, 300), ylim = c(-65, 50))+
      ggtitle(paste(zen$ID[i], a, sep = " "))
    print(p1)
    # paste out
    out_name <- paste0(map_whole_m, gsub(" ", "_", gsub(".csv", "", twilights[i])),"_", "iter_", a, ".png")
    ggsave(out_name, width = 2000, height = 1500, units = "px", dpi = 300)
    dev.off()
  }
}





#### 5. PLOTTING ALL ITERATIONS TOGETHER FOR EACH INDIVIDUAL ####


# load in geolocator files
f_direx <- "./Data_outputs/SGAT full tracks/"
files <- dir(f_direx)

# download world map
mapworld <- map_data("world", wrap = c(0, 360), ylim = c(-75,75))

f_l <-list()
for (i in 1:length(files)) {
  f_l[[i]] <- read.csv(paste0(f_direx, files[i]))
}
all <- do.call(rbind, f_l)
all$ID <- as.factor(as.character(all$ID))
all$Lon2 <- all$Lon.mean
all$Lon2[all$Lon2<0] <- all$Lon2[all$Lon2<0] +360

# plot for each individual
map_whole <- "./Plots/"

for (i in 1:nlevels(all$ID)) {
  print(i)
  id <- subset(all, ID == levels(all$ID)[i],)
  id$Iter <- as.factor(as.character(id$Iter))
  p1 <- ggplot(id[id$Lat.mean<90,], aes(Lon2, Lat.mean)) + 
    geom_polygon(data = mapworld, aes(x = long, y = lat, group = group), col = "grey50", fill = "grey70") +
    geom_point(aes(colour = Iter), size = 0.5) + geom_path(aes(colour = Iter)) + 
    theme_bw() + coord_fixed(xlim = c(120, 300), ylim = c(-65, 60))+
    ggtitle(levels(all$ID)[i])
  print(p1)
  # paste out
  out_name <- paste0(map_whole, gsub(" ", "_", levels(all$ID)[i]), ".png")
  ggsave(out_name, width = 2000, height = 1500, units = "px", dpi = 300)
  dev.off()
}
