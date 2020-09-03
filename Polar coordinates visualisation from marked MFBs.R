# script to load from FIJI roi the WT1 positive and negative, and cv, for use in polar co-ordinate determination
# needs to have a field with the central vein in the middle

# Tim Kendall, University of Edinburgh
# Tim.Kendall@ed.ac.uk

# August 2020

# load packages
require(RImageJROI)
require(ggplot2)
require(spatstat)
require(rJava)
require(rChoiceDialogs)
require(extrafont)
require(svglite)

# create function (polar_theta) to get polar from cartesian
polar_theta <- function(x,y) { 
  z <- x + 1i * y 
  res <- 90 - Arg(z) / pi * 180 
  res %% 360 
}

w <- owin (c(0,1392), c(0,1040)) # needs to set to get appropriate centroid

# function to check is ves of class ijroi is clockwise after read.ijroi
clockwise <- function(point_list) {
  point_temp <- rbind(point_list, point_list[1,])
  total = rep(0,nrow(point_temp)-1)
  for (i in 1:nrow(point_temp)-1) {
    total[i] <- (point_temp[i+1,1] - point_temp[i,1])*(point_temp[i+1,2] + point_temp[i,2])
    }
  sum(total)
}

##################################################################################################################
# for one image

# read roi files
neg <- read.ijroi('n.roi')
pos <- read.ijroi('p.roi')
ves <- read.ijroi('v.roi')

# ves needs converting to spatstat form
ves <- ij2spatstat(ves)

# convert ves to psp
ves_psp <- as.psp(ves, window = w)

# convert psp (vessel) to owin, then find centre
owin_ves_psp <- as.owin(ves_psp)
centroid_ves_psp <- centroid.owin(owin_ves_psp)

# put cordinates of each cell in single df
coordinates_mfb <- rbind(neg$coords, pos$coords)
colnames(coordinates_mfb) <- c('x', 'y')
coordinates_mfb <- as.data.frame(coordinates_mfb)

# For each point, n, set to xn,yn - xcentroid_ves_psp,ycentroid_ves_psp = (xn-xC, yn-yC); new x,y
coordinates_mfb$x_centre <- coordinates_mfb$x - centroid_ves_psp$x
coordinates_mfb$y_centre <- coordinates_mfb$y - centroid_ves_psp$y

# use polar_theta to generate theta with respect to centroid_ves_psp
coordinates_mfb$theta <- polar_theta(coordinates_mfb$x_centre, coordinates_mfb$y_centre)

# find the peak density of theta
density_theta <- density(coordinates_mfb$theta)
Peak_theta <- density_theta$x[which.max(density_theta$y)]

# use this peak theta to adjust all thetas, setting peak to '90'
Theta_adjust <- (360-Peak_theta) + 90
coordinates_mfb$theta_adjusted <- coordinates_mfb$theta + Theta_adjust

# for those now over 360, convert back to 0-360 range
coordinates_mfb$thetan_adjusted_corrected <- coordinates_mfb$theta_adjusted

for (i in 1:nrow(coordinates_mfb)) {
  if (coordinates_mfb$theta_adjusted[i] > 360) {
      coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i] - 360
  }
  else {
    coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i]
  }
}

# test visualisation of one image

# histogram of thetas
ggplot(coordinates_mfb, aes(x=theta_adjusted)) + 
  geom_histogram(breaks = seq(0, 720, 5), colour = "grey")

#histogram of '90' adjusted thetas
ggplot(coordinates_mfb, aes(x=thetan_adjusted_corrected)) + 
  geom_histogram(breaks = seq(0, 360, 5), colour = "grey")

# circular plots of histogram with density
ggplot(coordinates_mfb, aes(x=thetan_adjusted_corrected)) + geom_histogram(aes(y=..density..), breaks = seq(0, 360, 5), colour = "grey") +
  geom_density(alpha=.2, fill='#57B2C7') + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Density") + ggtitle("Events by direction") +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30))

# circular plot without density
ggplot(coordinates_mfb, aes(x = thetan_adjusted_corrected)) + geom_histogram(breaks = seq(0, 360, 10), colour = "grey") + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Events by direction") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30))

# density of re-referenced to 90 v original
ggplot(coordinates_mfb, aes(x=thetan_adjusted_corrected)) +  
  geom_density(alpha=.2, fill='#57B2C7')

# the density at peak then gets split, so need to plot original density
ggplot(coordinates_mfb, aes(x=theta_adjusted)) +  geom_density(alpha=.2, fill='#57B2C7')

###############################################################################################################################

# as a script for a folder of subfolders with p, n, v.roi

# select folder containing rois
initial.dir<-getwd() # current working directory
data_dir <- jchoose.dir(caption = "Select individual animal directory containing p,n,v ROIs")

# get list of separate image folder
image_dirs <- list.dirs(path = data_dir)
image_dirs <- image_dirs[-1] # remove the parent directory

# function to get polar coordinates from the rois, input is a folder
roi_to_polar <- function(x){
  setwd(x)
  
  # get name of each image based on image directory
  image_number <- gsub("(.*/)(.*$)", "\\2", x)
  
  # read roi files
  neg <- read.ijroi('n.roi')
  pos <- read.ijroi('p.roi')
  ves <- read.ijroi('v.roi')
  
  # check if ves is clockwise or anticlockwise and convert if needed
  if (clockwise(ves$coords) < 0) {
    ves$coords <- apply(ves$coords, 2, rev)
  }
  
  # ves needs converting to spatstat form
  ves <- ij2spatstat(ves)
  
  # convert ves to psp
  ves_psp <- as.psp(ves, window = w)
  
  # convert psp (vessel) to owin, then find centre
  owin_ves_psp <- as.owin(ves_psp)
  centroid_ves_psp <- centroid.owin(owin_ves_psp)
  
  # put cordinates of each cell in single df
  coordinates_mfb <- rbind(neg$coords, pos$coords)
  colnames(coordinates_mfb) <- c('x', 'y')
  coordinates_mfb <- as.data.frame(coordinates_mfb)
  
  # For each point, n, set to xn,yn - xcentroid_ves_psp,ycentroid_ves_psp = (xn-xC, yn-yC); new x,y
  coordinates_mfb$x_centre <- coordinates_mfb$x - centroid_ves_psp$x
  coordinates_mfb$y_centre <- coordinates_mfb$y - centroid_ves_psp$y
  
  # use polar_theta to generate theta with respect to centroid_ves_psp
  coordinates_mfb$theta <- polar_theta(coordinates_mfb$x_centre, coordinates_mfb$y_centre)
  
  # find the peak density of theta
  density_theta <- density(coordinates_mfb$theta)
  Peak_theta <- density_theta$x[which.max(density_theta$y)]
  
  # use this peak theta to adjust all thetas, setting peak to '90'
  Theta_adjust <- (360-Peak_theta) + 90
  coordinates_mfb$theta_adjusted <- coordinates_mfb$theta + Theta_adjust
  
  # for those now over 360, convert back to 0-360 range
  coordinates_mfb$thetan_adjusted_corrected <- coordinates_mfb$theta_adjusted
  
  for (i in 1:nrow(coordinates_mfb)) {
    if (coordinates_mfb$theta_adjusted[i] > 360) {
      coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i] - 360
    }
    else {
      coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i]
    }
  }
  
  # add image name as factor to data frame
  coordinates_mfb$image_number <- as.factor(image_number)
  
  return(coordinates_mfb)
  
}

# package all commands in one function
animal_images_polar <- function(){
  # select folder containing rois
  initial.dir<-getwd() # current working directory
  data_dir <- jchoose.dir(caption = "Select individual animal directory containing p,n,v ROIs")
  
  # get list of separate image folder
  image_dirs <- list.dirs(path = data_dir)
  image_dirs <- image_dirs[-1] # remove the parent directory
  
  # create data frame
  df <- lapply(image_dirs, roi_to_polar)
  
  return(df)
  
  setwd(initial.dir)
}


# run on an animal specific parent folder
WT3_1_allMFBs <- animal_images_polar()
WT3_2_allMFBs <- animal_images_polar()
WT3_3_allMFBs <- animal_images_polar()
WT3_4_allMFBs <- animal_images_polar()
WT3_5_allMFBs <- animal_images_polar()
WT3_6_allMFBs <- animal_images_polar()

############################################################################################################
# plots

# function to create plots from each image entry in list
image_polar_plots <- function(image_z){
  # density of re-referenced to 90 v original
  density_plot <- ggplot(image_z, aes(x=thetan_adjusted_corrected)) +  
    geom_density(alpha=.2, fill='#57B2C7')
  
  # circular plot without density
  circular_hist <- ggplot(image_z, aes(x = thetan_adjusted_corrected)) + geom_histogram(breaks = seq(0, 360, 10), colour = "grey") + coord_polar(start = 0) + theme_minimal() + 
    scale_fill_brewer() + ylab("Count") + ggtitle("Events by direction") + 
    scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30))
  
  # get image name
  name <- as.character(image_z$image_number[1])
  
  # save plots
  density_name <- paste("density_", name, ".svg", sep="")
  hist_name <- paste("hist_", name, ".svg", sep="")
  
  ggsave(filename = density_name, plot = density_plot, device = "svg")
  ggsave(filename = hist_name, plot = circular_hist, device = "svg")
  
  } 

# wrapper function to create plots
plot_creator <- function(y){
  require(svglite)
  require(ggplot2)
  
  initial.dir<-getwd() # current working directory
  data_dir <- jchoose.dir(default=initial.dir, caption = "Select directory to save plots to")
  setwd(data_dir)
  
  for (i in 1:length(y)){
    df1 <- as.data.frame(y[[i]])
    image_polar_plots(df1)
  }
  
  setwd(initial.dir)
}

# then run plot_creator(output_of_anima_images_polar_object)

# example plots
svglite("WT3_1_10 density.svg")
ggplot(WT3_1_allMFBs[[2]], aes(x=thetan_adjusted_corrected)) + theme_minimal() + geom_density(alpha=.2, fill='#57B2C7') +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
dev.off()

svglite("WT3_1_10 histogram.svg")
ggplot(WT3_1_allMFBs[[2]], aes(x = thetan_adjusted_corrected)) + geom_histogram(breaks = seq(0, 360, 10), colour = "grey") + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Events by direction") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30)) +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
dev.off()

svglite("WT3_2_5_new density.svg")

ggplot(WT3_2_allMFBs_new[[6]], aes(x=thetan_adjusted_corrected)) + theme_minimal() +
  geom_density(alpha=.2, fill='#57B2C7') +
  xlab("Aligned MFB angle (degrees)") +
  ylab("MFB population density") +
  theme (legend.key = element_blank(),
         text=element_text(family="Calibri"),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))

dev.off()

svglite("WT3_2_5_new histogram.svg")

ggplot(WT3_2_allMFBs_new[[6]], aes(x = thetan_adjusted_corrected)) + geom_histogram(breaks = seq(0, 360, 10), colour = "grey") + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Count/angle") +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30)) +
  theme (legend.key = element_blank(),
         text=element_text(family="Calibri"),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))

dev.off()

# polar histogram prior to correction

svglite("WT3_2_5_uncorrected histogram.svg")

ggplot(WT3_2_allMFBs_new[[6]], aes(x = theta)) + geom_histogram(breaks = seq(0, 360, 10), colour = "grey") + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Count/angle") +
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30)) +
  theme (legend.key = element_blank(),
         text=element_text(family="Calibri"),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))

# example of plotted points

w <- owin (c(0,1392), c(0,1040)) # set window

neg <- read.ijroi('n.roi')
pos <- read.ijroi('p.roi')
ves <- read.ijroi('v.roi')
ves <- ij2spatstat(ves)
#ves_psp <- as.psp(ves, window = w)
p_spat <- ij2spatstat(pos, window=w)
n_spat <- ij2spatstat(neg, window=w)
marks(p_spat) <- NULL
marks(n_spat) <- NULL
all_spat <- superimpose(p_spat, n_spat)

all_spat <- rescale(all_spat, 1/0.323) # rescale into microns
#ves_psp <- rescale(ves_psp, 1/0.323) # rescale into microns
ves <- rescale(ves, 1/0.323) # rescale into microns

centroid_ves_psp <- centroid.owin(ves)

all_spat_plot <- as.data.frame(all_spat)
ves_psp_plot <- as.data.frame(ves_psp)

svglite('WT3_2_00005_centroid.svg')

ggplot() + geom_point(data=all_spat_plot, aes(x,y)) +
  geom_segment(data=ves_psp_plot, mapping=aes(x=x0, y=y0, xend=x1, yend=y1), arrow=NULL, size=1) +
  theme_bw() + coord_fixed(ratio = 1) +
  geom_hline(yintercept = centroid_ves_psp$y) +
  geom_vline(xintercept = centroid_ves_psp$x) +
  xlab(expression(paste("Distance (",mu, "m)"))) +
  ylab(expression(paste("Distance (",mu, "m)"))) +
  scale_x_continuous(expand=c(0,0), limits = c(0,1392/(1/0.323))) +
  scale_y_continuous(expand=c(0,0), limits = c(0,1040/(1/0.323))) +
  theme(text=element_text(family="Calibri"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

dev.off()
  




#########################################################################################################################

# to plot on a single plot

big_WT3_1 <- do.call(rbind,WT3_1_allMFBs) # combine into a single df
big_WT3_2 <- do.call(rbind,WT3_2_allMFBs)
big_WT3_3 <- do.call(rbind,WT3_3_allMFBs)
big_WT3_4 <- do.call(rbind,WT3_4_allMFBs)
big_WT3_5 <- do.call(rbind,WT3_5_allMFBs)
big_WT3_6 <- do.call(rbind,WT3_6_allMFBs)

# add identifier for animal as we will combine into a single df
big_WT3_1$animal <- as.factor('1')
big_WT3_2$animal <- as.factor('2')
big_WT3_3$animal <- as.factor('3')
big_WT3_4$animal <- as.factor('4')
big_WT3_5$animal <- as.factor('5')
big_WT3_6$animal <- as.factor('6')

# combine all into a single df
All_CCl4_animals <- rbind(big_WT3_1, big_WT3_2, big_WT3_3, big_WT3_4, big_WT3_5, big_WT3_6)

# for all shown separately
ggplot(big_WT3_1, aes(x=thetan_adjusted_corrected, fill=image_number)) +  
  geom_density(alpha=.2)

ggplot(big_WT3_1, aes(x = thetan_adjusted_corrected, fill=image_number)) + geom_histogram(breaks = seq(0, 360, 10)) + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer(palette = 'Spectral') + ylab("Count") + ggtitle("Events by direction") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30))

# for combined
ggplot(big_WT3_1, aes(x=thetan_adjusted_corrected)) + geom_density(alpha=.2, fill='#57B2C7')

ggplot(big_WT3_1, aes(x = thetan_adjusted_corrected)) + geom_histogram(breaks = seq(0, 360, 10)) + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer() + ylab("Count") + ggtitle("Events by direction") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30))

# for the summary of each animal together
svglite("Density polar_all.svg")
ggplot(All_CCl4_animals, aes(x=thetan_adjusted_corrected, fill=animal)) + geom_density(alpha=.2) +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
dev.off()

svglite("Polar histogram_all.svg")
ggplot(All_CCl4_animals, aes(x = thetan_adjusted_corrected, fill=animal)) + geom_histogram(breaks = seq(0, 360, 10)) + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer(palette = 'Spectral') + ylab("Count") + ggtitle("Events by direction") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30)) +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))
dev.off()

# to access the data from all central veins of all animals
all_hist <- ggplot(All_CCl4_animals, aes(x=thetan_adjusted_corrected)) + 
  geom_histogram(breaks = seq(0, 360, 5), alpha=.5) + theme_minimal() + scale_fill_brewer()
count_all <- print(all_hist)
raw_count_all_CCl4 <- count_all[["data"]][[1]]
  
# access count data used for histograms
overlay_hist <- ggplot(big_WT3_1, aes(x=thetan_adjusted_corrected)) + 
  geom_histogram(breaks = seq(0, 360, 5), alpha=.5) + theme_minimal() + scale_fill_brewer()
count_hist <- print(overlay_hist)
raw_count <- count_hist[["data"]][[1]]

per <- 180 # estimated period based on density plots
# fit linear model
reslm <- lm(raw_count$density ~ sin(2*pi/per*raw_count$x)+cos(2*pi/per*raw_count$x))
summary(reslm)
# plot the fitted model
rg <- diff(range(raw_count$density))
plot(raw_count$density~raw_count$x,ylim=c(min(raw_count$density)-0.1*rg,max(raw_count$density)+0.1*rg))
lines(fitted(reslm)~raw_count$x,col=4,lty=2)
cor(raw_count$density, predict(reslm))

# for the total data
reslm_all <- lm(raw_count_all_CCl4$density ~ sin(2*pi/per*raw_count_all_CCl4$x)+cos(2*pi/per*raw_count_all_CCl4$x))
summary(reslm_all)
rg_all <- diff(range(raw_count_all_CCl4$density))
plot(raw_count_all_CCl4$density~raw_count_all_CCl4$x,ylim=c(min(raw_count_all_CCl4$density)-0.1*rg,max(raw_count_all_CCl4$density)+0.1*rg))
lines(fitted(reslm_all)~raw_count_all_CCl4$x,col=4,lty=2)

# in ggplot
svglite("All_density_fitted.svg")
ggplot(raw_count_all_CCl4, aes(x=x, y=density)) + geom_point() + ylim(0, max(raw_count_all_CCl4$density)) + theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  geom_line(aes(x=x, y=fitted(reslm_all))) 
dev.off()

# nonlinear least squares model, less good
nls_densityWT3_1 <- nls(raw_count$density ~ sin(2*pi/per*raw_count$x)+cos(2*pi/per*raw_count$x), data = raw_count, start = list(per=180))
summary(nls_densityWT3_1)
cor(raw_count$density, predict(nls_densityWT3_1))

##################################################################

# with new spatstat handling of as.psp, owin change

# function to get polar coordinates from the rois, input is a folder
roi_to_polar_new <- function(x){
  setwd(x)
  
  # get name of each image based on image directory
  image_number <- gsub("(.*/)(.*$)", "\\2", x)
  
  # read roi files
  neg <- read.ijroi('n.roi')
  pos <- read.ijroi('p.roi')
  ves <- read.ijroi('v.roi')
  
  # check if ves is clockwise or anticlockwise and convert if needed
  if (clockwise(ves$coords) < 0) {
    ves$coords <- apply(ves$coords, 2, rev)
  }
  
  # ves needs converting to spatstat form
  ves_psp <- ij2spatstat(ves)
  
  # convert ves to psp - deprecated
  # ves_psp <- as.psp(ves, window = w)
  
  # convert psp (vessel) to owin, then find centre
  # owin_ves_psp <- as.owin(ves_psp) - deprecated
  centroid_ves_psp <- centroid.owin(ves_psp)
  
  # put cordinates of each cell in single df
  coordinates_mfb <- rbind(neg$coords, pos$coords)
  colnames(coordinates_mfb) <- c('x', 'y')
  coordinates_mfb <- as.data.frame(coordinates_mfb)
  
  # For each point, n, set to xn,yn - xcentroid_ves_psp,ycentroid_ves_psp = (xn-xC, yn-yC); new x,y
  coordinates_mfb$x_centre <- coordinates_mfb$x - centroid_ves_psp$x
  coordinates_mfb$y_centre <- coordinates_mfb$y - centroid_ves_psp$y
  
  # use polar_theta to generate theta with respect to centroid_ves_psp
  coordinates_mfb$theta <- polar_theta(coordinates_mfb$x_centre, coordinates_mfb$y_centre)
  
  # find the peak density of theta
  density_theta <- density(coordinates_mfb$theta)
  Peak_theta <- density_theta$x[which.max(density_theta$y)]
  
  # use this peak theta to adjust all thetas, setting peak to '90'
  Theta_adjust <- (360-Peak_theta) + 90
  coordinates_mfb$theta_adjusted <- coordinates_mfb$theta + Theta_adjust
  
  # for those now over 360, convert back to 0-360 range
  coordinates_mfb$thetan_adjusted_corrected <- coordinates_mfb$theta_adjusted
  
  for (i in 1:nrow(coordinates_mfb)) {
    if (coordinates_mfb$theta_adjusted[i] > 360) {
      coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i] - 360
    }
    else {
      coordinates_mfb$thetan_adjusted_corrected[i] <- coordinates_mfb$theta_adjusted[i]
    }
  }
  
  # add image name as factor to data frame
  coordinates_mfb$image_number <- as.factor(image_number)
  
  return(coordinates_mfb)
  
}

# package all commands in one function
animal_images_polar_new <- function(){
  # select folder containing rois
  initial.dir<-getwd() # current working directory
  data_dir <- jchoose.dir(caption = "Select individual animal directory containing p,n,v ROIs")
  
  # get list of separate image folder
  image_dirs <- list.dirs(path = data_dir)
  image_dirs <- image_dirs[-1] # remove the parent directory
  
  # create data frame
  df <- lapply(image_dirs, roi_to_polar_new)
  
  return(df)
  
  setwd(initial.dir)
}

# run on an animal specific parent folder
WT3_1_allMFBs_new <- animal_images_polar_new()
WT3_2_allMFBs_new <- animal_images_polar_new()
WT3_3_allMFBs_new <- animal_images_polar_new()
WT3_4_allMFBs_new <- animal_images_polar_new()
WT3_5_allMFBs_new <- animal_images_polar_new()
WT3_6_allMFBs_new <- animal_images_polar_new()

# to plot on a single plot

big_WT3_1_new <- do.call(rbind,WT3_1_allMFBs_new) # combine into a single df
big_WT3_2_new <- do.call(rbind,WT3_2_allMFBs_new)
big_WT3_3_new <- do.call(rbind,WT3_3_allMFBs_new)
big_WT3_4_new <- do.call(rbind,WT3_4_allMFBs_new)
big_WT3_5_new <- do.call(rbind,WT3_5_allMFBs_new)
big_WT3_6_new <- do.call(rbind,WT3_6_allMFBs_new)

# add identifier for animal as we will combine into a single df
big_WT3_1_new$animal <- as.factor('1')
big_WT3_2_new$animal <- as.factor('2')
big_WT3_3_new$animal <- as.factor('3')
big_WT3_4_new$animal <- as.factor('4')
big_WT3_5_new$animal <- as.factor('5')
big_WT3_6_new$animal <- as.factor('6')

# combine all into a single df
All_CCl4_animals_new <- rbind(big_WT3_1_new, big_WT3_2_new, big_WT3_3_new, big_WT3_4_new, big_WT3_5_new, big_WT3_6_new)

# for the summary of each animal together
svglite("Density polar_all_new.svg")

ggplot(All_CCl4_animals_new, aes(x=thetan_adjusted_corrected, fill=animal)) +
  geom_density(alpha=.7) +
  scale_fill_brewer(palette = 'Blues', direction=-1) + 
  xlab("Aligned MFB angle (degrees)") +
  ylab("MFB population density") +
  theme_minimal() +
  theme (legend.key = element_blank(),
         text=element_text(family="Calibri"),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))

dev.off()

svglite("Polar histogram_all.svg")

ggplot(All_CCl4_animals_new, aes(x = thetan_adjusted_corrected, fill=animal)) + 
  geom_histogram(breaks = seq(0, 360, 10), alpha=0.7) + coord_polar(start = 0) + theme_minimal() + 
  scale_fill_brewer(palette = 'Blues', direction=-1) + 
  ylab("Count/angle (degree)") + 
  scale_x_continuous("", limits = c(0, 360), breaks = seq(0, 360, 30), labels = seq(0, 360, 30)) +
  theme (legend.key = element_blank(),
         text=element_text(family="Calibri"),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12))

dev.off()

# to access the data from all central veins of all animals
all_hist_new <- ggplot(All_CCl4_animals_new, aes(x=thetan_adjusted_corrected)) + 
  geom_histogram(breaks = seq(0, 360, 5), alpha=.5) + theme_minimal() + scale_fill_brewer()
count_all_new <- print(all_hist_new)
raw_count_all_CCl4_new <- count_all_new[["data"]][[1]]

# to fit sine
# in ggplot

reslm_all_new <- lm(raw_count_all_CCl4_new$density ~ sin(2*pi/per*raw_count_all_CCl4_new$x)+cos(2*pi/per*raw_count_all_CCl4_new$x))
summary(reslm_all_new)

svglite("All_density_fitted_new.svg")

ggplot(raw_count_all_CCl4_new, aes(x=x, y=density)) + geom_point() + ylim(0, max(raw_count_all_CCl4_new$density)) + theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  xlab("Aligned MFB angle (degrees)") +
  ylab("Total MFB population density") +
  geom_line(aes(x=x, y=fitted(reslm_all_new))) 

dev.off()

# damped sine
reslm_damped <- lm(raw_count_all_CCl4_new$density ~ exp(-0.01*(2*pi/per*raw_count_all_CCl4_new$x))*sin(2*pi/per*raw_count_all_CCl4_new$x)+cos(2*pi/per*raw_count_all_CCl4_new$x))
summary(reslm_damped)

ggplot(raw_count_all_CCl4_new, aes(x=x, y=density)) + geom_point() + ylim(0, max(raw_count_all_CCl4_new$density)) + theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  xlab("Aligned MFB angle (degrees)") +
  ylab("Total MFB population density") +
  geom_line(aes(x=x, y=fitted(reslm_damped))) 
