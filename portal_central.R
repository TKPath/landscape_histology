# portal-central relationship and lobule area
# starts with a roi from FIJI, annotated portal tracts (pt) or central veins (cv) from each image, separate folder for each pair

# Tim Kendall, University of Edinburgh
# Tim.Kendall@ed.ac.uk

# August 2020

require(spatstat)
require(ggplot2)
require(svglite)
require(RImageJROI)
require(extrafont)
require(rJava)
require(rChoiceDialogs)

# all shared the same window
window_pt <- owin(c(0,1464), c(0,1956))

# function to get nearest p-c distances, input is a folder returning a list of distances for each image

human_portal_central <- function(x){
  setwd(x)
  
  # get name of each image based on image directory
  image_number <- gsub("(.*/)(.*$)", "\\2", x)
  
  pt <- read.ijroi('pt.roi')
  cv <- read.ijroi('cv.roi')
  
  # combine cells into a single ppp and return
  pt_spat <- ij2spatstat(pt, window=window_pt)
  cv_spat <- ij2spatstat(cv, window=window_pt)
  
  marks(pt_spat) <- 'pt'
  marks(cv_spat) <- 'cv'
  
  plot(pt_spat, main = image_number)
  plot(cv_spat, main = image_number)
  
  # For each cv, calculate distance to nearest 6 pts and return
  p_c_distances <- nndist(cv_spat, pt_spat, k=1:6)
  p_c_distances$mean <- 
  return(p_c_distances)
}

# wrapper funtion to use human_portal_central
Portal_Central_wrapper_function <- function(){
  window_pt <- owin(c(0,1464), c(0,1956))
  
  # select folder containing rois
  initial.dir <- getwd() # current working directory
  data_dir <- jchoose.dir(caption = "Select parent folder containing case folders")
  
  # get list of separate image folder
  case_dirs <- list.dirs(path = data_dir)
  case_dirs <- case_dirs[-1] # remove the parent directory
  
  # create data frame
  df <- lapply(case_dirs, human_portal_central)
  
  return(df)
  
  setwd(initial.dir)
}

Portal_Central_human <- Portal_Central_wrapper_function()
# check plots

PC_means <- lapply(Portal_Central_human, rowMeans) # each value is the mean p-c distance of each cv to 6 nearest pts

case_dirs <- list.dirs(path = data_dir, full.names = FALSE)
case_dirs <- case_dirs[-1] # remove the parent directory
names(PC_means) <- case_dirs

PC_means_df <- data.frame(PC = unlist(PC_means), ID = rep(names(PC_means), sapply(PC_means,length)))
# set diagnostic types
PC_means_df$diagnosis <- 'Normal'
PC_means_df$diagnosis[PC_means_df$ID %in% c('a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10')] <- 'Obstructed' 
PC_means_df$diagnosis <- as.factor(PC_means_df$diagnosis)

# need to convert from pixels to actual distance to allow inter-image comparisons in the end
# for thumbnails created from ndpi via split, .275 pixels/micron
PC_means_df$PCmicron <- PC_means_df$PC / 0.275

ggplot(PC_means_df, aes(x=ID, y=PCmicron)) + geom_boxplot() + theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))

# for each lobule (cv) calculate area using mean (of nearest 6 pts)
PC_means_df$lobule_area <- ((3*sqrt(3))/2)*(PC_means_df$PCmicron^2) # in square microns

# area in square mm
ggplot(PC_means_df, aes(x=ID, y=lobule_area/1000000)) +
  stat_boxplot(geom ='errorbar') + geom_boxplot() +
  ylab(expression(paste("Lobule area (",mm^2, ")"))) + scale_x_discrete(labels = 1:10) + xlab("Case") +
  theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))

ggplot(PC_means_df, aes(x=ID, y=lobule_area/1000000, fill = diagnosis)) +
  stat_boxplot(geom ='errorbar') + geom_boxplot() +
  ylab(expression(paste("Lobule area (",mm^2, ")"))) + scale_x_discrete(labels = 1:10) + xlab("Case") +
  theme_minimal() +
  theme (legend.key = element_blank(), text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))

# get mean lobule area for each case
require(plyr)

obstructed_set <- ddply(PC_means_df, c("ID", "diagnosis"), summarise,
                        mean = mean(lobule_area), sd = sd(lobule_area),
                        sem = sd(lobule_area)/sqrt(length(lobule_area)))

svglite('Mean lobule size normal obstructed.svg')

ggplot(obstructed_set, aes(x=diagnosis, y=mean/1000000, fill = diagnosis)) +
  stat_boxplot(geom ='errorbar') + geom_boxplot() + geom_jitter(width = 0.25) +
  ylab(expression(paste("Lobule area (",mm^2, ")"))) + xlab("Diagnosis") +
  theme_minimal() + scale_fill_manual('Diagnosis',values=c("lightsalmon4","olivedrab4")) +
  theme (legend.key = element_blank(), legend.title = element_blank(),
         text=element_text(family="Calibri"), axis.text = element_text(size = 12), axis.title = element_text(size = 12))

dev.off()

ddply(obstructed_set, "diagnosis", summarise, Mean = mean(mean), sd = sd(mean), sem = sd(mean)/sqrt(length(mean)))

shapiro.test(subset(obstructed_set$mean, obstructed_set$diagnosis == "Normal"))
shapiro.test(subset(obstructed_set$mean, obstructed_set$diagnosis == "Obstructed"))

t.test(mean ~ diagnosis, data = obstructed_set)

####################################################

# to plot illustrative portal central

window_pt <- owin(c(0,1464), c(0,1956))

pt <- read.ijroi('Portal central set/a1/pt.roi')
cv <- read.ijroi('Portal central set/a1/cv.roi')

# combine cells into a single ppp and return
pt_spat <- ij2spatstat(pt, window=window_pt)
cv_spat <- ij2spatstat(cv, window=window_pt)

marks(pt_spat) <- 'pt'
marks(cv_spat) <- 'cv'

all_spat <- superimpose(pt_spat, cv_spat)
all_spat_df <- as.data.frame(all_spat$x)
all_spat_df <- cbind(all_spat_df, all_spat$y)
all_spat_df <- cbind(all_spat_df, all_spat$marks)
colnames(all_spat_df) <- c('x', 'y', 'marks')

svglite('a1.svg')

ggplot(aes(x,y,group=marks), data=all_spat_df) +
  geom_point(aes(colour=marks)) + coord_fixed(ratio = 1) +
  scale_colour_manual(values=c('blue', 'orange')) +
  theme_bw() + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.title=element_blank(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

dev.off()
