# to generate landscape metrics from classified raster images

# Tim Kendall, University of Edinburgh
# Tim.Kendall@ed.ac.uk

# August 2020

single_metric_choice <- function(x, metric_choice){
  require(raster)
  require(rgdal)
  library(purrr)
  library(landscapemetrics)
  
  geo_input <- raster(x) # load image
  crs(geo_input) <- '+proj=utm +zone=30 +datum=WGS84' # apply projection, utm, zone for london is #30, arbitrary choice
  
  metrics_agg <- calculate_lsm(geo_input,
                               type = metric_choice,
                               level = c("class", "landscape")) # selected metrics of class and landscape level
  
  return(metrics_agg) # return the object
}

landscape_liver_choice <-  function(){ # overall function starting on folder of categorical raster images and single_metric function defining the categories
  require(rstudioapi)
  require(svDialogs)
  require(pbapply)
  
  input_dir <- selectDirectory(caption = 'Select input directory',
                               label = 'Select') # get input folder that contains classified images
  
  metric_choice <- select.list(choices = c('complexity metric', 'diversity metric', 'area and edge metric',
                                           'aggregation metric', 'shape metric', 'core area metric'),
                               multiple = TRUE, title = 'Select metrics to calculate', graphics = TRUE) # user selects a metric
  
  raw_images <- list.files(path=input_dir, all.files = FALSE, full.names = TRUE, pattern = "\\.tif$",
                           recursive = TRUE, include.dirs = FALSE, no.. = FALSE) # get list of images, full paths
  
  short_name <- gsub("(.*input/)(.*)(\\.svs)(.*$)", "\\2", raw_images, ignore.case = TRUE)  # get image name without extension; will need changing depend on how files are names
  
  names(raw_images) <- short_name # apply the preserved names attribute to carry forward to output list
  
  metric_list <- pblapply(raw_images, single_metric_choice, metric_choice = metric_choice) # apply single_metric_choice function to each image in list
  
  for (i in 1:length(metric_list)){
    metric_list[[i]]$image_name <- names(metric_list[i]) # add unique image id to each metric to allow combination
  }
  
  require(dplyr)
  metric_list_combined <- bind_rows(metric_list) # combine into a single tibble
  
  metric_list_combined$new_metric <- paste(metric_list_combined$level, metric_list_combined$class,
                                           metric_list_combined$metric, sep = "_") # level, class and metric as a combination need to be joined as a new column
  
  write.csv(metric_list_combined, paste(metric_choice, "csv", sep = ".")) # write data to csv
  
  require(reshape2)
  metric_wide <- dcast(metric_list_combined, image_name ~ new_metric, value.var = "value") # convert to wide form
  
  metric_wide$source_id <- gsub("(Geo_classified_)(.*)(_tiles_)(.*)", "\\2",
                                metric_wide$image_name, ignore.case = TRUE)  # add root image id
  
  write.csv(metric_wide, paste(metric_choice, "_wide", ".csv", sep = "")) # write data to csv
  
  return(metric_wide) # return the full wide list
}

###########################

# run for each metric type of interest e.g.

complexity_wide <- landscape_liver_choice() # select complexity metric
aggregation_wide <- landscape_liver_choice() # select aggregation metric - slow, run overnight. elapsed=13h 43m 38s.
diversity_wide <- landscape_liver_choice() # select diversity metric
area_wide <- landscape_liver_choice() # area and edge metric
core_wide <- landscape_liver_choice() # core area metric
shape_wide <- landscape_liver_choice() # shape metric