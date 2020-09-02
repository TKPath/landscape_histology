# to use the xy saved centroids from thyroid follicle outputs for spatstat
# Tim Kendall, University of Edinburgh
# Tim.Kendall@ed.ac.uk

# August 2020


# for each directory of centroid csv files, make a list linking filename to diagnosis
require(tidyverse)
require(RColorBrewer)
require(ggplot2)
require(spatstat)
require(svglite)
require(extrafont)
require(viridis)

data_dir <- rstudioapi::selectDirectory() # select directory, repeat for each 3

raw <- list.files(path=data_dir, all.files = FALSE, full.names = TRUE, pattern = "\\.csv$",
                  recursive = TRUE, include.dirs = FALSE, no.. = FALSE)


key_tib <- tibble(diagnosis = 'normal', filename = raw)
key_temp1 <- tibble(diagnosis = 'hashimotos', filename = raw)

key_tib <- bind_rows(key_tib, key_temp1)
key_tib$ID <- gsub("(.*crops/)(.*)(\\.svs)(.*$)", "\\2", key_tib$filename) # get ID, will need to be changed depending on filenames

# function to run through list of files extracting points and converting to ppp

xy_to_ppp <- function(x){
  require(spatstat)
  window_thyroid <- owin(c(0,1773), c(0,1850)) # set common window for thyroid crops, alter as needed
  raw_csv <- read.csv(x, header = TRUE, row.names = 1) # read csv
  ID <- gsub("(.*crops/)(.*)(\\.svs)(.*$)", "\\2", x) # get ID, alter as needed
  raw_ppp <- ppp(x = raw_csv$X, raw_csv$Y, window = window_thyroid) 
  marks(raw_ppp) <- ID
  return(raw_ppp)
}

# test normal versus Hashimotos based on follicle patterning

hashi_ppp <- lapply(key_tib$filename, xy_to_ppp)
hashi_ppp_rescale <- lapply(hashi_ppp, rescale, s = 505.83, unitname = "mm") # convert from pixels to mm, 505.83 px/mm here, alter as needed
names(hashi_ppp_rescale) <- key_hashi$ID

hashi_solist <- as.solist(hashi_ppp_rescale) # create solist

# construct hyperframe from ppp set

hashi_hyper <- hyperframe() # initialise hyperframe
hashi_hyper$cases <- hashi_solist # add thyroid ppps
hashi_hyper$ID <- key_hashi$ID # add corresponding IDs
hashi_hyper$diagnosis <- key_hashi$diagnosis # add corresponding diagnosis

# plot the frames for exploratory check
plot(hashi_hyper, arrange=FALSE) 

# other simple descriptors for the ppps

# simple intensity
intensity <- sapply(hashi_ppp_rescale, intensity)

# Clark-Evans test of CSR
ce_test <- sapply(hashi_ppp, clarkevans, correction = 'Donnelly')

# Hopskel-Skellam index
hop_test <- sapply(hashi_ppp, hopskel)

# average nearest neighbour distance
nn_mean <- sapply(hashi_ppp, function(x) mean(nndist(x)))

# combine these simple metrics
csr_indices <- tibble(ID = key_hashi$ID, diagnosis = key_hashi$diagnosis, intensity = intensity, clarkEvans = ce_test,
                      hopskel = hop_test, nearestNeighbour = nn_mean)

###################### plots ##########################
# brewer dark2 pallete for thyroid
require(RColorBrewer)
require(ggplot2)
require(extrafont)

# plot example ppp

example_ppp_df <- tibble(x = hashi_ppp[['GTEX-1J1R8-0126']]$x, y = hashi_ppp[['GTEX-1J1R8-0126']]$y)

svglite('output/ppp_thyroid_normal.svg', width = 6, height = 6)
ggplot(aes(x,y), data = example_ppp_df) + geom_point(colour = pal_thyroid[[1]], size = 2) + coord_fixed(ratio = 1) +
  theme_bw() + 
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.title=element_blank(), axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())
dev.off()

# plot of simple metric based on group
csr_long <- csr_indices %>% gather(metric, value, intensity:nearestNeighbour, factor_key = TRUE) # prepare for ggplot
csr_long$diagnosis <- as.factor(csr_long$diagnosis)

svglite::svglite('output/thyroid_clarkHopskel.svg', width = 5, height = 3)
ggplot(filter(csr_long, metric %in% c('clarkEvans', 'hopskel')), aes(x = metric, y = value, fill = diagnosis)) +
  geom_boxplot(alpha = 0.8, outlier.color = NA, show.legend = FALSE) +
  stat_boxplot(geom = 'errorbar') +
  geom_point(position = position_jitterdodge(), show.legend = FALSE) + theme_minimal() +
  scale_fill_brewer(palette = 'Dark2') +
  theme(text = element_text(family = "Calibri"), axis.title = element_text(size = 15)) +
  ylim(0, 2.25) +
  labs(x = 'Metric', y = 'Metric value') +
  scale_x_discrete(labels=c("clarkEvans" = "Clark & Evans Aggregation Index", "hopskel" = "Hopkins-Skellam Index"))
dev.off()

svglite::svglite('output/thyroid_nn.svg', width = 5, height = 3)
ggplot(filter(csr_long, metric == 'nearestNeighbour'), aes(x = metric, y = value, fill = diagnosis)) +
  geom_boxplot(alpha = 0.8, outlier.color = NA) +
  stat_boxplot(geom = 'errorbar') +
  scale_fill_brewer(palette = 'Dark2', name = 'Diagnosis', labels = c('Normal', 'Hashimoto')) +
  geom_point(position = position_jitterdodge(), show.legend = FALSE) + theme_minimal() +
  ylim(0, 0.125) +
  theme(text = element_text(family = "Calibri"), axis.title = element_text(size = 15)) +
  labs(x = NULL, y = 'Nearest neighbour distance (mm)') 
  #scale_x_discrete(labels=c('nearestNeighbour' = 'Mean Nearest Neighbour distance'))
dev.off()

svglite::svglite('output/thyroid_intensity.svg', width = 5, height = 3)
ggplot(filter(csr_long, metric == 'intensity'), aes(x = metric, y = value, fill = diagnosis)) +
  geom_boxplot(alpha = 0.8, outlier.color = NA) +
  stat_boxplot(geom = 'errorbar') +
  scale_fill_brewer(palette = 'Dark2', name = 'Diagnosis', labels = c('Normal', 'Hashimoto')) +
  geom_point(position = position_jitterdodge(), show.legend = FALSE) + theme_minimal() +
  theme(text = element_text(family = "Calibri"), axis.title = element_text(size = 15),
        axis.text.x = element_blank()) +
  labs(x = NULL, y = 'Intensity (point/mm)') +
  #scale_x_discrete(labels=c('intensity' = 'Intensity')) +
  ylim(0,100)
dev.off()

# calculate and plot separate case functions

key_hashi_only <- key_tib %>% filter(diagnosis == 'hashimotos') # subset
hashi_only_ppp <- lapply(key_hashi_only$filename, xy_to_ppp)
hashi_only_ppp_scale <- solapply(hashi_only_ppp, rescale, s = 505.83, unitname = "mm") # scale

Lest_hashi_only <- lapply(hashi_only_ppp_scale, Lest)
Fest_hashi_only <- lapply(hashi_only_ppp_scale, Fest)
Gest_hashi_only <- lapply(hashi_only_ppp_scale, Gest)
Jest_hashi_only <- lapply(hashi_only_ppp_scale, Jest)

key_normal_only <- key_tib %>% filter(diagnosis == 'normal') # subset
normal_only_ppp <- lapply(key_normal_only$filename, xy_to_ppp)
normal_only_ppp_scale <- solapply(normal_only_ppp, rescale, s = 505.83, unitname = "mm") # scale

Lest_normal_only <- lapply(normal_only_ppp_scale, Lest)
Fest_normal_only <- lapply(normal_only_ppp_scale, Fest)
Gest_normal_only <- lapply(normal_only_ppp_scale, Gest)
Jest_normal_only <- lapply(normal_only_ppp_scale, Jest)

svglite::svglite('output/thyroid_Fest_individuals.svg', width = 5, height = 3) # change for G, L, J
ggplot() + # normal is cyan, pal_thyroid[1]
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[1]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[2]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[3]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[4]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[5]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[6]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[7]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[8]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[9]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_normal_only[[10]], inherit.aes = TRUE, colour = pal_thyroid[1]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[1]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[2]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[3]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[4]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[5]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[6]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[7]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[8]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[9]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  geom_line(aes(x=r, y=km), data = Fest_hashi_only[[10]], inherit.aes = TRUE, colour = pal_thyroid[2]) +
  theme_minimal() +
  xlab(expression(paste("Radius (mm)"))) +
  ylab("F(r)") +
  xlim(0,0.3) +
  theme(text = element_text(family = "Calibri"), axis.title = element_text(size = 15)) 
dev.off()

############# statistical comparisons ############################


t.test(value ~ diagnosis, filter(csr_long, metric == 'intensity'))

t.test(value ~ diagnosis, filter(csr_long, metric == 'nearestNeighbour'))

t.test(value ~ diagnosis, filter(csr_long, metric == 'clarkEvans')) 

t.test(value ~ diagnosis, filter(csr_long, metric == 'hopskel')) 

test_Lest_hashi <- studpermu.test(hashi_hyper, 
                                  cases ~ diagnosis, summaryfunction = Lest, use.Tbar = FALSE)

test_Fest_hashi <- studpermu.test(hashi_hyper, 
                                  cases ~ diagnosis, summaryfunction = Fest, use.Tbar = FALSE) 

test_Gest_hashi <- studpermu.test(hashi_hyper, 
                                  cases ~ diagnosis, summaryfunction = Gest, use.Tbar = FALSE) 

test_Jest_hashi <- studpermu.test(hashi_hyper, 
                                  cases ~ diagnosis, summaryfunction = Jest, use.Tbar = FALSE) 




