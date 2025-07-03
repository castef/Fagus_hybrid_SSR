### Allenwiller and Waldi phenology ###

library(plyr)
library(xlsx)
library(dplyr)
#library(gplots)
library(ggpubr)
library(readr)
library(ggplot2)
library(lubridate) ## working with data formats
library(tidyverse)
#library("graph4lg")
library("sf")  # for rasters
library("sp") # for rasters
library("raster") # for rasters
library("tidyverse") 
library("terra") # for rasters
library(stringr) # to add the zeros in Waldi dataset
library(jsonlite)
library(data.table)


rm(list=ls())
#path = "C:/Users/Camilla/Dropbox (Old (1))/Dropbox/WSL_PhD/Projects/Hybridization/Phenology/"
setwd(dir="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean/Phenology/")
path = "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Data_clean/Phenology/"

# custom colors for subspecies
mycols <- c("orientalis" = "#482173FF", "sylvatica" = "#ffcc00", "hybrid" = "#25858EFF")
colSet <- scale_colour_manual(values = mycols)
fillSet <- scale_fill_manual(values = mycols)

######## dataset preparation and cleaning #############################
all_data = read.csv(paste0(path,"Spring_phenology_all_sites_raw.csv"))

### correcting the old ID in Waldi
add_zeros = function(id) {
  # extract the letters and digits
  matches <- str_match(id, "^([A-Z]+)([0-9]+)$")
  # return original id if it does not match the pattern
  if (is.na(matches[1,1])) {
    return(id)  
  }
  letters <- matches[1,2]
  digits <- matches[1,3]
  # pad digits with leading zeros until length is 3
  padded_digits <- str_pad(digits, width = 3, side = "left", pad = "0")
  # return the combined string
  return(paste0(letters, padded_digits))
}

# Subset the data and apply the function
all_data[all_data$site == "Waldi",] =  all_data[all_data$site == "Waldi",] %>%
  mutate(tree.id = sapply(tree.id, add_zeros))
# remove offspring
all_data = subset(all_data, generation == "adult")
# convert the variables into numeric and factor for later analysis
all_data[,6:9] = sapply(all_data[,6:9],as.numeric)
all_data[, 6:9][is.na(all_data[6:9])] = 0  ## converting all NA into 0
all_data$species = factor(as.character(all_data$species))
# calculating DOY and extract the year 
all_data$doy = as.numeric(as.character(as.POSIXlt(all_data$date, format = "%d.%m.%Y")$yday))
all_data$year = year(as.POSIXlt(all_data$date, format = "%d.%m.%Y"))
## calculate average stage for the data where is not calculated and leave the current data where is calculated
all_data$average.stage = as.numeric(as.character(all_data$average.stage))
all_data = all_data %>%
  mutate(average.stage = ifelse(
    is.na(average.stage),
    (1 * stage1 + 2 * stage2 + 3 * stage3 + 4 * stage4) / 100,
    average.stage
  ))

# sample size for each category
all_data_sum = all_data %>% group_by(site, year, species) %>% summarise(N = n())

ggplot(all_data_sum, aes(y = N, x = species, fill = species)) + geom_bar (stat = "identity") + theme_bw()+
  facet_grid(site~year ) + fillSet + geom_text(aes(label=N), vjust=0.1 , size = 4) 

# visual overview of the data - time in vertical to better see the shift in the season
## orientalis point in 2022 are not very visible

ggplot(all_data, aes(as.numeric(as.character(doy)), average.stage, fill = species))+
  geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 2) + 
  fillSet+ colSet+ theme_bw()+
  stat_smooth(aes(group=species, fill = species, col = species))+ 
  ylab("Average Phenological Stage")+
  xlab("DOY (Day of Year)")+
  facet_grid(site~year )


ggplot(data = all_data, aes(x=as.factor(doy), y=average.stage, fill=species))+
  geom_boxplot(position = position_dodge())+ 
  fillSet+ colSet+ theme_bw()+
  scale_y_continuous( breaks = seq(0, 200, 50))+
  ylab("Average Phenological Stage")+
  xlab("DOY (Day of Year)")+
  facet_grid(site~year )


##### cleaning dataset
# average stage < 1?
hist(all_data$average.stage)
subset(all_data, average.stage<1)
# mistake: transform all average.stage <1 to 1
all_data[which(all_data$average.stage < 1), "average.stage"] = 1

# check if there are trees where the averge stage is not increasing with doy (due to errors)
errors <- all_data %>%
  group_by(site, species, tree.id, year) %>%
  arrange(doy) %>%
  mutate(
    next_stage = lead(average.stage),
    next_doy = lead(doy),
    is_error = ifelse(!is.na(next_stage) & next_stage < average.stage, TRUE, FALSE)
  ) %>%
  filter(is_error == TRUE) %>%
  dplyr::select(site, species, tree.id, year, doy, next_doy, average.stage, next_stage)
errors


## check visually the errors
error_trees <- errors %>%
  dplyr::select(site, species, tree.id, year) %>%
  distinct()

error_data <- all_data %>%
  inner_join(error_trees, by = c("site", "species", "tree.id", "year"))

ggplot(error_data, aes(as.numeric(as.character(doy)), average.stage, col = tree.id, label = tree.id)) +
  geom_line() +geom_text(position=position_jitter(width=0,height=0), size = 3)+
  facet_grid(year~species) +
  theme(legend.position = "none")

# removing wrong observations 
error_tabi <- paste(error_trees$site, error_trees$species, error_trees$tree.id, error_trees$year)
cleaned_data <- all_data[!paste(all_data$site, all_data$species, all_data$tree.id, all_data$year) %in% error_tabi, ]

dim(all_data)
dim(cleaned_data)

ggplot(cleaned_data, aes(as.numeric(as.character(doy)), average.stage, col = tree.id,label = tree.id))+
  geom_line()+  
  #geom_text(position=position_jitter(width=0,height=0), size = 3)+
  facet_grid(year~species)+
  theme(legend.position = "none")

# plot single trees to check everything makes sense or if there are outliers
ggplot(data = cleaned_data, aes(x = as.factor(doy), y = average.stage, fill = as.factor(species))) +
  geom_boxplot() +
  fillSet +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  theme_bw() +
  ylab("Average Phenological Stage") +
  xlab("DOY") +
  facet_grid( year~ site)

# checking closely outlier
subset(cleaned_data, year== 2023 & site == "Waldi" & doy >= 110 & average.stage <= 2)

# removing
cleaned_data <- cleaned_data[!(cleaned_data$year == 2023 & cleaned_data$site == "Waldi" & cleaned_data$doy >= 110 & cleaned_data$average.stage <= 2), ] 
fwrite(cleaned_data, "Spring_phenology_all_sites_cleaned.csv")

######## getting T data from ERA5 and calculating GDD and CD #########
#loading data from Waldi
json_file <-"C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Phenology/T_data/era5_waldi_2.json"
temp_data_station <- fromJSON(json_file, flatten=TRUE)
temp_matrix=as.data.frame(temp_data_station$hourly$temperature_2m)
names(temp_matrix)="value"
temp_matrix$time=temp_data_station$hourly$time
temp_matrix$date=as.Date(temp_data_station$hourly$time)
res_final=temp_matrix %>% group_by(date) %>% summarise(temp_min=min(value),temp_max=max(value),temp_mean=mean(value))
res_final$lat=temp_data_station$latitude
res_final$long=temp_data_station$longitude
res_final$site= "Waldi"
res_final_w = res_final

#loading data from Allenwiller
json_file <-"C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Phenology/T_data/era5_allenwiller_2.json"
temp_data_station <- fromJSON(json_file, flatten=TRUE)
temp_matrix=as.data.frame(temp_data_station$hourly$temperature_2m)
names(temp_matrix)="value"
temp_matrix$time=temp_data_station$hourly$time
temp_matrix$date=as.Date(temp_data_station$hourly$time)
res_final=temp_matrix %>% group_by(date) %>% summarise(temp_min=min(value),temp_max=max(value),temp_mean=mean(value))
res_final$lat=temp_data_station$latitude
res_final$long=temp_data_station$longitude
res_final$site= "Allenwiller"
res_final_a = res_final

#loading data from Leiselheim
json_file <- "C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Phenology/T_data/era5_leiselheim_2.json"
temp_data_station <- fromJSON(json_file, flatten=TRUE)
temp_matrix=as.data.frame(temp_data_station$hourly$temperature_2m)
names(temp_matrix)="value"
temp_matrix$time=temp_data_station$hourly$time
temp_matrix$date=as.Date(temp_data_station$hourly$time)
res_final=temp_matrix %>% group_by(date) %>% summarise(temp_min=min(value),temp_max=max(value),temp_mean=mean(value))
res_final$lat=temp_data_station$latitude
res_final$long=temp_data_station$longitude
res_final$site= "Leiselheim"
res_final_l = res_final

# merge all datasets
res_final = rbind(res_final_a,res_final_l,res_final_w)
#fwrite(res_final,"temperatures_all_sites.csv")

# load T data
t_dat = fread("temperatures_all_sites.csv", header = T, sep = ",")
# change data format and calculate DOY
t_dat$date = format(ymd(t_dat$date), "%d.%m.%Y")
t_dat$date2 = as.Date(t_dat$date, "%d.%m.%Y") ## this format is required for the for loop after
t_dat$year = year(as.POSIXlt(t_dat$date, format = "%d.%m.%Y"))
t_dat$doy = as.POSIXlt(t_dat$date, format = "%d.%m.%Y")$yday

# GDD for each day starting at the beginning of the year (if GDD is < 0, GDD = 0 (< 0 growth does not progress))
t_dat = t_dat %>% group_by(year, site) %>% mutate(GDD = pmax(0, temp_mean - 5)) # T base = 5 Degrees

# Cumulated GDD
t_dat = t_dat %>% group_by(year, site) %>% mutate(Cum_GDD = cumsum(GDD)) 

# Cumulative Chilling Days
# defining growing season (from November of previous year, till November current year)
t_dat = t_dat %>% group_by(year, site) %>%
  mutate(season = case_when(month(date2) >= 11 ~ year(date2) + 1, 
                            month(date2) < 11 ~ year(date2), 
                            date2 >= 2023-11-01 ~ 2024)) %>%
  ungroup()%>%
  # start counting Chilling days from 1rst novermber of the previous year
  group_by(site,season) %>%
  # count as 1 if the day was chilling day (< 5°C)
  mutate(below5 = ifelse(temp_mean < 5, 1, 0),
         # counting number of chilling days 
         CD = cumsum(below5))  


#fwrite(t_dat,"T_GDD_CD_all_sites.csv")

# visual check of annual trends per site
t_dat = fread("T_GDD_CD_all_sites.csv")

# checking T trend across years
ggplot(t_dat, aes(x=doy )) +
  geom_line( aes(y=temp_mean,  group = as.factor(year),color= year), linewidth = 0.7) + 
  scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~site, ncol = 1)

# GDD and CD across sites and years
ggplot(t_dat, aes(x=doy )) +
  geom_area( aes(y=CD*10, group = as.factor(season), fill = season), position = "identity",linewidth = 0.9,stat = "align", alpha = 0.3) + # multiply by 10 to get the same range than GDD
  scale_y_continuous(name = "GDD", sec.axis = sec_axis(~./10, name="Accumulated Chilling Days"))+
  scale_fill_gradient(low = "yellow", high = "red")+ theme_bw()+
  geom_line( aes(y=Cum_GDD,  group = as.factor(year),color= year), linewidth = 0.9) + 
  scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~site, ncol = 1)


GDD_CD_sites = ggplot(subset(t_dat,site =="Waldi"), aes(x=date2 )) +
  geom_area( aes(y=CD*10),position = "identity",linewidth = 0.9, alpha = 0.3) + # multiply by 10 to get the same range than GDD
  scale_y_continuous(name = "Growing Degree Days", sec.axis = sec_axis(~./10, name="Accumulated Chilling Days"))+
  geom_line( aes(y=Cum_GDD), linewidth = 0.9) + theme_bw()+
  labs(x = "Year")+
  scale_color_gradient(low = "yellow", high = "red")+
  facet_wrap(~site, ncol = 1)
GDD_CD_sites

png(filename ="C:/Users/stefanin/Dropbox/WSL_PhD/Projects/Hybridization/Phenology/GDD_CD_sites.png", width = 1500, height = 1000, res =200)
plot(GDD_CD_sites)
dev.off()

## add GDD and CD to the data  
# loading phenological data
p_dat = fread("Spring_phenology_all_sites_cleaned.csv")
# loading climate data
t_dat = read.csv("T_GDD_CD_all_sites.csv")
all_dat = merge(p_dat, t_dat[, c("site", "year", "doy", "Cum_GDD","CD")],  by.x = c("year", "doy", "site"),  by.y = c("year", "doy", "site"), all.x = T)
fwrite(all_dat, "Spring_phenology_all_sites_final.csv")

####  comparing average stage ####
all_dat <- fread("Spring_phenology_all_sites_final.csv")

# check sample size for each date
sample_sizes <- all_dat %>%
  group_by(site, year, date, species) %>%
  summarise(sample_size = n(), .groups = "drop")
# avergae per year
avg_samples <- sample_sizes %>%
  group_by(site, year, species) %>%
  summarise(avg_sample = mean(sample_size), .groups = "drop")
avg_labels <- avg_samples %>%
  pivot_wider(names_from = species, values_from = avg_sample, names_prefix = "N") %>%
  mutate(label = paste0(
    "Ns = ", round(Nsylvatica), "\n",
    "No = ", round(Norientalis)
  ))

# over years DOY (waldi)
  ggplot(data = subset(all_dat, site == "Waldi"), 
         aes(x = as.factor(doy), y = average.stage, fill = species)) +
  geom_point(position = position_jitter(h = 0.1, w = 0.1), shape = 21, alpha = 0.5, size = 2) + 
  #geom_text(data = subset(avg_labels, site == "Waldi"), aes(x = Inf, y = -Inf, label = label),inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5)+
  fillSet + colSet + theme_bw() +
  stat_smooth(aes(group = species, fill = species, col = species)) + 
  guides(color = "none") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  ylab("Average Phenological Stage") +
  xlab("DOY") +
  facet_grid(year ~ site) 

# stat smooth on Waldi  average stage ~ DOY
doy_line <- ggplot(data = subset(all_dat, site == "Waldi"), 
                   aes(x = as.factor(doy), y = average.stage, fill = species)) +
  geom_point(position = position_jitter(h = 0.1, w = 0.1), shape = 21, alpha = 0.5, size = 2) + 
  #geom_text(data = subset(avg_labels, site == "Waldi"), aes(x = Inf, y = -Inf, label = label),inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5)+
  fillSet + colSet + theme_bw() +
  stat_smooth(aes(group = species, fill = species, col = species)) + 
  guides(color = "none") +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  ylab("Average Phenological Stage") +
  xlab("DOY") +
  facet_grid(year ~ site) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", 
        strip.text = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        plot.title = element_text(hjust = 0.5))
doy_line
  
# single boxplot allenwiller
ggplot(data = subset(all_dat,site == "Allenwiller"), aes(x = as.factor(doy), y = average.stage, fill = species)) +
  #geom_jitter( aes(color = species), size = 0.8)+
  #geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.9, size = 2) + 
  #stat_smooth(aes(group=species, fill = species, col = species))+ 
  geom_boxplot(alpha = 0.8) +
  fillSet + colSet + guides(color = "none")+
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
  theme_bw() +
  ylab("Average Phenological Stage") +
  xlab("DOY") +
  facet_grid( year~ site) +
  geom_text(data = subset(avg_labels, site == "Allenwiller"), aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5)


# remove Laiselheim
all_dat <- all_dat[site != "Leiselheim"]

# compare statistically average stage 
# check distribution
hist(subset(all_dat, site == "Allenwiller")$average.stage)
hist(subset(all_dat, site == "Waldi")$average.stage)

# since there are many observations for stage 4, remove the late spring observations to avoid bias
all_dat_subset = all_dat[all_dat$doy < 125]

# boxplot all sites, separately for each year
ggplot(data = all_dat_subset, aes(x = species, y = average.stage, fill = species, col = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(color = "black", position = position_jitterdodge(jitter.width = 0.2), size = 1.5) +
  colSet + fillSet + theme_bw() + 
  ylab("Average Phenological Stage") +
  xlab("") + 
  guides(fill = "none", col = "none")+
  facet_grid(year~site)+
  stat_compare_means(label.x = 1.3, label.y = 4.2) 


## Allenwiller
all_boxplot <- ggplot(data = subset(all_dat_subset,site == "Allenwiller"), aes(x = factor(doy), y = average.stage, fill = species, col = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(color = "black", position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.5) +
  colSet + fillSet + theme_bw() +  
  facet_grid(site ~year)+
  #labs(title = "Allenwiller")+ 
  ylab("Average Phenological Stage") +
  xlab("Doy") + 
  guides(col = "none")+
  stat_compare_means(label.x = 1.3, label.y = 4.2,  aes( label =sprintf("p = %5.3f", as.numeric(..p.format..)) ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.text = element_text(size = 10),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "black", fill = "white"), 
        plot.title = element_text(hjust = 0.5)
  )

all_boxplot

## Waldi all data
wal_boxplot <- ggplot(data = subset(all_dat_subset,site == "Waldi" ), aes(x = factor(doy), y = average.stage, group = interaction(doy, species), fill = species, col = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(color = "black",position = position_jitterdodge(jitter.width = 0.2),size = 1, alpha = 0.5) +
  colSet + fillSet + theme_bw() + 
  #labs(title = "Waldi")+ 
  facet_grid(site ~year)+
  ylab("Average Phenological Stage") +
  xlab("Doy") + 
  guides(col = "none")+
  stat_compare_means(label.x = 1.3, label.y = 4.2,  aes( label =sprintf("p = %5.3f", as.numeric(..p.format..)) ))+
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right", 
        strip.text = element_text(size = 10),
        legend.key = element_rect(colour = "white"),
        legend.key.spacing.y = unit(0.4, 'cm'), 
        strip.background = element_rect(colour = "black", fill = "white"), 
        plot.title = element_text(hjust = 0.5)
  )
  
wal_boxplot


## statistical comparison of average stage
library(purrr)
library(broom)

significance_level <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ "ns"
  )
}
# all statistics for all observations
test_results <- all_dat_subset %>%
  group_by(site, year, doy) %>%
  filter(n_distinct(species) == 2) %>%  # test only if both speices are present
  nest() %>%
  mutate(
    wilcox = map(data, ~ wilcox.test(average.stage ~ species, data = .x) %>% tidy()),
    ttest = map(data, ~ t.test(average.stage ~ species, data = .x) %>% tidy()),
    mean_by_species = map(data, ~ .x %>%
                            group_by(species) %>%
                            summarise(mean_stage = mean(average.stage),
                                      n = n(),
                                      .groups = "drop"))
  ) %>%
  unnest(cols = c(wilcox, ttest), names_sep = "_") %>%
  mutate(
    sylvatica = map_chr(mean_by_species, ~ as.character(.x$species[1])),
    mean_sylvatica = map_dbl(mean_by_species, ~ .x$mean_stage[1]),
    n_sylvatica = map_dbl(mean_by_species, ~ .x$n[1]),
    hohenackeriana = map_chr(mean_by_species, ~ as.character(.x$species[2])),
    mean_hohenackeriana = map_dbl(mean_by_species, ~ .x$mean_stage[2]),
    n_hohenackeriana = map_dbl(mean_by_species, ~ .x$n[2]),
    sig_wilcox = significance_level(wilcox_p.value),
    sig_ttest = significance_level(ttest_p.value)
  ) %>%
  dplyr::select(site, year, doy,
                sylvatica, mean_sylvatica, n_sylvatica,
         hohenackeriana, mean_hohenackeriana, n_hohenackeriana,
         wilcox_statistic = wilcox_statistic,
         wilcox_p = wilcox_p.value,
         sig_wilcox,
         ttest_statistic = ttest_statistic,
         ttest_p = ttest_p.value,
         sig_ttest)
test_results

#write.csv(test_results, "Phenology/Spring_phenology_tests.csv")

#### linear interpolation to estimate DOY and GDD of stage #####
all_dat <- fread("Spring_phenology_all_sites_final.csv")
# waldi only
wal <- subset(all_dat, site == "Waldi")

stages <- c(2, 3)
output <- data.frame()

tab_for_loop <- unique(wal[, c("species", "year")])

for (i in 1:nrow(tab_for_loop)) {
  sp <- tab_for_loop$species[i]
  yr <- tab_for_loop$year[i]
  
  dat_subset <- wal %>% 
    filter(species == sp, year == yr)
  
  for (stage_target in stages) {
    
    # get DOYs just before and just after the target stage
    dat_grouped <- dat_subset %>%
      group_by(doy) %>%
      summarize(mean_stage = mean(average.stage, na.rm = TRUE), .groups = "drop")
    
    below <- dat_grouped %>%
      filter(mean_stage < stage_target) %>%
      arrange(desc(mean_stage)) %>%
      slice(1)
    
    above <- dat_grouped %>%
      filter(mean_stage > stage_target) %>%
      arrange(mean_stage) %>%
      slice(1)
    
    # skip if not enough data
    if (nrow(below) == 0 | nrow(above) == 0) next
    
    doys_to_use <- c(below$doy, above$doy)
    
    # subset data keeping onyl these dOys
    model_data <- dat_subset %>% filter(doy %in% doys_to_use)
    
    # skip if less than 3 tot obsevrations
    #if (nrow(model_data) < 3) next
    
    # fit model
    lm_doy <- lm(doy ~ average.stage, data = model_data)
    lm_gdd <- lm(Cum_GDD ~ average.stage, data = model_data)
    
    # predict stage
    new_data <- data.frame(average.stage = stage_target)
    pred_doy_ci <- predict(lm_doy, newdata = new_data, interval = "confidence")
    pred_gdd <- predict(lm_gdd, newdata = new_data)
    
    # extract
    pred_doy <- pred_doy_ci[1, "fit"]
    doy_lwr  <- pred_doy_ci[1, "lwr"]
    doy_upr  <- pred_doy_ci[1, "upr"]
    
    # model checks
    r2_doy <- summary(lm_doy)$r.squared
    pval_doy <- summary(lm_doy)$coefficients[2, 4]
    
    r2_gdd <- summary(lm_gdd)$r.squared
    pval_gdd <- summary(lm_gdd)$coefficients[2, 4]
    
    row_out <- data.frame(
      species = sp,
      year = yr,
      stage = stage_target,
      doy_before = below$doy,
      doy_after = above$doy,
      predicted_DOY = pred_doy,
      doy_lwr = doy_lwr,
      doy_upr = doy_upr,
      predicted_GDD = pred_gdd,
      r_squared_doy = r2_doy,
      p_value_doy = pval_doy,
      r_squared_gdd = r2_gdd,
      p_value_gdd = pval_gdd,
      n_obs = nrow(model_data)
    )
    
    output <- rbind(output, row_out)
  }
}
output

#fwrite(output, "Phenology/Waldi_phenology_linear_regression_results.csv")

# visual check
# calculate difference
doy_diffs <- output %>%
  filter(species %in% c("sylvatica", "orientalis")) %>%
  dplyr::select(species, year, stage, predicted_DOY) %>%
  pivot_wider(names_from = species, values_from = predicted_DOY) %>%
  mutate(doy_diff = sylvatica - orientalis,
         label = paste0("ΔDOY = ", round(doy_diff, 1)))
doy_diffs

ggplot(subset(output, year %in% c(2021, 2023, 2025)), 
       aes(x = predicted_DOY, y = stage, color = species)) +
  geom_point(size = 2) +
  geom_line() +
  ylim(c(1, 4)) +
  colSet + fillSet + theme_bw() +
  facet_wrap(~ year, ncol = 4) +
  labs(title = "Predicted DOY of Phenological Transitions", x = "DOY", y = "Stage") +
  geom_text(data = doy_diffs %>% filter(year %in% c(2021, 2023, 2025)),
            aes(x = Inf, y =0.2+stage, label = label),
            inherit.aes = FALSE, hjust = 1.3, vjust = 0.5, size = 3.5) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right", 
    strip.text = element_text(size = 10),
    legend.key = element_rect(colour = "white"),
    legend.key.spacing.y = unit(0.4, 'cm'), 
    strip.background = element_rect(colour = "black", fill = "white"), 
    plot.title = element_text(hjust = 0.5)
  )


# difference between the two subspecies
  species_diff <- output %>%
    dplyr::select(species, year, stage, predicted_DOY, predicted_GDD) %>%
    pivot_wider(
      names_from = species,
      values_from = c(predicted_DOY, predicted_GDD),
      names_sep = "_"
    ) %>%
    mutate(
      delta_DOY = predicted_DOY_sylvatica - predicted_DOY_orientalis,
      delta_GDD = predicted_GDD_sylvatica - predicted_GDD_orientalis
    ) %>%
    dplyr::select(year, stage, delta_DOY, delta_GDD)
  species_diff
 
  
  