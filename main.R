library(tidyverse)
library(ggpubr)

# read in longevity monitoring data straight from google sheet
data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQlJKKlDeLMZYH32YtGRqpnL9kJFLfiUOokQZ51kvvRTgvvx4WwpemWiwCnz6hlMarYmRViWOQxVbHn/pub?gid=1585669731&single=true&output=csv")

# convert timestamp column to R date/time object
data$Timestamp <- strptime(data$Timestamp, "%d/%m/%Y %H:%M:%S")

# seems to automatically detect daylight savings/non-daylight savings timezones
data$Timestamp

# if you wanted to convert directly to number of seconds since 1970
# as.numeric(data$Timestamp[1,])
# and to convert back to POSIXct
# as.POSIXct(as.numeric(data$Timestamp[1,]), origin = "1970-01-01")

colnames(data) <- c("timestamp", "site", "species", "individual", "observer", 
                    "buds", "flowers", "finished", "notes", "pop_flowering", 
                    "plant_flowering", "no_birth", "no_death",
                    "keep", "births", "deaths", "notes2")

# remove unwanted rows from monitoring data
data <- data %>%
  dplyr::filter(keep == 1)

# calculate weighted mean date of "births" i.e. bud -> flower
# weight is number of flowers per inflorescence monitored
longevity <- data %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, births))

# calculate weighted mean date of "deaths" i.e. flower -> finished flower
deathtime <- data %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, deaths))

# calculate number of flowers monitored per plant
no_flowers <- data %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(deaths))

longevity <- longevity %>%
  dplyr::left_join(deathtime, by = "individual") %>%
  dplyr::left_join(no_flowers, by = "individual")
colnames(longevity) <- c("individual", "birthtime", "deathtime", "no_flowers")
rm(deathtime, no_flowers)

longevity$longevity_days <- difftime(longevity$deathtime, longevity$birthtime, units = "days")

# get trait data measured for each individual
trait_data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
colnames(trait_data) <- c("timestamp", "site", "observer", "species", "individual", 
                          "infloheight_m", "habit", "height_m", "inflosize_cm", 
                          "flowersperplant", "budsperplant", "florallength_cm", 
                          "floraldiam_cm", "colour", "symmetry", "tube", "notes")

longevity <- longevity %>%
  dplyr::left_join(trait_data, by = "individual") %>%
  dplyr::select(individual:longevity_days, infloheight_m:tube)

# symmetry not recorded in trait_data for all individuals, get it by species

species_sym <- trait_data %>%
  dplyr::select(4,15) %>%
  dplyr::distinct()
colnames(species_sym) <- c("species", "symmetry_all")
species_individual <- data %>%
  dplyr::select(species, individual) %>%
  dplyr::distinct()
species_sym <- species_individual %>%
  dplyr::left_join(species_sym, by = "species")
rm(species_individual)

# join symmetry to longevity
longevity <- longevity %>%
  dplyr::left_join(species_sym, by = "individual")

# save this data as csv to check later
readr::write_csv(longevity, "data_output/floral_longevity_output.csv")

# calculate mean per species
mean_longevity <- longevity %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean_long = mean(longevity_days, na.rm = TRUE))

# calculate standard deviation per species
sd <- longevity %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(sd_long = sd(longevity_days, na.rm = TRUE))

mean_longevity <- mean_longevity %>%
  dplyr::left_join(sd, by = "species")
rm(sd)


# boxplot of longevity by symmetry
ggplot(data = longevity, aes(x = longevity_days, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot.pdf", width = 9, height = 5)

# t-test of longevity by symmetry
ttest <- t.test(longevity$longevity_days[longevity$symmetry == "zygomorphic"], 
                longevity$longevity_days[longevity$symmetry != "zygomorphic"])
ttest

# WHAT IF I ONLY CHECKED SASSAFRAS GULLY FLOWERS EVERY 2 DAYS???
# does it make significant difference to longevity estimates??
# read in data to check this idea

data2 <- readr::read_csv("data_input/sassafras_longevity_if2days.csv")
 
# convert timestamp column to R date/time object
data2$timestamp <- strptime(data2$timestamp, "%d/%m/%Y %H:%M")

# remove unwanted rows from monitoring data
data2 <- data2 %>%
  dplyr::filter(keep2 == 1)

# calculate weighted mean date of "births" i.e. bud -> flower
# weight is number of flowers per inflorescence monitored
longevity2 <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, births2))

# calculate weighted mean date of "deaths" i.e. flower -> finished flower
deathtime <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, deaths2))

# calculate number of flowers monitored per plant
no_flowers <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(deaths2))

longevity2 <- longevity2 %>%
  dplyr::left_join(deathtime, by = "individual") %>%
  dplyr::left_join(no_flowers, by = "individual")
colnames(longevity2) <- c("individual", "birthtime2", "deathtime2", "no_flowers2")
rm(deathtime, no_flowers)

longevity2$longevity_days2 <- difftime(longevity2$deathtime2, longevity2$birthtime2, units = "days")

# join symmetry to longevity
longevity2 <- longevity2 %>%
  dplyr::left_join(species_sym, by = "individual")

# calculate mean per species
mean_longevity2 <- longevity2 %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean_long2 = mean(longevity_days2, na.rm = TRUE))

# calculate standard deviation per species
sd2 <- longevity2 %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(sd_long2 = sd(longevity_days2, na.rm = TRUE))

mean_longevity2 <- mean_longevity2 %>%
  dplyr::left_join(sd2, by = "species")
rm(sd2)

longevity_comp_mean <- mean_longevity2 %>%
  dplyr::left_join(mean_longevity, by = "species")

# looks like it makes SOME difference, but not huge?

longevity_comp_ind <- longevity %>%
  dplyr::select(1:5) %>%
  dplyr::right_join(longevity2, by = "individual")

# calculate differences in longevity

longevity_comp_mean$diff_mean <- longevity_comp_mean$mean_long - longevity_comp_mean$mean_long2
longevity_comp_mean$diff_sd <- longevity_comp_mean$sd_long - longevity_comp_mean$sd_long2  
longevity_comp_ind$diff_long <- longevity_comp_ind$longevity_days - longevity_comp_ind$longevity_days2

#write output to discuss with H and R

write_csv(longevity_comp_mean, "data_output/compare_longevity_1or2daymonitoring_means.csv")
write_csv(longevity_comp_ind, "data_output/compare_longevity_1or2daymonitoring.csv")
