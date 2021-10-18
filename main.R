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
                    "plant_flowering", "keep", "births", "deaths", "notes2")

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
species_sym <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
species_sym <- species_sym %>%
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

# boxplot of longevity by symmetry
ggplot(data = longevity, aes(x = longevity_days, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot.png", width = 9, height = 5)

# t-test of longevity by symmetry
ttest <- t.test(longevity$longevity_days[longevity$symmetry == "zygomorphic"], 
                longevity$longevity_days[longevity$symmetry != "zygomorphic"])
ttest

