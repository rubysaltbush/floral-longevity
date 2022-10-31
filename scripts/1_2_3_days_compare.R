# WHAT IF I ONLY CHECKED SASSAFRAS GULLY FLOWERS EVERY 2 DAYS???
# does it make significant difference to longevity estimates??
# read in data to check this idea

data <- readr::read_csv("data_input/sassafras_longevity_if2days.csv")

# convert timestamp column to R date/time object
data$timestamp <- strptime(data$timestamp, "%d/%m/%Y %H:%M")

# split into 1, 2 or 3 days monitoring data
data1 <- data %>%
  dplyr::filter(keep1 == 1) %>%
  dplyr::select(timestamp:deaths1)

# remove unwanted rows from monitoring data
data2 <- data %>%
  dplyr::filter(keep2 == 1) %>%
  dplyr::select(timestamp:notes, keep2:deaths2)

data3 <- data %>%
  dplyr::filter(keep3 == 1) %>%
  dplyr::select(timestamp:notes, keep3:deaths3)

# calculate weighted mean date of "births" i.e. bud -> flower
# weight is number of flowers per inflorescence monitored
longevity1 <- data1 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, births1))

longevity2 <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, births2))

longevity3 <- data3 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, births3))

# calculate weighted mean date of "deaths" i.e. flower -> finished flower
deathtime1 <- data1 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, deaths1))

deathtime2 <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, deaths2))

deathtime3 <- data3 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(timestamp, deaths3))

# calculate number of flowers monitored per plant
no_flowers1 <- data1 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(deaths1))

no_flowers2 <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(deaths2))

no_flowers3 <- data3 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(deaths3))

# join all together
longevity1 <- longevity1 %>%
  dplyr::left_join(deathtime1, by = "individual") %>%
  dplyr::left_join(no_flowers1, by = "individual")
colnames(longevity1) <- c("individual", "birthtime1", "deathtime1", "no_flowers1")
rm(deathtime1, no_flowers1, data1)

longevity2 <- longevity2 %>%
  dplyr::left_join(deathtime2, by = "individual") %>%
  dplyr::left_join(no_flowers2, by = "individual")
colnames(longevity2) <- c("individual", "birthtime2", "deathtime2", "no_flowers2")
rm(deathtime2, no_flowers2, data2)

longevity3 <- longevity3 %>%
  dplyr::left_join(deathtime3, by = "individual") %>%
  dplyr::left_join(no_flowers3, by = "individual")
colnames(longevity3) <- c("individual", "birthtime3", "deathtime3", "no_flowers3")
rm(deathtime3, no_flowers3, data3)

# then calculate longevity in days as the difference between the weighted mean
# birth time and the weighted mean death time of all flowers
longevity1$longevity_days1 <- difftime(longevity1$deathtime1, longevity1$birthtime1, units = "days")

longevity2$longevity_days2 <- difftime(longevity2$deathtime2, longevity2$birthtime2, units = "days")

longevity3$longevity_days3 <- difftime(longevity3$deathtime3, longevity3$birthtime3, units = "days")

# join individual longevities together to compare for individuals
longevity_comp_ind <- longevity1 %>%
  dplyr::left_join(longevity2, by = "individual") %>%
  dplyr::left_join(longevity3, by = "individual")
rm(longevity1, longevity2, longevity3)

# get species and symmetry for individuals
# get trait data measured for each individual
trait_data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
colnames(trait_data) <- c("timestamp", "site", "observer", "species", "individual", 
                          "infloheight_m", "habit", "height_m", "inflosize_cm", 
                          "flowersperplant", "budsperplant", "florallength_cm", 
                          "floraldiam_cm", "colour", "symmetry", "tube", "notes")

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
rm(species_individual, trait_data)

# join symmetry to longevity
longevity_comp_ind <- longevity_comp_ind %>%
  dplyr::left_join(species_sym, by = "individual")

# calculate mean and standard deviation per species
longevity_comp_mean <- longevity_comp_ind %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean_long1 = mean(longevity_days1, na.rm = TRUE),
                   sd_long1 = sd(longevity_days1, na.rm = TRUE),
                   mean_long2 = mean(longevity_days2, na.rm = TRUE),
                   sd_long2 = sd(longevity_days2, na.rm = TRUE),
                   mean_long3 = mean(longevity_days3, na.rm = TRUE),
                   sd_long3 = sd(longevity_days3, na.rm = TRUE))

# looks like it makes SOME difference between 1 to 2 days, bigger difference
# for 3 days


# calculate differences in longevity
# for species means
longevity_comp_mean$diff_mean_1_2 <- longevity_comp_mean$mean_long1 - longevity_comp_mean$mean_long2
longevity_comp_mean$diff_sd_1_2 <- longevity_comp_mean$sd_long1 - longevity_comp_mean$sd_long2  
longevity_comp_mean$diff_mean_1_3 <- longevity_comp_mean$mean_long1 - longevity_comp_mean$mean_long3
longevity_comp_mean$diff_sd_1_3 <- longevity_comp_mean$sd_long1 - longevity_comp_mean$sd_long3  
# and individuals
longevity_comp_ind$diff_long_1_2 <- longevity_comp_ind$longevity_days1 - longevity_comp_ind$longevity_days2
longevity_comp_ind$diff_long_1_3 <- longevity_comp_ind$longevity_days1 - longevity_comp_ind$longevity_days3

# very small difference between 1 and 2 days, but quite a large difference for 3 days :(

#write output to discuss with H and R

write_csv(longevity_comp_mean, "data_output/compare_longevity_1or2daymonitoring_means.csv")
write_csv(longevity_comp_ind, "data_output/compare_longevity_1or2daymonitoring.csv")
