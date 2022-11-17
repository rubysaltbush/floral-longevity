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

# error checking - births=deaths?
check <- data %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(births) == sum(deaths))

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

# add species symmetry to means
species_sym <- species_sym %>%
  dplyr::select(species, symmetry_all) %>%
  dplyr::distinct()

mean_longevity <- mean_longevity %>%
  dplyr::left_join(species_sym, by = "species")
rm(species_sym)


# boxplot of longevity by symmetry INDIVIDUALS
ggplot(data = longevity, aes(x = longevity_days, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot.pdf", width = 9, height = 5)

# t-test of longevity by symmetry
ttest <- t.test(longevity$longevity_days[longevity$symmetry_all == "zygomorphic"], 
                longevity$longevity_days[longevity$symmetry_all != "zygomorphic"])
ttest
# NO LONGER SIGNIFICANT WITH DIANELLA ADDED IN


# boxplot of longevity by symmetry
ggplot(data = mean_longevity, aes(x = mean_long, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Species mean floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot_speciesmean.pdf", width = 9, height = 5)

# t-test of species mean longevity by symmetry
ttest_species <- t.test(mean_longevity$mean_long[mean_longevity$symmetry_all == "zygomorphic"], 
                        mean_longevity$mean_long[mean_longevity$symmetry_all != "zygomorphic"])
ttest_species

# yup I was inflating my power by including all individuals rather than species :(

# out of curiosity are SD different?
ttest_speciesSD <- t.test(mean_longevity$sd_long[mean_longevity$symmetry_all == "zygomorphic"], 
                          mean_longevity$sd_long[mean_longevity$symmetry_all != "zygomorphic"])
ttest_speciesSD
# nope! they're not, ach well

# 3 zygomorphic outliers, one actinomorphic outlier, what if I removed these?

mean_longevity_nooutliers <- mean_longevity %>%
  dplyr::filter(!(species %in% c("Hakea sericea", "Grevillea buxifolia", "Hybanthus vernonii", "Acacia oxycedrus")))

# boxplot of longevity by symmetry
ggplot(data = mean_longevity_nooutliers, aes(x = mean_long, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Species mean floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot_speciesmean_nooutliers.pdf", width = 9, height = 5)

# t-test of species mean longevity by symmetry
ttest_species <- t.test(mean_longevity_nooutliers$mean_long[mean_longevity_nooutliers$symmetry_all == "zygomorphic"], 
                        mean_longevity_nooutliers$mean_long[mean_longevity_nooutliers$symmetry_all != "zygomorphic"])
ttest_species

# aha yes so it is the outliers, how justified is removing them though?
# WITH DIANELLA ADDED IN SIGNIFICANCE VERY MARGINAL
# for Hybanthus I would think that cleistogamy might be playing a role?? very low insect visitation

# seems like symmetry not the go here, possible that zygomorphic flowers EITHER
# flower for a short time bc good quality visitation OR flower for a long time
# because low frequency of visitation, and you'd need visitation data to know
# this for sure. If longevity is an evolved trait though, what determines it?

# Could I publish a paper that concludes longevity ≠ symmetry? Is this worth
# pursuing?

# Song et al (2022) suggest that longevity is more influenced by components of 
# male than female fitness - so, flowers putting more effort into pollen export
# will last longer than flowers more focussed on pollen receipt. How could you test
# this idea?

# quick test between longevity and plant height out of curiosity, not that
# I've sampled a large range of plant heights

plot(height_m ~ longevity_days, data = longevity)
plot(infloheight_m ~ longevity_days, data = longevity)

# and out of curiosity, longevity and number flowers per plant (proxy for plant attractiveness??)
plot(flowersperplant ~ longevity_days, data = longevity)
plot(budsperplant ~ longevity_days, data = longevity)
longevity$display <- longevity$flowersperplant + longevity$budsperplant
plot(display ~ longevity_days, data = longevity)
# not much there really
plot(florallength_cm ~ longevity_days, data = longevity)
plot(floraldiam_cm ~ longevity_days, data = longevity)
longevity$flordimen_cm <- as.numeric(longevity$floraldiam_cm)*as.numeric(longevity$florallength_cm)
plot(flordimen_cm ~ longevity_days, data = longevity)
# looks like some signal here for flower size? might be spurious tho

# other q's - does floral longevity = leaf longevity?
# really need to think about this all more at a later date

