# script to read in floral longevity field data and process raw data
# into mean and standard deviation floral longevity per species
# coupled with floral symmetry

fieldlong <- cache_csv("data_output/mean_longevity_Sydney_fieldwork.csv", 
                       function() {

# read in longevity monitoring data straight from google sheet
data <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQlJKKlDeLMZYH32YtGRqpnL9kJFLfiUOokQZ51kvvRTgvvx4WwpemWiwCnz6hlMarYmRViWOQxVbHn/pub?gid=1585669731&single=true&output=csv")

# convert timestamp column to R date/time object
data$Timestamp <- strptime(data$Timestamp, "%d/%m/%Y %H:%M:%S")

# seems to automatically detect daylight savings/non-daylight savings timezones
head(data$Timestamp)

# if you wanted to convert directly to number of seconds since 1970
# as.numeric(data$Timestamp[1,])
# and to convert back to POSIXct
# as.POSIXct(as.numeric(data$Timestamp[1,]), origin = "1970-01-01")

colnames(data) <- c("timestamp", "site", "species", "individual", "observer", 
                    "buds", "flowers", "finished", "notes", "pop_flowering", 
                    "plant_flowering", "no_birth", "no_death",
                    "keep", "births", "deaths", "notes2")

# remove unwanted rows from monitoring data
# these have been individually marked, unwanted rows include flowers that got
# eaten before monitoring finished or that stayed in bud stage for a long time etc
data <- data %>%
  dplyr::filter(keep == 1)

# error checking - births=deaths?
check <- data %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(sum(births) == sum(deaths))
print(paste(sum(check$`sum(births) == sum(deaths)`), "individuals of", nrow(check), "individuals have equal births and deaths"))
rm(check)

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
# 4 individuals with only 1 flower monitored, 9 with 2 flowers
paste("on average", round(mean(longevity$no_flowers), 1), "flowers monitored per plant")
# "on average 11.3 flowers monitored per plant"
paste("standard error", 
      round((sd(longevity$no_flowers)/sqrt(length(longevity$no_flowers))), 1))
# "standard error 0.7"

# calculate mean longevity per plant in days
longevity$longevity_days <- difftime(longevity$deathtime, longevity$birthtime, units = "days")

# get trait data measured for each individual
trait_data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
colnames(trait_data) <- c("timestamp", "site", "latitude", "longitude", 
                          "observer", "species", "individual", 
                          "infloheight_m", "habit", "height_m", "inflosize_cm", 
                          "flowersperplant", "budsperplant", "florallength_cm", 
                          "floraldiam_cm", "colour", "symmetry", "tube", "notes",
                          "no_carpels",	"no_ovules",	"no_stamens")

longevity <- longevity %>%
  dplyr::left_join(trait_data, by = "individual") %>%
  dplyr::select(individual:longevity_days, infloheight_m:colour, tube, 
                no_carpels:no_stamens)

# symmetry not recorded in trait_data for all individuals, get it by species
species_sym_and_site <- trait_data %>%
  dplyr::select(species, site, latitude, longitude, symmetry) %>%
  dplyr::distinct()
species_individual <- data %>%
  dplyr::select(species, individual) %>%
  dplyr::distinct()
species_sym_and_site <- species_individual %>%
  dplyr::left_join(species_sym_and_site, by = "species")
rm(species_individual, trait_data)

# join symmetry to longevity
longevity <- longevity %>%
  dplyr::left_join(species_sym_and_site, by = "individual")

# remove specifics from actinomorphic symmetry category
longevity$symmetry <- gsub(" \\(incl. rotational and spiral\\)", "", longevity$symmetry)

# boxplot of longevity by symmetry INDIVIDUALS
ggplot(data = longevity, aes(x = longevity_days, y = symmetry, fill = symmetry)) +
  geom_boxplot() +
  scale_fill_discrete(type = my_colours$symmetry) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/field_data_individuals_symmetry_longevity_boxplot.pdf", width = 9, height = 5)

# calculate mean longevity per species
fieldlong <- longevity %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(mean_long_days = mean(longevity_days, na.rm = TRUE))

# calculate standard deviation per species
sd <- longevity %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(sd_long = sd(longevity_days, na.rm = TRUE))

# calculate n per species
n <- longevity %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(n = n())

# join all together
fieldlong <- fieldlong %>%
  dplyr::left_join(sd, by = "species") %>%
  dplyr::left_join(n, by = "species")
rm(sd, n, longevity)

# add species symmetry to means
species_sym_and_site <- species_sym_and_site %>%
  dplyr::select(species, site, latitude, longitude, symmetry) %>%
  dplyr::distinct()

fieldlong <- fieldlong %>%
  dplyr::left_join(species_sym_and_site, by = "species")
rm(species_sym_and_site, data)

# boxplot of longevity by symmetry SPECIES
ggplot(data = fieldlong, aes(x = mean_long_days, y = symmetry, fill = symmetry)) +
  geom_boxplot() +
  scale_y_discrete(labels = c("actinomorphic", "zygomorphic")) +
  scale_fill_viridis_d(alpha = 0.6, direction = -1) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Species mean floral longevity (days)") +
  ylab("")
ggsave("figures/field_data_speciesmean_symmetry_longevity_boxplot.pdf", width = 7, height = 5)

# t-test of species mean longevity by symmetry
t.test(fieldlong$mean_long_days[fieldlong$symmetry == "zygomorphic"], 
       fieldlong$mean_long_days[fieldlong$symmetry != "zygomorphic"])
# Welch Two Sample t-test
# 
# data:  fieldlong$mean_long_days[fieldlong$symmetry == "zygomorphic"] and fieldlong$mean_long_days[fieldlong$symmetry != "zygomorphic"]
# t = -0.51429 days, df = 31.946, p-value = 0.6106 days
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.457111 days  2.063379 days
# sample estimates:
#   Time differences in days
# mean of x mean of y 
# 5.737704  6.434570 
# no difference in mean longevity by symmetry for these 34 species

# output species mean longevity with symmetry
readr::write_csv(fieldlong, "data_output/mean_longevity_Sydney_fieldwork.csv")

})

# possible that zygomorphic flowers EITHER flower for a short time bc good 
# quality visitation OR flower for a long time because low frequency of 
# visitation, and you'd need visitation data to tease this out in field study

# given phylogenetic signal in both floral symmetry and floral longevity, 
# likely need a much larger sample of multiple families and genera to explore
# evolutionary relationship between these two traits

