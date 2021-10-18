library(tidyverse)
library(ggpubr)

data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQlJKKlDeLMZYH32YtGRqpnL9kJFLfiUOokQZ51kvvRTgvvx4WwpemWiwCnz6hlMarYmRViWOQxVbHn/pub?gid=1585669731&single=true&output=csv")

data$Timestamp <- strptime(data$Timestamp, "%d/%m/%Y %H:%M:%S")

data$Timestamp[1,]

as.numeric(data$Timestamp[1,])

as.POSIXct(as.numeric(data$Timestamp[1,]), origin = "1970-01-01")

#make new time column with date as seconds since 1970
data$time_sec <- as.numeric(data$Timestamp)

data2 <- data
colnames(data2) <- c("timestamp", "site", "species", "individual", "observer", "buds", "flowers", "finished", "notes", "pop_flowering", "plant_flowering", "keep", "births", "deaths", "notes2", "time_sec")

data2 <- data2 %>%
  dplyr::filter(keep == 1)

longevity <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(time_sec, births))

deathtime <- data2 %>%
  dplyr::group_by(individual) %>%
  dplyr::summarise(weighted.mean(time_sec, deaths))

longevity <- longevity %>%
  dplyr::left_join(deathtime, by = "individual")
colnames(longevity) <- c("individual", "birthtime", "deathtime")
rm(deathtime)

longevity$longevity <- longevity$deathtime - longevity$birthtime
longevity$longevity_days <- longevity$longevity/(60*60*24)

trait_data <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
trait_data <- trait_data %>%
  dplyr::select(5:17)
colnames(trait_data) <- c("individual", "infloheight_m", "habit", "height_m", "inflosize_cm", "flowersperplant", "budsperplant", "florallength_cm", "floraldiam_cm", "colour", "symmetry", "tube", "notes")

longevity <- longevity %>%
  dplyr::left_join(trait_data, by = "individual")

species_sym <- read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vQ8l1sSP14c_ofDymqna9mTqeE6KK1scNGt6YBTCMGhvqeh884eyW1JwNQVvdL-znAuXGxjcOh-sA-t/pub?gid=1824852273&single=true&output=csv")
species_sym <- species_sym %>%
  dplyr::select(4,15) %>%
  dplyr::distinct()
colnames(species_sym) <- c("species", "symmetry_all")
data3 <- data2 %>%
  dplyr::select(species, individual) %>%
  dplyr::distinct()
species_sym <- data3 %>%
  dplyr::left_join(species_sym, by = "species")
rm(data3)

longevity <- longevity %>%
  dplyr::left_join(species_sym, by = "individual")

ggplot(data = longevity, aes(x = longevity_days, y = symmetry_all, fill = symmetry_all)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot.png", width = 9, height = 5)

ttest <- t.test(longevity$longevity_days[longevity$symmetry_all == "zygomorphic"], 
                longevity$longevity_days[longevity$symmetry_all != "zygomorphic"])

