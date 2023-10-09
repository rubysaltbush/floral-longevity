# script to compare longevity between zygomorphic and actinomorphic taxa
# just using straight up means, no phylogeny at first

# boxplot of longevity by symmetry for all individuals
ggplot(data = sym_long, aes(x = mean_long_days, y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_discrete(type = my_colours$symmetry) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot_individuals.pdf", width = 9, height = 5)

# check distribution of longevity for actin vs zyg taxa
hist(sym_long$mean_long_days[sym_long$sym_species == "zygomorphic"])
hist(sym_long$mean_long_days[sym_long$sym_species == "actinomorphic"])
# look very similar

# t-test of longevity by symmetry
t.test(sym_long$mean_long_days[sym_long$sym_species == "zygomorphic"], 
       sym_long$mean_long_days[sym_long$sym_species == "actinomorphic"])

# p = 0.0001, actinomorphic mean 3.85, zygomorphic 4.67
# zygomorphic flowers longer lived which makes sense if fewer visitors
# BUT shall have to see if this difference remains when phylogeny considered
# and genera subsampled to reduce taxonomic bias

# what are the longest lived flowers?
sym_long %>%
  dplyr::arrange(dplyr::desc(mean_long_days)) %>%
  dplyr::select(Accepted_name, Accepted_family, sym_species, mean_long_days) %>%
  dplyr::distinct() %>%
  head()
# two orchids, a Lentibulariaceae, Siparuna muricata and and Ericaceae

# what are the shortest lived flowers?
sym_long %>%
  dplyr::arrange(mean_long_days) %>%
  dplyr::select(Accepted_name, Accepted_family, sym_species, mean_long_days) %>%
  dplyr::distinct() %>%
  head()
# a real mix of families, some with pseudanthia - Asteraceae, Araceae
# and some with solo flowers - an iris, Talinaceae, Orobanchaceae, Cistaceae

