# script to compare longevity between zygomorphic and actinomorphic taxa
# just using straight up means, no phylogeny at first

# boxplot of longevity by symmetry INDIVIDUALS
ggplot(data = sym_long, aes(x = mean_long_days, y = sym_species, fill = sym_species)) +
  geom_boxplot() +
  scale_fill_viridis_d(alpha = 0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggpubr::theme_pubr(legend = "none") +
  xlab("Floral longevity (days)") +
  ylab("")
ggsave("figures/symmetry_longevity_boxplot.pdf", width = 9, height = 5)

# for Hervé, check distribution of longevity for actin vs zyg taxa
hist(sym_long$mean_long_days[sym_long$sym_species == "zygomorphic"])
hist(sym_long$mean_long_days[sym_long$sym_species == "actinomorphic"])

# t-test of longevity by symmetry
ttest <- t.test(sym_long$mean_long_days[sym_long$sym_species == "zygomorphic"], 
                sym_long$mean_long_days[sym_long$sym_species == "actinomorphic"])
ttest
# p = 0.0001, actinomorphic mean 3.85, zygomorphic 4.67
# zygomorphic flowers longer lived which makes sense if fewer visitors
# BUT shall have to see if this difference remains when phylogeny considered
# and genera subsampled to reduce taxonomic bias
rm(ttest)
