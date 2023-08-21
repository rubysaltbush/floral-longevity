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

# hmm okay, definitely quite a range! shall see how this changes as I score symmetry

# t-test of longevity by symmetry
ttest <- t.test(sym_long$mean_long_days[sym_long$sym_species == "zygomorphic"], 
                sym_long$mean_long_days[sym_long$sym_species != "zygomorphic"])
ttest
# p = 0.06, with actinomorphic mean 3.72 and zygomorphic 4.3
# shall see how this changes as I add in taxa, and consider phylogeny

# now p = 0.007, actinomorphic mean 3.84, zygomorphic 4.48
rm(ttest)
