################################################################################
# map_plot unit tests
################################################################################
library("phyloseq"); library("testthat"); library("ggplot2")
data("GlobalPatterns")
# Subset to small dataset for quicker testing
GP <- prune_taxa(tail(names(sort(taxa_sums(GlobalPatterns))), 50), GlobalPatterns)

# Pretend GP doesn't have sample_data or tax_table
GP.tax <- tax_table(GP)
GP.sd  <- sample_data(GP)
GP.tr  <- phy_tree(GP)
# GP <- phyloseq(otu_table(GP), GP.tr)
GP.otu <- otu_table(GP)

# Try ordination
GP.ord <- ordinate(GP.otu, "DCA")
